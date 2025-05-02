#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

use strict;
use warnings;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Slice;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::Utils::AncestralAllelesUtils;
use Bio::EnsEMBL::Variation::VariationFeature;

use Compress::Zlib;
use FileHandle;
use Getopt::Long;
use List::Util qw(first);
use Pod::Usage qw(pod2usage);
$| = 1;

my $args = scalar @ARGV;
my $config = {};
my $gvf_line;
GetOptions(
    $config,
    'help|h',

    'gvf_file=s',
    'vcf_file=s',
    'species=s',
    'registry=s',
    'fasta_file=s',
    'ancestral_allele|aa',
    'ancestral_allele_file=s',
    'global_maf',
    'evidence',
    'clinical_significance',
    'structural_variations|svs',
    'incl_consequences',
    'protein_coding_details',
    'sift',
    'polyphen',
    'somatic',
    'set_name=s',
    'individual=s',

) or pod2usage(1);

pod2usage(1) if ($config->{'help'} || !$args);

main($config);

sub main {
  my $config = shift;
   
  if ($config->{'fasta_file'}) {
    my $fai_index = Bio::DB::HTS::Faidx->new($config->{'fasta_file'});
    die("Unable to get FASTA index") if (!$fai_index);
    $config->{fasta_index} = $fai_index;
  } 

    init_db_connections($config);
    init_data($config);
    if( -e $config->{vcf_file} && !-z $config->{vcf_file}) {
      unlink($config->{vcf_file}) or die $config->{vcf_file} . " already exists and can not delete it\n";
    }
    my $fh_vcf = FileHandle->new('> ' . $config->{vcf_file});
    $config->{fh} = $fh_vcf;
    print_header($config);
    read_gvf_file($config);
    my $fh = $config->{fh};
    $fh->close();
}

sub init_db_connections {
    my $config = shift;
    my $registry = 'Bio::EnsEMBL::Registry';
    my $registry_file = $config->{registry};
    die "Could not find registry_file $registry_file" unless (-e $registry_file);
    $registry->load_all($registry_file);
    my $ontology = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
    my $vdba = $registry->get_DBAdaptor($config->{species}, 'variation');
    my $cdba = $registry->get_DBAdaptor($config->{species}, 'core');
    $config->{vdba} = $vdba;
    $config->{cdba} = $cdba;
    $config->{slice_adaptor} = $cdba->get_SliceAdaptor;
    $config->{ontology_adaptor} = $ontology if ($ontology); 
}

sub init_data {
    my $config = shift;

    my $vdba = $config->{vdba};
    my $dbh = $vdba->dbc->db_handle;

    my $source_adaptor = $vdba->get_SourceAdaptor;   
    my $sources = $source_adaptor->fetch_all;
   
    my $data_type = 'variation';
    if ($config->{structural_variations}) {
        $data_type = 'structural_variation';
    }
    my $source_to_desc = {};
    # Get variant source and source description for VCF header
    # Try getting sources by data_types from source table first 
    foreach my $source (@$sources) {
      my @source_data_types = @{$source->get_all_data_types};
      if (grep {$_ eq $data_type} @source_data_types) {
        my $name = $source->name;
        my $desc = $source->description;
        my $version = $source->version;
        $name = ($version ? "$name\_$version" : $name);
        $desc =~ s|<.+?>||g;
        $source_to_desc->{$name} = $desc;
      }
    }
    # If the data_type hasen't been set in the source table
    # extract all sources from the variation or structural_variation table 
    if (scalar keys %$source_to_desc == 0) {
      my $sth = $dbh->prepare(qq{
          select distinct s.name, s.version, s.description from source s, $data_type v
          where s.source_id = v.source_id;
      }, {mysql_use_result => 1});
      my $source_to_desc = {};
      my ($source_name, $source_version, $source_desc);
      $sth->execute();
      $sth->bind_columns(\($source_name, $source_version, $source_desc));
      while ($sth->fetch) {
          my $source = ($source_version ? "$source_name\_$source_version" : $source_name);
          $source_desc =~ s|<.+?>||g;
          $source_to_desc->{$source} = $source_desc;
      }
      $sth->finish(); 
    }
    $config->{source_to_desc} = $source_to_desc;

    $config->{abbr_to_evidence_value} = {
        'E_Multiple_observations' => 'Multiple_observations',
        'E_Freq' => 'Frequency',
        'E_Hapmap' => 'HapMap',
        'E_1000G' => '1000Genomes',
        'E_Cited' => 'Cited',
        'E_ESP' => 'ESP',
        'E_Phenotype_or_Disease' => 'Phenotype_or_Disease', 
        'E_TOPMed' => 'TOPMed',
        'E_gnomAD' => 'gnomAD',
        'E_ExAC' => 'ExAC',
    };
    $config->{evidence_value_to_abbr} = {
        'Multiple_observations' => 'E_Multiple_observations',
        'Frequency' => 'E_Freq',
        'HapMap' => 'E_Hapmap',
        '1000Genomes' => 'E_1000G',
        'Cited' => 'E_Cited',
        'ESP' => 'E_ESP',
        'Phenotype_or_Disease' => 'E_Phenotype_or_Disease',
        'TOPMed' => 'E_TOPMed',
        'gnomAD' => 'E_gnomAD',
        'ExAC' => 'E_ExAC',
    };
    $config->{clin_significance_to_abbr} = {
        'uncertain significance' => 'CLIN_uncertain_significance',
        'not provided' => 'CLIN_not_provided',
        'benign' => 'CLIN_benign',
        'likely benign' => 'CLIN_likely_benign',
        'likely pathogenic' => 'CLIN_likely_pathogenic',
        'pathogenic' => 'CLIN_pathogenic',
        'drug response' => 'CLIN_drug_response',
        'histocompatibility' => 'CLIN_histocompatibility',
        'other' => 'CLIN_other',
        'confers sensitivity' => 'CLIN_confers_sensitivity',
        'risk factor' => 'CLIN_risk_factor',
        'association' => 'CLIN_association',
        'protective' => 'CLIN_protective',
    };
    $config->{abbr_to_clin_significance} = {
        'CLIN_uncertain_significance' => 'uncertain significance',
        'CLIN_not_provided' => 'not provided',
        'CLIN_benign' => 'benign',
        'CLIN_likely_benign' => 'likely benign',
        'CLIN_likely_pathogenic' => 'likely pathogenic',
        'CLIN_pathogenic' => 'pathogenic',
        'CLIN_drug_response' => 'drug response',
        'CLIN_histocompatibility' => 'histocompatibility',
        'CLIN_other' => 'other',
        'CLIN_confers_sensitivity' => 'confers sensitivity',
        'CLIN_risk_factor' => 'risk factor',
        'CLIN_association' => 'association',
        'CLIN_protective' => 'protective',
    };
    $config->{svs_gvf2vcf} = {};
    if ($config->{structural_variations}) {
      # Get all available structural variant classes and descriptions to be stored in the VCF header
      my $sv_var_class_names = {};
      my $sth = $dbh->prepare(qq/select distinct a.value from attrib a, structural_variation_feature svf where svf.class_attrib_id = a.attrib_id;/);
      $sth->execute() or die $sth->errstr;
      while (my $row = $sth->fetchrow_arrayref) {
        $sv_var_class_names->{$row->[0]} = 1;
      }
      $sth->finish;
      my $ontology = $config->{ontology_adaptor};
      foreach my $name (keys %$sv_var_class_names) {
        $config->{svs_gvf2vcf}->{$name} = $name;
        my $definition = $name;
        if ($ontology) {
          my $terms = $ontology->fetch_all_by_name($name, 'SO');
          foreach my $term (@$terms) {
            $definition = $term->definition;
            # trim definition remove everything between [], remove " and \ characters
            $definition =~ s/\[.*\]|\\|"//g;
          }
        } 
        $config->{header_sv_class}->{$name} = $definition;
      }
    }

    my @vep_consequence_info =  qw/Allele Consequence Feature_type Feature/;
    if ($config->{protein_coding_details}) {
      push @vep_consequence_info, 'Amino_acids';
    }
    if ($config->{sift}) {
      push @vep_consequence_info, 'SIFT';
    }
    if ($config->{polyhen}) {
      push @vep_consequence_info, 'PolyPhen';
    }

    $config->{vep} = \@vep_consequence_info;

    if ($config->{ancestral_allele_file}) {
      my $ancestral_fai_index = Bio::DB::HTS::Faidx->new($config->{ancestral_allele_file});
      my $ancestral_allele_utils = Bio::EnsEMBL::Variation::Utils::AncestralAllelesUtils->new(-fasta_db => $ancestral_fai_index);
      die("Unable to get ancestral FASTA index") if (!$ancestral_fai_index);         
      $config->{ancestral_allele_utils} = $ancestral_allele_utils;  
    }
}

sub read_gvf_file {
  my $config = shift;
  my $gvf_file = $config->{gvf_file};

  if ($gvf_file =~ /gz$/) {
    my $fh_gvf = gzopen($gvf_file, "rb") or die "Error reading $gvf_file: $gzerrno\n";
    while ($fh_gvf->gzreadline($_) > 0) {
        chomp;
        my $line = $_;
        next if ($line =~ /^##/);
        parse_gvf_line($config, $line);
    }
    die "Error reading $gvf_file: $gzerrno\n" if $gzerrno != Z_STREAM_END;
    $fh_gvf->gzclose(); 
  } else {
    my $fh_gvf = FileHandle->new($gvf_file, 'r') or die "Error reading $gvf_file $!\n";
    while (<$fh_gvf>) {
        chomp;
        my $line = $_;
        next if ($line =~ /^##/);
        parse_gvf_line($config, $line);
    }
    $fh_gvf->close(); 
  }
}

sub parse_gvf_line {
    my $config = shift;
    my $line = shift;

    my $gvf_line = get_gvf_line(\$line);
    my $vcf_line = {};

    $vcf_line->{'QUAL'} = '.';
    $vcf_line->{'FILTER'} = '.';
    my @dbxref = split(':', $gvf_line->{Dbxref}, 2);
    my ($variation_id, $db) = ($dbxref[1], $dbxref[0]);
    $vcf_line->{'ID'} = $variation_id;
    push @{$vcf_line->{INFO}}, $db;

    if ($config->{structural_variations}) {
        add_svs_annotation($config, $vcf_line, $gvf_line);
        print_vcf_line($config, $vcf_line);
        return;
    }

    add_position_and_alleles($config, $vcf_line, $gvf_line);

    my $type_of_sequence_alteration = $gvf_line->{type};
    push @{$vcf_line->{INFO}}, "TSA=$type_of_sequence_alteration";

    if ($gvf_line->{evidence_values}) {
        push @{$vcf_line->{INFO}}, map { $config->{evidence_value_to_abbr}->{$_} } split(',', $gvf_line->{evidence_values});
    }

    if ($gvf_line->{clinical_significance}) {
        push @{$vcf_line->{INFO}}, map { $config->{clin_significance_to_abbr}->{$_} } split(',', $gvf_line->{clinical_significance});
    }

    if ($gvf_line->{global_minor_allele_frequency}) {
        add_gmaf($vcf_line, $gvf_line);
    }

    if ($gvf_line->{Genotype}) {
        add_genotypes($vcf_line, $gvf_line);
    }

    if ($config->{incl_consequences} && (!defined $gvf_line->{Variant_effect}) ) {
        push @{$vcf_line->{INFO}}, 'VE=intergenic_variant';
    }

    my $attributes = {
      'Variant_effect' => 'VE',
      'variant_peptide' => 'VarPep',
      'sift_prediction' => 'Sift',
      'polyphen_prediction' => 'Polyphen',
      'reference_peptide' => 'RefPep',
      'Variant_freq' => 'AF',
      'ancestral_allele' => 'AA',
     };

    while (my ($attribute, $key) = each %$attributes) {
        if (defined $gvf_line->{$attribute}) {
            my $value = $gvf_line->{$attribute};
            add_info($vcf_line, $key, $value);
        }
    }

    if ($config->{incl_consequences}) {
      parse_consequence_info($gvf_line, $vcf_line); 
    }
    
    print_vcf_line($config, $vcf_line);

}

sub get_gvf_line {
    my $line = shift;
    my $gvf_line = {};
    my ($seq_id, $source, $type, $start, $end, $score, $strand, $phase, $attrib) = split(/\t/, $$line);
    $gvf_line->{seq_id} = $seq_id;
    $gvf_line->{source} = $source;
    $gvf_line->{type} = $type;
    $gvf_line->{start} = $start;
    $gvf_line->{end} = $end;
    $gvf_line->{strand} = $strand;

    foreach my $pair (split(';', $attrib)) {
      my ($key, $value) = split('=', $pair);
      if ($key && $value) {
        $gvf_line->{$key} = $value;
      }
    }
    return $gvf_line;
}

sub add_position_and_alleles {
    my ($config, $vcf_line, $gvf_line) = @_;

    my $sa = $config->{slice_adaptor};
    my $ref = $gvf_line->{Reference_seq};
    my $alt = $gvf_line->{Variant_seq};
    my $start = $gvf_line->{start};
    my $end = $gvf_line->{end};
    my $id = $gvf_line->{Dbxref};

    my $seq_region_name = $gvf_line->{seq_id};
    if ($ref eq '.' || $ref eq '~') {
      $ref = _get_seq($config, $seq_region_name, $start, $end);
    }

    my $slice = _get_Slice($config, $seq_region_name);

    $alt =~ s/,/\//g;
    my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
        -start   => $start,
        -end     => $end,
        -slice => $slice,
        -strand  => 1,
        -allele_string => "$ref/$alt",
        -map_weight  => 1,
    );

    my $vcf_record = $vf->to_VCF_record;
    my ($vcf_start, $vcf_ref, $vcf_alt) = ($vcf_record->[1], $vcf_record->[3], $vcf_record->[4]);
    # This can change the start and length of the reference allele
    # We need to correct the ancestral allele accordingly
    my $ancestral_allele = $gvf_line->{ancestral_allele};
    if ($ref ne $vcf_ref || $start != $vcf_start) {
      if ($ancestral_allele && $config->{ancestral_allele_file}) {
        my $ancestral_allele_utils = $config->{ancestral_allele_utils};
        my $ancestral_allele_end = $start + length($vcf_ref) - 1;
        my $vcf_ancestral_allele = $ancestral_allele_utils->assign($seq_region_name, $start, $ancestral_allele_end);
        $gvf_line->{ancestral_allele} = $vcf_ancestral_allele;
      }
    }

    $vcf_line->{'#CHROM'} = $gvf_line->{seq_id};
    $vcf_line->{POS} = $vcf_start;
    $vcf_line->{REF} = $vcf_ref;
    $vcf_line->{ALT} = $vcf_alt;
    $vcf_line->{look_up} = ();
}

sub _get_seq {
  my ($config, $seq_region_name, $start, $end) = @_;
  my $seq;
  my $region = "$seq_region_name:$start-$end";
  if ($config->{fasta_index}) {
    my $fasta_index = $config->{fasta_index};
    $seq = $fasta_index->get_sequence_no_length($region);
  } else {
    my $sa = $config->{slice_adaptor};
    my $slice = $sa->fetch_by_toplevel_location($region);
    $seq = $slice->seq();
  }
  if (!defined $seq) {
    warn "Couldn't get sequence for region $region\n";
    my $repeat = ($end - $start + 1);
    if ($repeat <= 0) {
      warn "Region $region is smaller than or equal to 0\n";
      return '';
    }
    $seq = 'N' x $repeat;
  }
  return $seq;
}

sub _get_Slice {
  my ($config, $seq_region_name) = @_;

  my $slice = $config->{slice}->{$seq_region_name};
  if (!defined $slice) {
    my $sa = $config->{slice_adaptor};
    my $slice = $sa->fetch_by_toplevel_location($seq_region_name);
    $config->{slice}->{$seq_region_name} = $slice;
  }
  return $config->{slice}->{$seq_region_name};
}

sub add_gmaf {
    my ($vcf_line, $gvf_line) = @_;
    # global_minor_allele_frequency=@ 0.435 950
    # global_minor_allele_frequency=0 0.4666 1019
    my $minor_allele = '';
    my ($index, $frequency, $count) = split("\\|", $gvf_line->{global_minor_allele_frequency});
    if ($frequency && $count) {
        if ($index == 0) {
            $minor_allele = $vcf_line->{REF};
        } else {
            my @alleles = split(',', $gvf_line->{Variant_seq});
            $minor_allele = $alleles[$index - 1];
        }
        push @{$vcf_line->{INFO}}, ("MA=$minor_allele", "MAF=$frequency", "MAC=$count");
    }
}

sub add_genotypes {
    my ($vcf_line, $gvf_line) = @_;
    my $genotypes = '';
    my @gvf_alleles = split(',', $gvf_line->{Variant_seq});
    my @gvf_genotype_idxs = split(':', $gvf_line->{Genotype});
    my $look_up = $vcf_line->{look_up};
    my @vcf_alleles = ();
    my $ref = $vcf_line->{REF};
    my $alt = $vcf_line->{ALT};
    push @vcf_alleles, $ref;
    push @vcf_alleles, split(',', $alt);
    my @vcf_gtypes_idxs = ();

    my %look_up_pos = ();
    my $i = 0;
    for my $a (@vcf_alleles) {
        if ($look_up->{$a}) {
            $look_up_pos{$look_up->{$a}} = $i;
        } else {
            $look_up_pos{$a} = $i;
        }
        $i++;
    }
    for my $gt_idx (@gvf_genotype_idxs) {
        my $allele = $gvf_alleles[$gt_idx];
        my $new_idx = $look_up_pos{$allele};
        push @vcf_gtypes_idxs, "$new_idx";
    }
    $genotypes = join('/', @vcf_gtypes_idxs);
    $vcf_line->{FORMAT} = 'GT';
    $vcf_line->{SAMPLE} = $genotypes;
}

sub parse_consequence_info {
  my $gvf_line = shift;
  my $vcf_line = shift;

  my @alleles = split(',', $gvf_line->{Variant_seq});  
  my $id = $vcf_line->{ID};
  my $allele2consequence = {};
  my $is_intergenic = 0;
  my @values = split(',', ($gvf_line->{Variant_effect} || ''));
  #Format=Consequence|Index|Feature_type|Feature_id
  if (@values) {
    foreach my $value (@values) {
      my @split_values = split(' ', $value);
      if (scalar @split_values == 4) {
        my ($consequence, $index, $feature_type, $feature_id) = @split_values;
        my $allele = $alleles[$index];
        $allele2consequence->{$allele}->{$feature_id}->{Consequence} = $consequence;
        $allele2consequence->{$allele}->{$feature_id}->{Feature_type} = $feature_type;
        $allele2consequence->{$allele}->{$feature_id}->{Allele} = $allele;
        $allele2consequence->{$allele}->{$feature_id}->{Feature} = $feature_id;
      } else {
        print STDERR "4 values expected instead got: ", join(', ', @split_values), " for $id\n";
      }
    } 
  } else {
    $is_intergenic = 1;
  } 
   
  @values = split(',', ($gvf_line->{sift_prediction} || ''));
  #Format=Index|Sift_qualitative_prediction|Sift_numerical_value|Feature_id
  #tolerated(0.15)

  if (@values) {
    foreach my $value (@values) {
      my @split_values = split(' ', $value);
      if (scalar @split_values == 4) {
        my ($index, $prediction, $numerical_value, $feature_id) = @split_values;
        my $allele = $alleles[$index];
        $allele2consequence->{$allele}->{$feature_id}->{SIFT} = "$prediction($numerical_value)";
      } else {
        print STDERR "4 values expected instead got: ", join(', ', @split_values), " for $id\n";
      }
    }
  }

  @values = split(',', ($gvf_line->{polyphen_prediction} || ''));
  if (@values) {
    foreach my $value (@values) {
      my @split_values = split(' ', $value);
      if (scalar @split_values == 4) {
        my ($index, $prediction, $numerical_value, $feature_id) = @split_values;
        my $allele = $alleles[$index];
        $allele2consequence->{$allele}->{$feature_id}->{PolyPhen} = "$prediction($numerical_value)";
      } else {
        print STDERR "4 values expected instead got: ", join(', ', @split_values), " for $id\n";
      }
    }
  }
  my $ref_pep = $gvf_line->{reference_peptide};
  @values = split(',', ($gvf_line->{variant_peptide} || ''));
  #Format=Index|Amino_acid|Feature_id
  if (@values) {
    foreach my $value (@values) {
      my @split_values = split(' ', $value);
      if (scalar @split_values == 3) {
        my ($index, $amino_acid, $feature_id) = @split_values;
        my $allele = $alleles[$index];
        $allele2consequence->{$allele}->{$feature_id}->{Amino_acids} = "$ref_pep/$amino_acid";
      } else {
        print STDERR "3 values expected instead got: ", join(', ', @split_values), " for $id\n";
      }
    }
  }
  my @vep_consequence_info = @{$config->{vep}};
  my @consequences = ();
  if ($is_intergenic) {
    foreach my $allele (@alleles) {
      my @fields = ();
      foreach my $info_field (@vep_consequence_info) {
        if ($info_field eq 'Allele') {
          push @fields, $allele;
        } elsif ($info_field eq 'Consequence') {
          push @fields, 'intergenic_variant';
        } else {
          push @fields, '';
        }
      }    
      push @consequences, join('|', @fields);
    }
  } else {
    foreach my $allele (keys %$allele2consequence) {
      foreach my $feature_id (keys %{$allele2consequence->{$allele}}) {
        my @fields = ();
        foreach my $info_field (@vep_consequence_info) {
          my $field_value = $allele2consequence->{$allele}->{$feature_id}->{$info_field} || '';
          push @fields, $field_value;
        }
        push @consequences, join('|', @fields);
      }
    }
  }
  add_info($vcf_line, 'CSQ', join(',', @consequences));
}

sub add_info {
    my ($vcf_line, $key, $value) = @_;
    if ($value) {
      if ($value =~ m/\s/) {
          my $values = join( ',', map { join( '|', split(' ', $_) ) } split(',', $value) );
          push @{$vcf_line->{INFO}}, "$key=$values";
      } else {
          push @{$vcf_line->{INFO}}, "$key=$value";
      }
    }
}

sub add_svs_annotation {
    my ($config, $vcf_line, $gvf_line) = @_;

    my $seq_region_name = $gvf_line->{seq_id};
    my $start = $gvf_line->{start};
    my $end = $gvf_line->{end};
    my $pos = $start - 1 ;
    if ($pos == 0) {
        $pos = 1;
    }
    my $base = _get_seq($config, $seq_region_name, $pos, $pos);
    $vcf_line->{'#CHROM'} = $gvf_line->{seq_id};
    $vcf_line->{REF} = $base;
    $vcf_line->{POS} = $pos;
    my $gvf_svs_type = $gvf_line->{type};
    my $vcf_svs_type = $config->{svs_gvf2vcf}->{$gvf_svs_type};
    $vcf_line->{ALT} = '<' . $vcf_svs_type . '>';
    die "No vcf type for $gvf_svs_type" unless (defined $vcf_svs_type);
    push @{$vcf_line->{INFO}}, ("SVTYPE=$vcf_svs_type", "END=$end");
    if ($gvf_line->{study_accession}) {
        my $study_acc = $gvf_line->{study_accession};
        push @{$vcf_line->{INFO}}, "ST_ACC=$study_acc";
    }
    if ($gvf_line->{Parent}) {
        my $parent = $gvf_line->{Parent};
        push @{$vcf_line->{INFO}}, "Parent=$parent";
    }
    if ($gvf_line->{End_range} || $gvf_line->{Start_range}) {
         push @{$vcf_line->{INFO}}, "IMPRECISE";
    }
}

sub print_vcf_line {
    my ($config, $vcf_line) = @_;
    my $fh = $config->{fh};
    my @oblig_abbrs = ('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER');
    my @add_abbrs = ('FORMAT', 'SAMPLE');

    my @output = ();
    my $all_fields_defined = 1;
    foreach my $abbr (@oblig_abbrs) {
        $all_fields_defined = 0 unless (defined $vcf_line->{$abbr});
        push @output, $vcf_line->{$abbr};
    }
    my $info = join(';', @{$vcf_line->{INFO}}) || '';
    push @output, $info;
    if ($vcf_line->{FORMAT}) {
        foreach my $abbr (@add_abbrs) {
            push @output, $vcf_line->{$abbr};
        }
    }
    if ($all_fields_defined) {
        print $fh join("\t", @output), "\n";
    } else {
        print STDERR join("\t", @output), "\n";
    }
}

sub print_header {
    my $config = shift;

    my $fh = $config->{fh};
    my $sa = $config->{slice_adaptor};
    my $slices = $sa->fetch_all('toplevel', undef, 0, 1);
    my $mca = $slices->[0]->adaptor->db->get_MetaContainerAdaptor;
    my $schema_version = $mca->get_schema_version;
    my $species_name = $mca->get_scientific_name;
    $species_name =~ s/ /_/g;
    $species_name = lc $species_name;
    my ( $sec, $min, $hr, $mday, $mon, $year ) = localtime;
    $year += 1900;    # correct the year
    $mon++;           # correct the month
    my $vcf_file_date = sprintf "%4d%02d%02d", $year, $mon, $mday;
    my $url = '';
    my $ensembl_url = 'https://www.ensembl.org';
    my ($division) = @{$mca->list_value_by_key('species.division')};
    if (!$division) {
      $url = $ensembl_url;
    }
    elsif ($division =~ /^EnsemblVertebrates$/i) {
      $url = 'https://e'.$schema_version.'.ensembl.org/'.$species_name;
    } else {
      $division =~ s/^Ensembl//;
      $division = lc $division;
      $url = 'https://'.$division.'.ensembl.org/'.$species_name;
    }
    my $vcf_source_field = "ensembl;version=$schema_version;url=$url";
    my $reference_info = "https://ftp.ensembl.org/pub/release-$schema_version/fasta/$species_name/dna/";

    # Meta-information
    print $fh join("\n", '##fileformat=VCFv4.1', "##fileDate=$vcf_file_date", "##source=$vcf_source_field", "##reference=$reference_info"), "\n";

    while (my ($source, $desc) = each %{$config->{source_to_desc}}) {
        print $fh "##INFO=<ID=$source,Number=0,Type=Flag,Description=\"$desc\">\n";
    } 

    unless ($config->{structural_variations}) {
        my $abbr = 'TSA';
        my $desc = 'Type of sequence alteration. Child of term sequence_alteration as defined by the sequence ontology project.';
        print $fh "##INFO=<ID=$abbr,Number=1,Type=String,Description=\"$desc\">\n";
    }

    if ($config->{evidence}) {
        my $link = $ensembl_url . '/info/genome/variation/prediction/variant_quality.html#evidence_status';
        while (my($abbr, $evidence_value) = each %{$config->{abbr_to_evidence_value}}) {
            print $fh "##INFO=<ID=$abbr,Number=0,Type=Flag,Description=\"$evidence_value.$link\">\n";
        }
    }

    if ($config->{clinical_significance}) {
        my $link = $ensembl_url . '/info/genome/variation/phenotype/phenotype_annotation.html#clin_significance';
        while (my($abbr, $clin_significance) = each %{$config->{abbr_to_clin_significance}}) {
            print $fh "##INFO=<ID=$abbr,Number=0,Type=Flag,Description=\"$clin_significance.$link\">\n";
        }
    }

    if ($config->{global_maf}) {
        print $fh "##INFO=<ID=MA,Number=1,Type=String,Description=\"Minor Allele\">\n";
        print $fh "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor Allele Frequency\">\n";
        print $fh "##INFO=<ID=MAC,Number=1,Type=Integer,Description=\"Minor Alelele Count\">\n";
    }

    if ($config->{ancestral_allele}) {
        print $fh "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n";
    }

    my $consequence_desc = {};
    if ($config->{incl_consequences}) {
        $consequence_desc->{'VE'} = {
            'Number' => '.',
            'Type' => 'String',
            'Desc' => 'Variant effect of a variant overlapping a sequence feature as computed by the ensembl variant effect pipeline. Format=Consequence|Index|Feature_type|Feature_id. Index indentifies for which variant sequence the effect is described for.',
        };
    }
    if ($config->{protein_coding_details}) {
        $consequence_desc->{'VarPep'} = {
            'Number' => '.',
            'Type' => 'String',
            'Desc' => 'Variant peptide that is translated as a result of a missense variant. Format=Index|Amino_acid|Feature_id. The index identifies the missense variant. The amino acid translated with the missense variant. The feature id for the feature overlapping the variant.',
        };
        $consequence_desc->{'RefPep'} = {
            'Number' => '.',
            'Type' => 'String',
            'Desc' => 'Amino acid translated with reference allele.',
        };
    }
    if ($config->{sift}) {
        $consequence_desc->{'Sift'} = {
            'Number' => '.',
            'Type' => 'String',
            'Desc' => 'Prediction for effect of missense variant on protein function as computed by Sift. Format=Index|Sift_qualitative_prediction|Sift_numerical_value|Feature_id. The index identifies the missense variant. Qualitative prediction is tolerated or deleterious. The numerical value is the normalized probability that the amino acid change is tolerated, so scores nearer 0 are more likely to be deleterious.',
        };
    }
    if ($config->{polyphen}) {
        $consequence_desc->{'Polyphen'} = {
            'Number' => '.',
            'Type' => 'String',
            'Desc' => 'Prediction for effect of missense variant on protein function as computed by Polyphen (human only). Format=Index|Polyphen_qualitative_prediction|Polyphen_numerical_value|Feature_id. The index identifies the missense variant. Qualitative prediction (one of probably damaging, possibly damaging, benign or unknown). Numerical value which is the probability that a substitution is damaging, so values nearer 1 are more confidently predicted to be deleterious.',
        };
    }

    my @abbr_order = ('VE', 'VarPep', 'RefPep', 'Sift', 'Polyphen');
    foreach my $abbr (@abbr_order) {
        if ($consequence_desc->{$abbr}) {
            my $type = $consequence_desc->{$abbr}->{Type};
            my $number = $consequence_desc->{$abbr}->{Number};
            my $desc = $consequence_desc->{$abbr}->{Desc};
            print $fh "##INFO=<ID=$abbr,Number=$number,Type=$type,Description=\"$desc\">\n";
        }
    }

    if ($config->{incl_consequences}) { 
      my @vep_fields = @{$config->{vep}};
      my $format = join('|', @vep_fields);
      print $fh "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl's Variant Effect Pipeline. Format=$format\">\n";
    }

    if ($config->{structural_variations}) {
        print $fh join("\n", (
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
            '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
            '##INFO=<ID=ST_ACC,Number=.,Type=String,Description="Study Accession. Study dbVar ID (estd or nstd)">',
            '##INFO=<ID=Parent,Number=.,Type=String,Description="The structural variant id. It identifies the region of variation.">',
        )), "\n";
        foreach my $name (sort keys %{$config->{header_sv_class}}) {
          my $definition = $config->{header_sv_class}->{$name};
          print $fh "##ALT=<ID=$name,Description=\"$definition\">\n";
        }
    }

    my @header = ('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO');
    print $fh join("\t", @header);
    if ($config->{individual}) {
        print $fh "\tFORMAT\tSAMPLE\n";
    } else {
        print $fh "\n";
    }
}

__END__

=head1 NAME

load_dbsnp.pl

=head1 DESCRIPTION

The script is run as part of the GVF and VCF dumping pipeline. The pipeline generates the list of
required arguments for the script which are required to convert a specific GVF file into a VCF file. 


=head1 SYNOPSIS

gvf2vcf.pl [arguments]

=head1 OPTIONS

=over 4

=item B<--help>

Displays this documentation

=item B<--gvf_file FILE>

GVF file which will be converted into a VCF file

=item B<--vcf_file FILE>

New VCF file

=item B<--species >

Species for which to genereate VCF file

=item B<--registry FILE>

Registry file which provides database connections
to core and variation databases from which the GVF
file was dumped. Database connections are
required for populating meta information in the VCF
file.

=item B<--fasta_file FILE>
Provide fasta file for sequence look ups which
otherwise would be achieved by database queries
which are slower. It is recommended to provide
a fasta file for human file conversions.

=item B<--ancestral_allele|aa>

Parse ancestral allele from GVF file and
store in new VCF file

=item B<--ancestral_allele_file FILE>

Provide the ancestral genome FASTA file for looking up
ancestral alleles. This is required for indels where
the allele string or position of the reference have changed
after formatting alleles to match the VCF format

=item B<--global_maf>

Parse minor allele, minor allele frequency and minor allele count
from GVF file and store in new VCF file

=item B<--evidence>

Parse evidence attributes (Multiple observations, Frequency, Cited,
Phenotype or Disease, 1000 Genomes, ESP, ExAC, gnomAD, TOPMed) from GVF file and store in new VCF file

=item B<--clinical_significance>

Parse clinical significance attributes from GVF file and
sotre in VCF file. Available attributes are described here:
https://www.ensembl.org/info/genome/variation/phenotype/phenotype_annotation.html#clin_significance

=item B<--structural_variations|svs>

Parse structural variants from GVF and store in VCF

=item B<--incl_consequences>

Parse variant consequences from GVF and store in VCF

=item B<--protein_coding_details>

Parse protein consequences from GVF and store in VCF

=item B<--sift>

Parse SIFT annotations from GVF and store in VCF

=item B<--polyphen>

Parse PolyPhen annotations from GVF and store in VCF

=item B<--individual STRING>

Individual name for which genotypes are stored in the GVF file.
The name will be stored in the VCF header line.

=back

=head1 CONTACT
  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.
  Questions may also be sent to the Ensembl help desk at
  <https://www.ensembl.org/Help/Contact>.
=cut

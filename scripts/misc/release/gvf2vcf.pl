#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
  <helpdesk.org>.

=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Slice;
use Bio::EnsEMBL::Utils::Sequence qw(expand reverse_comp);
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

use Getopt::Long;
use FileHandle;
use Compress::Zlib;
use List::Util qw(first);

$| = 1;

my $config = {};
my $gvf_line;
GetOptions(
    $config,
    'gvf_file=s',
    'vcf_file=s',
    'species=s',
    'registry=s',

    'ancestral_allele|aa',
    'global_maf',
    'evidence',
    'clinical_significance',
    'validation_status',

    'structural_variations|svs',

    'incl_consequences',
    'protein_coding_details',
    'sift',
    'polyphen',

    'individual=s',
    'population=s',
    'variation_id',
    'allele_string',
    'set_name',
    'somatic',
    'debug',

) or die "Error: Failed to parse command-line args. Try --help for usage instructions\n";

main($config);

sub main {
    my $config = shift;
    init_db_connections($config);
    init_data($config);
    my $fh_vcf = FileHandle->new('> ' . $config->{vcf_file});
    $config->{fh} = $fh_vcf;
    print_header($config);
    parse_gvf_file($config);
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
    my $table = 'variation';
    if ($config->{structural_variations}) {
        $table = 'structural_variation';
    }
    my $dbh = $vdba->dbc->db_handle;
    my $sth = $dbh->prepare(qq{
        select distinct s.name, s.version, s.description from source s, $table v
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

    $config->{source_to_desc} = $source_to_desc;

    $config->{abbr_to_evidence_value} = {
        'E_Multiple_observations' => 'Multiple_observations',
        'E_Freq' => 'Frequency',
        'E_Hapmap' => 'HapMap',
        'E_1000G' => '1000Genomes',
        'E_Cited' => 'Cited',
        'E_ESP' => 'ESP',
        'E_Phenotype_or_Disease' => 'Phenotype_or_Disease', 
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
            $term->definition =~ m/^"(.*)\."\s\[.*\]$/;
            $definition = $1; 
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
}


sub parse_gvf_file {
    my $config = shift;

    my $gvf_file = $config->{gvf_file};
    my $fh_gvf = gzopen($gvf_file, "rb") or die "Error reading $gvf_file: $gzerrno\n";

    while ($fh_gvf->gzreadline($_) > 0) {
        chomp;
        my $line = $_;
        next if ($line =~ /^##/);
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
            next;
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

        if (0) {
          foreach my $key (qw/Reference_seq Variant_seq variation_id allele_string/) {
            add_info($vcf_line, $key, $gvf_line->{$key});
          }
        }
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

    die "Error reading $gvf_file: $gzerrno\n" if $gzerrno != Z_STREAM_END;
    $fh_gvf->gzclose(); 
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
#    my %attributes = map { split('=', $_) } split(';', $attrib);
#    foreach my $key (keys %attributes) {
#        $gvf_line->{$key} = $attributes{$key};
#    }
    return $gvf_line;
}

sub add_position_and_alleles {
    my ($config, $vcf_line, $gvf_line) = @_;

    my $sa = $config->{slice_adaptor};
    my $ref = $gvf_line->{Reference_seq};
    my $alt = $gvf_line->{Variant_seq};
    my $start = $gvf_line->{start};
    my $end = $gvf_line->{end};

    my $seq_region_name = $gvf_line->{seq_id};
    if ($ref eq '.' || $ref eq '~') {
        my $slice = $sa->fetch_by_toplevel_location("$seq_region_name:$start-$end");
        $ref = $slice->seq();
    }
    my $pos = $start;
    my @alts = split(',', $alt);

    my @alleles = ();
    my %allele_lengths = ();

    push @alleles, $ref;
    # if gvf contains genotypes for heterozygous site variant_seq contains all alleles
    push @alleles, grep {$_ ne $ref} @alts;

    foreach my $allele (@alleles) {
        $allele =~ s/\-//g;
        $allele_lengths{ length($allele) } = 1;
    }
    # in/del/unbalanced
    my $look_up = ();
    if (scalar keys %allele_lengths > 1) {
        # we need ref base before the variation; default to N
        #print STDERR "$ref/$alt $seq_region_name $start $end ", $gvf_line->{Dbxref}, "\n";
        if ($ref ne '-') {
          $pos--;
        }
        my $slice = $sa->fetch_by_toplevel_location("$seq_region_name:$pos-$pos");
        my $prev_base = $slice->seq() || 'N';
        for my $i (0..$#alleles) {
            $alleles[$i] =~ s/\-//g;
            $look_up->{ $prev_base . $alleles[$i] } = $alleles[$i];
            $alleles[$i] = $prev_base . $alleles[$i];
        }
        $ref = shift @alleles;
        #print STDERR "$pos $ref ", join('/', @alleles), "\n";
    } else {
        shift @alleles;
    }
    $alt = join(',', @alleles);

    $vcf_line->{'#CHROM'} = $gvf_line->{seq_id};
    $vcf_line->{POS} = $pos;
    $vcf_line->{REF} = $ref;
    $vcf_line->{ALT} = $alt;
    $vcf_line->{look_up} = $look_up;
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
    my $sa = $config->{slice_adaptor};
    my $slice = $sa->fetch_by_toplevel_location("$seq_region_name:$pos-$pos");
    my $base = $slice->seq() || 'N';
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
    my ( $sec, $min, $hr, $mday, $mon, $year ) = localtime;
    $year += 1900;    # correct the year
    $mon++;           # correct the month
    my $vcf_file_date = sprintf "%4d%02d%02d", $year, $mon, $mday;
    my $vcf_source_field = "ensembl;version=$schema_version;url=http://e$schema_version\.ensembl.org/$species_name";
    my $reference_info = "ftp://ftp.ensembl.org/pub/release-$schema_version/fasta/$species_name/dna/";

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
        my $link = 'http://www.ensembl.org/info/docs/variation/data_description.html#evidence_status';
        while (my($abbr, $evidence_value) = each %{$config->{abbr_to_evidence_value}}) {
            print $fh "##INFO=<ID=$abbr,Number=0,Type=Flag,Description=\"$evidence_value.$link\">\n";
        }
    }

    if ($config->{clinical_significance}) {
        my $link = 'http://www.ensembl.org/info/genome/variation/data_description.html#clin_significance';
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

    if ($config->{population}) {
       print $fh '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">', "\n";
    }

    if ($config->{structural_variations}) {
        print $fh join("\n", (
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description=“Type of structural variant”>',
            '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=“Imprecise structural variation”>',
            '##INFO=<ID=END,Number=1,Type=Integer,Description=“End position of the variant described in this record”>',
            '##INFO=<ID=ST_ACC,Number=.,Type=String,Description=“Study Accession. Study dbVar ID (estd or nstd)”>',
            '##INFO=<ID=Parent,Number=.,Type=String,Description=“The structural variant id. It identifies the region of variation.”>',
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



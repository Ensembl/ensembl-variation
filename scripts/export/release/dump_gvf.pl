#!/usr/bin/env perl
=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut


use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Iterator;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);
use Bio::EnsEMBL::Utils::Sequence qw(expand reverse_comp);
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

use Getopt::Long;
use FileHandle; 
use List::Util qw(first);

$| = 1;

my $config = {};
my $gvf_line;
GetOptions(
$config,
  'help|h',
  'debug',
  'gvf_file=s',
  'species=s',
  'registry=s',

  'seq_region_ids_file=s',

  'seq_region_id=s',
  'seq_region_name=s',
  'slice_piece_name=s',
  'slice_piece_start=s',
  'slice_piece_end=s',
  'is_slice_piece',
  'slice_piece_size=s',

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

  'failed',
  'somatic',

  'print_variation_sets=s', 
  'set_name=s', # phenotype_associated clinically_associated

  'individual=s',

  'population=s',
  'input_gvf=s',
  'frequencies_dir=s', # cache_files_dir?
  'short_name=s',

  'allele_string',
  'variation_id',
) or die "Error: Failed to parse command-line args.\n";

check_arguments($config);
init_db_connections($config);
init_variation_set($config) if (defined $config->{set_name});
init_failed_variations($config) if (defined $config->{failed});
init_sample_data($config) if (defined $config->{individual} || defined $config->{population});
init_consequence_data($config) if ($config->{incl_consequences});
init_slices($config);
$config->{fh} = FileHandle->new($config->{gvf_file}, 'w');
print_header($config);
if ($config->{structural_variations}) {
    dump_svs_data($config);
} else {
    dump_data($config);
}
my $fh = $config->{fh};
$fh->close();

sub check_arguments {
    my $config = shift;
    foreach my $arg (qw/registry species gvf_file/) {
        die "Argument --$arg required, try --help for usage instructions\n" unless $config->{$arg};
    } 
  $config->{slice_piece_size} ||= 1e6;
}

sub init_db_connections {
  my $config = shift;
  my $registry = 'Bio::EnsEMBL::Registry';
  my $registry_file = $config->{registry};
  die "Could not find registry_file $registry_file" unless (-e $registry_file);
  $registry->load_all($registry_file);

  my $vdba = $registry->get_DBAdaptor($config->{species}, 'variation') || usage('Cannot connect to variation db.');
  if ($config->{failed}) {
    $vdba->include_failed_variations(1);
  }
  my $cdba = $registry->get_DBAdaptor($config->{species}, 'core') || usage('Cannot connect to core db.');
  $config->{vdba} = $vdba;
  $config->{cdba} = $cdba;
  $config->{slice_adaptor} = $cdba->get_SliceAdaptor;
  $config->{individual_gt_adaptor} = $vdba->get_IndividualGenotypeAdaptor;
  $config->{individual_adaptor} = $vdba->get_IndividualAdaptor;
  $config->{vf_adaptor}  = $vdba->get_VariationFeatureAdaptor;
  $config->{tv_adaptor}  = $vdba->get_TranscriptVariationAdaptor;
  $config->{variation_set_adaptor} = $vdba->get_VariationSetAdaptor;
  $config->{svf_adaptor} = $vdba->get_StructuralVariationFeatureAdaptor;
  $config->{sv_adaptor} = $vdba->get_StructuralVariationAdaptor;
}

sub init_variation_set {
  my $config = shift;
  my $variation_set_name = $config->{set_name};
  if ($variation_set_name eq 'phenotype_associated') {
    $variation_set_name = 'All phenotype/disease-associated variants';
  }
  if ($variation_set_name eq 'clinically_associated') {
    $variation_set_name = 'Clinically associated variants';
  }
  my $vsa = $config->{variation_set_adaptor};
  my $variation_set = $vsa->fetch_by_name($variation_set_name);
  die "Wrong variation set name: $variation_set_name" unless($variation_set);
  $config->{variation_set} = $variation_set;
}

sub init_failed_variations {
  my $config = shift;
  my ($variation_id, $failed_desc_id, $desc);
  my $failed_variations = {};

  my $vdba = $config->{vdba};
  my $dbh = $vdba->dbc->db_handle;

  my $sth = $dbh->prepare(qq{
    SELECT fv.variation_id, fd.description
    FROM failed_variation fv, failed_description fd
    WHERE fv.failed_description_id = fd.failed_description_id;
  });
  $sth->execute();
  $sth->bind_columns(\($variation_id, $desc));
  while ($sth->fetch()) {
    $failed_variations->{$variation_id}->{$desc} = 1;
  }
  $sth->finish();
  $config->{failed_variations} = $failed_variations;
}

sub init_sample_data {
    my $config = shift;
    if ($config->{individual}) {
        my $name = $config->{individual};
        my $ia = $config->{individual_adaptor};
        my $individuals = $ia->fetch_all_by_name($name);
        die "More than one individual for name $name" if (scalar @$individuals > 1);
        die "No individual for name $name" if (scalar @$individuals == 0);
        $config->{_individual} = $individuals->[0];
    }
}

sub init_consequence_data {
    my $config = shift;
    my @protein_func_pred_tools = ();
    if ($config->{sift} || $config->{polyphen}) {
        $config->{protein_coding_details} = 1;
        push @protein_func_pred_tools, 'sift' if ($config->{sift});
        push @protein_func_pred_tools, 'polyphen' if ($config->{polyphen});
    }
    $config->{protein_func_pred_tools} = \@protein_func_pred_tools;
}

sub init_slices {
  my $config = shift;
  my $sa = $config->{slice_adaptor};
  my $slices = [];
  if ($config->{seq_region_ids_file}) {
    my $seq_region_file = $config->{seq_region_ids_file};
    die "Could not find seq_region_file $seq_region_file" unless (-e $seq_region_file);
    my $fh = FileHandle->new($seq_region_file, 'r');
    while (<$fh>) {
      chomp;
      my $slice = $sa->fetch_by_seq_region_id($_);
      push @$slices, $slice;
    } 
    $fh->close();
  } elsif ($config->{is_slice_piece}) {
    my $seq_name = $config->{seq_region_name};
    my $start = $config->{slice_piece_start};
    my $end = $config->{slice_piece_end};
    my $slice = $sa->fetch_by_region('toplevel', $seq_name, $start, $end, 1);
    push @$slices, $slice;
  } else {
    $slices = $sa->fetch_all('toplevel', undef, 0, 1);
  }
  $config->{slices} = $slices;
}

sub dump_svs_data {
    my $config = shift;

    my $slices = $config->{slices};
    my $svfa = $config->{svf_adaptor};
    my $sva = $config->{sv_adaptor};

    my $prev_svs;
    my $id_count = 0;
    
    foreach my $slice (@$slices) {
        print STDERR join(' ', $slice->seq_region_name, $slice->start, $slice->end), "\n";
        my $svfs = $svfa->fetch_all_by_Slice($slice, 1);
        for my $svf (@$svfs) {
            next if ($svf->seq_region_start <= $slice->start); # avoid duplicated lines caused by vf overlapping two slice pieces

            my $gvf_line = {};
            # ignore CNV probes
            next if $svf->var_class eq 'CNV_PROBE';

            $gvf_line->{seqid} = $svf->slice->seq_region_name;
            my $start = $svf->seq_region_start;
            my $end   = $svf->seq_region_end;
            if ($start > $end) {
                $start = $end;
            }
            $gvf_line->{start}  = $start;
            $gvf_line->{end}    = $end;
            $gvf_line->{strand} = $svf->strand == 1 ? '+' : ($svf->strand == -1 ? '-' : '.');
            $gvf_line->{type}   = $svf->class_SO_term;
            my $source          = $svf->source_name;
            $gvf_line->{source} = $source;
            $source .= '_' . $svf->structural_variation->source_version if defined $svf->structural_variation->source_version;
            $gvf_line->{attributes}->{Dbxref} = "$source:" . $svf->variation_name;
            $gvf_line->{score} = '.';
            $gvf_line->{phase} = '.';

            my $coords = join('-', $svf->seq_region_name, $svf->seq_region_start, $svf->seq_region_end);
            if (my $prev_coords = $prev_svs->{$svf->variation_name}) {
                warn "repeated SV: " . $svf->variation_name . " coords 1: $prev_coords, coords 2: $coords\n" if $prev_coords ne $coords;
                next;
            }
            my $sv = $svf->structural_variation;
            $gvf_line->{attributes}->{study_accession} = $svf->study->name if $svf->study;

            if ((defined $svf->inner_start) && (defined $svf->outer_start) && ($svf->inner_start != $svf->outer_start)) {
                $gvf_line->{attributes}->{Start_range} = join(',', $svf->outer_start, $svf->inner_start);
            }
            if ((defined $svf->inner_end) && (defined $svf->outer_end) && ($svf->inner_end != $svf->outer_end)) {
                $gvf_line->{attributes}->{End_range} = join(',', $svf->inner_end, $svf->outer_end);
            }
            if ($sv) {
                if (ref $sv eq 'Bio::EnsEMBL::Variation::SupportingStructuralVariation') {
                    if (my $parents = $sva->fetch_all_by_supporting_evidence($sv)) {
                        next unless(scalar @$parents);
                        $gvf_line->{attributes}->{Parent} = join(',', map { $_->variation_name } @$parents);
                    }
                }
            }

            # we can now have SVs that map to multiple locations so we can't use the
            # feature's own identifier and we have to use a file-wide count as for
            # normal variations
            $gvf_line->{attributes}->{ID} = ++$id_count;
            print_gvf_line($config, $gvf_line);
            $prev_svs->{$svf->variation_name} = $coords if $gvf_line;
        }
    }
}

sub dump_data {
    my $config = shift;
    my $vfa = $config->{vf_adaptor};
    my $slice_adaptor = $config->{slice_adaptor};
    my $slices = $config->{slices};
    my $vf_it;
    my $max_length = $config->{slice_piece_size};
    my $overlap = 1;

    foreach my $slice (@$slices) {
      my $full_slice = $slice_adaptor->fetch_by_seq_region_id($slice->get_seq_region_id);
      my $slice_end = $full_slice->seq_region_end;
        my $slice_pieces = split_Slices([$slice], $max_length, $overlap);
        foreach my $slice_piece (@$slice_pieces) {
          my $slice_piece_start = $slice_piece->seq_region_start;
          my $slice_piece_end = $slice_piece->seq_region_end;
            if ($config->{individual}) { 
                update_gts($config, $slice_piece);
            }
            if ($config->{somatic}) {
                my $vfs = $vfa->fetch_all_somatic_by_Slice($slice_piece);
                $vf_it = Bio::EnsEMBL::Utils::Iterator->new($vfs);
            } elsif ($config->{set_name}) {
                my $variation_set = $config->{variation_set};
                my $vfs = $variation_set->get_all_VariationFeatures_by_Slice($slice_piece);
                $vf_it = Bio::EnsEMBL::Utils::Iterator->new($vfs);
            } else {
                $vf_it = $vfa->fetch_Iterator_by_Slice($slice_piece);
            }

            if ($config->{incl_consequences} || $config->{protein_coding_details}) {
                my @vfs = ();
                my $gvf_lines = {};
                my $count = 0;
                while (my $vf = $vf_it->next) {
                  my $vf_start = $vf->seq_region_start;
                  my $vf_end = $vf->seq_region_end;
                  if ($vf_end < $vf_start) {
                    ($vf_start, $vf_end) = ($vf_end, $vf_start);
                  }
                    next if ($vf_start == $slice_piece_end && $vf_end >= $slice_piece_end && $slice_piece_end != $slice_end);
                    next if ($vf_end >= $slice_piece_end);
                    push @vfs, $vf;
                    $count++;
                    if ((!$vf_it->peek) || $count == 1000) {
                        if ($config->{incl_consequences}) {
                            $gvf_lines = add_variant_effect($config, \@vfs);
                        } elsif ($config->{population}) {
                            $gvf_lines = add_frequencies($config, \@vfs);
                        }
                        annotate_gvf_lines($config, $gvf_lines);
                        @vfs = ();
                        $count = 0;
                    }
                }
            } else {
                while (my $vf = $vf_it->next) {
                  my $vf_start = $vf->seq_region_start;
                  my $vf_end = $vf->seq_region_end;
#                   slice_piece_start = 1 slice_piece_end = 5 vf is insertion vf_start = 5 vf_end = 4 -> not included in next slice_piece slice_piece_start = 5 slice_piece_end = 9  
#                   slice_piece_start = 1 slice_piece_end = 5 vf is deletion vf_start = 4 vf_end = 5 -> included in next slice_piece slice_piece_start = 5 slice_piece_end = 9 -> need to be filtered out
                    next if ($vf_start == $slice_piece_end && $vf_end >= $slice_piece_end && $slice_piece_end != $slice_end);
                    next if ($vf_end >= $slice_piece_end);
                    my $gvf_line = {};
                    annotate_vf($config, $gvf_line, $vf);
                    print_gvf_line($config, $gvf_line) if ((scalar keys %$gvf_line) > 1);
                }
            }
        } 
    }
}

sub annotate_gvf_lines {
    my ($config, $gvf_lines) = @_;
   
    foreach my $vf_id (keys %$gvf_lines) {
        my $gvf_line = $gvf_lines->{$vf_id};
        my $vf = $gvf_lines->{$vf_id}->{vf};
        annotate_vf($config, $gvf_line, $vf);
        print_gvf_line($config, $gvf_line);
    }
}

sub annotate_vf {
    my ($config, $gvf_line, $vf) = @_;
    if ($config->{failed}) {
        my $failed_variations = $config->{failed_variations};
        return unless (exists $failed_variations->{$vf->{_variation_id}});
        my $failed_descs = join(',', keys %{ $failed_variations->{ $vf->{_variation_id} } });
        $gvf_line->{attributes}->{ensembl_failure_reason} = $failed_descs;
    }

    $gvf_line->{type} = $vf->class_SO_term;
    $gvf_line->{score} = '.';
    $gvf_line->{phase} = '.';

    add_source($gvf_line, $vf);
    add_coords($gvf_line, $vf);
    add_alleles($gvf_line, $vf);

    if ($config->{individual}) {
        return unless (add_genotype($config, $gvf_line, $vf));
    }

    if ($config->{ancestral_allele}) {
        add_ancestral_allele($gvf_line, $vf);
    }
    if ($config->{evidence}) {
        add_evidence($gvf_line, $vf);
    }
    if ($config->{clinical_significance}) {
        add_clinical_significance($gvf_line, $vf);
    }
    if ($config->{global_maf}) {
        add_global_maf($gvf_line, $vf);
    }
    if ($config->{variation_id}) {
        my $variation = $vf->variation;
        $gvf_line->{attributes}->{variation_id} = $variation->dbID;
    } 

    if (defined $gvf_line->{attributes}->{invalid_variant_strings}) {
        $gvf_line = {};
    } else {
        $gvf_line->{attributes}->{ID} = ++$config->{id_count};
    }

}

sub add_source {
    my ($gvf_line, $vf) = @_;
    my $source          = $vf->source->name;
    $gvf_line->{source} = $source;
    $source .= '_' . $vf->source_version if defined $vf->source_version;
    $gvf_line->{attributes}->{Dbxref} = "$source:" . $vf->variation_name;
}

sub add_coords {
    my ($gvf_line, $vf) = @_;
    $gvf_line->{seqid} = $vf->slice->seq_region_name;
    my $start = $vf->seq_region_start;
    my $end = $vf->seq_region_end;
    if ($start > $end) {
        $start = $end;
    }
    $gvf_line->{start} = $start;
    $gvf_line->{end} = $end;
    $gvf_line->{strand} = $vf->strand == 1 ? '+' : ($vf->strand == -1 ? '-' : '.');
}

sub add_alleles {
    my ($gvf_line, $vf) = @_;
    my @alleles = split /\//, $vf->allele_string;
    map {expand(\$_)} @alleles;
    $gvf_line->{allele_string} = join('/', @alleles);
    if ($config->{allele_string}) { 
        $gvf_line->{attributes}->{allele_string} = join(',', @alleles);
    }
    my $ref_seq = shift @alleles unless @alleles == 1; # shift off the reference allele

    if ($vf->allele_string eq 'HGMD_MUTATION') {
        $gvf_line->{attributes}->{'comment'} = 'HGMD_MUTATION';
        @alleles = ();
        if ($vf->var_class eq 'deletion') {
            $ref_seq = '.'; # get slice?
            push @alleles, '-';
        } elsif ($vf->var_class eq 'insertion') {
            $ref_seq = '-';
            push @alleles, '.';
        } else {
            $ref_seq = '.'; # get slice?
            push @alleles, '.';
        }
    }
    # transform to valid variant and ref sequences

    my $Variant_seq_regex = qr/^([A-DGHKMNR-WY]+| # Any valid IUPAC Nucleotide
                                        ~\d*| # A ~ optionally followed by an integer
                                    [.\-!^\*] # Any [.-!^*]
                             )$/ix;  # from gvf_validator

    my $invalid_variant_seq = 0;
    my @invalid_variant_strings = ();
    for my $i (0..$#alleles) {
        unless ($alleles[$i] =~ $Variant_seq_regex) {
            push @invalid_variant_strings, $alleles[$i];
            $alleles[$i] = '.';
            $invalid_variant_seq = 1;
        }
    }
    if ($invalid_variant_seq) {
        $gvf_line->{attributes}->{invalid_variant_strings} = join(',', @invalid_variant_strings);
    }
    $gvf_line->{attributes}->{Variant_seq} = join ',', @alleles;

    my $Reference_seq_regex = qr/^([A-DGHKMNR-WY]+|~\d*|[.\-])$/i; # from gvf_validator
    my $sa = $config->{slice_adaptor};
    $ref_seq ||= '.';
    unless ($ref_seq =~ $Reference_seq_regex) {
        $gvf_line->{attributes}->{invalid_ref_seq} = $ref_seq;
        my $location = $gvf_line->{seqid} . ':' . $gvf_line->{start} . '-' . $gvf_line->{end};
        my $ref_seq_slice = $sa->fetch_by_toplevel_location($location);
        $ref_seq = $ref_seq_slice->seq;
    }
    $ref_seq = '~' if (CORE::length($ref_seq) > 50);
    $gvf_line->{attributes}->{Reference_seq} = $ref_seq;
}

sub add_ancestral_allele {
    my ($gvf_line, $vf) = @_;
    my $variation = $vf->variation;
    if (defined($variation->ancestral_allele)) {
        $gvf_line->{attributes}->{ancestral_allele} = $variation->ancestral_allele;
    }
}

sub add_evidence {
    my ($gvf_line, $vf) = @_;
    my $variation = $vf->variation;
    my $values = $variation->get_all_evidence_values();
    if (scalar @$values) {
        $gvf_line->{attributes}->{evidence_values} = join(',', @$values);
    }
}

sub add_clinical_significance {
    my ($gvf_line, $vf) = @_;
    my $variation = $vf->variation;
    my @states = @{$variation->get_all_clinical_significance_states};
    if (scalar @states) {
        $gvf_line->{attributes}->{clinical_significance} = join(',', @states);
    }
}

sub add_global_maf {
    my ($gvf_line, $vf) = @_;
    my $variation = $vf->variation;
    if ( defined( $variation->minor_allele_frequency)) {
        my @alleles = split /\//, $gvf_line->{allele_string};
        my $allele_idx = first { $alleles[$_] eq $variation->minor_allele } 0..$#alleles;
        if (defined($allele_idx)) {
            $gvf_line->{attributes}->{global_minor_allele_frequency} =
                join('|',
                    $allele_idx,
                    $variation->minor_allele_frequency,
                    $variation->minor_allele_count,);
        }
    }
}

sub add_variant_effect {
    my ($config, $vfs) = @_;

    my $tva = $config->{tv_adaptor};
    my $gvf_lines = {};
    foreach my $vf (@$vfs) {
        my $vf_id = $vf->dbID;
        $gvf_lines->{$vf_id}->{vf} = $vf;
    }

    foreach my $tv (@{$tva->fetch_all_by_VariationFeatures($vfs)}) {
        my $tv_stable_id = $tv->transcript_stable_id;
        my $vf = $tv->variation_feature;
        my $vf_id = $vf->dbID;
        $gvf_lines->{$vf_id}->{vf} = $vf;
        my @alleles = split /\//, $vf->allele_string;
        my $ref_seq = shift @alleles unless @alleles == 1; # shift off the reference allele
        
        if ($config->{protein_coding_details}) {
            my $ref_tva = $tv->get_reference_TranscriptVariationAllele;
            if (my $pep = $ref_tva->peptide) {
                $gvf_lines->{$vf_id}->{attributes}->{reference_peptide} = $pep;
            }
        }

        foreach my $tv_allele (@{$tv->get_all_alternate_TranscriptVariationAlleles}) {
            my $allele_idx = first { $alleles[$_] eq $tv_allele->variation_feature_seq } 0..$#alleles;
            if ($config->{incl_consequences}) {
                for my $oc (@{$tv_allele->get_all_OverlapConsequences}) {
                    next unless (defined $allele_idx);
                    push @{ $gvf_lines->{$vf_id}->{attributes}->{Variant_effect} ||= [] },
                        join(' ', $oc->SO_term,
                                  $allele_idx,
                                  $oc->feature_SO_term,
                                  $tv_stable_id,);
                }
            }
            if ($config->{protein_coding_details}) {
                if ($tv_allele->pep_allele_string) {
                    push @{ $gvf_lines->{$vf_id}->{attributes}->{variant_peptide} ||= [] },
                        join(' ',
                            $allele_idx,
                            $tv_allele->peptide,
                            $tv_stable_id,);

                    for my $tool (@{$config->{protein_func_pred_tools}}) {
                        my $pred_meth = "$tool\_prediction";
                        my $score_meth = "$tool\_score";
                        if (my $pred = $tv_allele->$pred_meth) {
                            $pred =~ s/\s/_/g;
                            push @{ $gvf_lines->{$vf_id}->{attributes}->{"$tool\_prediction"} ||= [] },
                                join(' ',
                                    $allele_idx,
                                    $pred,
                                    $tv_allele->$score_meth,
                                    $tv_stable_id,);
                        }
                    }
                }
            }
        } # end foreach transcript_variation_allele
    } # end foreach transcript_variation

    return $gvf_lines;
}

sub update_gts {
    my ($config, $slice) = @_;
    my $igta = $config->{individual_gt_adaptor};
    my $individual = $config->{_individual};
    my $igts = $igta->fetch_all_by_Slice($slice, $individual);
    my $gt_hash = {};
    for my $igt (@$igts) {
        my $key = join('_',
            $igt->slice->seq_region_name,
            $igt->seq_region_start,
            $igt->seq_region_end,
        );
        my $gts = $gt_hash->{$key} ||= [];
        my $seen = 0;
        if (@$gts) {
            for my $gt (@$gts) {
                if (($gt->allele1 eq $igt->allele1 && $gt->allele2 eq $igt->allele2) ||
                    ($gt->allele1 eq $igt->allele2 && $gt->allele2 eq $igt->allele1)) {
                    $seen++;
                    last;
                }
            }
        }
        push @$gts, $igt unless $seen;
    }
    $config->{gts} = $gt_hash;
}

sub add_genotype {
    my ($config, $gvf_line, $vf) = @_;

    my $gt_hash = $config->{gts};
    my $key = join('_', $vf->slice->seq_region_name, $vf->seq_region_start, $vf->seq_region_end,);
    my @gts = @{ $gt_hash->{$key} || [] };
    GT : for my $gt (@gts) {
        my $gt1 = $gt->allele1;
        my $gt2 = $gt->allele2;

        unless ($vf->strand == $gt->strand) {
            reverse_comp(\$gt1);
            reverse_comp(\$gt2);
        }
        my @alleles = split /\//, $vf->allele_string;
        my $ref_allele = shift @alleles;
        my $hom = $gt1 eq $gt2;
        my $gt1_is_ref = $gt1 eq $ref_allele;
        my $gt2_is_ref = $gt2 eq $ref_allele;

        my $gt1_in_alt;
        my $gt2_in_alt;

        if ($hom && $gt1_is_ref) {
            # homozygous for the reference, don't include
            # next VF;
            return 0;
        }
        my $variant_seq;

        for my $allele (@alleles) {
            $gt1_in_alt = 1 if $allele eq $gt1;
            $gt2_in_alt = 1 if $allele eq $gt2;
            last if ($gt1_in_alt && $gt2_in_alt);
        }

        if ($hom && $gt1_in_alt) {
            # homozygous for alt
            $variant_seq = "$gt1";
        } elsif ($gt1_is_ref && $gt2_in_alt) {
            # het for ref and alt
            $variant_seq = "$gt1,$gt2";
        } elsif ($gt2_is_ref && $gt1_in_alt) {
            # het for ref and alt
            $variant_seq = "$gt2,$gt1";
        } elsif ($gt1_in_alt && $gt2_in_alt) {
            # het for 2 alts
            $variant_seq = "$gt1,$gt2";
        } else {
            # the genotype sequence isn't in the allele string
            # first check if the classes match
            my $gt_incl_indel = ($gt1 eq '-' || $gt eq '-');
            my $vf_is_indel = ($vf->var_class eq 'deletion' || $vf->var_class eq 'insertion');
            next GT unless $gt_incl_indel == $vf_is_indel;
            my $vf_is_sub = $vf->var_class eq 'substitution';
            my $gt_is_sub = (length $gt1 > 1) || (length $gt2 > 1);
            next GT unless $vf_is_sub == $gt_is_sub;

            # sometimes there are multiple genotypes at the same location,
            # but this only seems to happen when the associated variation
            # differs

            next GT if ($gt->variation->name ne $vf->variation_name);

            # otherwise sometimes the variations match, but the alleles
            # don't match the genotype e.g. rs55844409 has allele string
            # T/C but Watson has a G/G genotype. We suspect these are dbSNP
            # bugs, but we still include them

            $variant_seq = $hom ? $gt1 : "$gt1,$gt2";
        }

        $gvf_line->{attributes}->{Variant_seq} = $variant_seq;
        $gvf_line->{attributes}->{Zygosity} = $hom ? 'homozygous' : 'heterozygous';
        my @variant_seqs = split ',', $variant_seq;
        my $index = 0;
        my %allele_index = map {$_ => $index++} @variant_seqs;
        my $genotype = $allele_index{$gt1} . ":" . $allele_index{$gt2};
        $gvf_line->{attributes}->{Genotype} = $genotype;
    }

    unless ($gvf_line->{attributes}->{Zygosity}) {
        # this vf isn't in this individual, so don't output it
        # next VF;
        $gvf_line = {};
        return 0;
    }
    return 1;
}

sub print_header {
    my $config = shift;
    
    my $gff_version = '##gff-version 3';
    my $gvf_version = '##gvf-version 1.07';

    # build up a date string in the format specified by the GFF spec
    my ( $sec, $min, $hr, $mday, $mon, $year ) = localtime;
    $year += 1900;    # correct the year
    $mon++;           # correct the month
    my $date = sprintf "%4d-%02d-%02d", $year, $mon, $mday;

    my $slices = $config->{slices};
    my $assembly = $slices->[0]->coord_system->version;
    my $mca = $slices->[0]->adaptor->db->get_MetaContainerAdaptor;
    my $tax_id = $mca->get_taxonomy_id;

    my $gff_header = "##file-date $date\n"
    . "##genome-build ensembl $assembly\n"
    . "##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$tax_id";

    my $schema_version = $mca->get_schema_version;
    my $species_name = $mca->get_scientific_name;
    $species_name =~ s/ /_/g;
    my $url = 'http://e'.$schema_version.'.ensembl.org/'.$species_name;
    my $gvf_header = "##feature-ontology http://song.cvs.sourceforge.net/viewvc/song/ontology/so.obo?revision=1.283\n"
    . "##data-source Source=ensembl;version=$schema_version;url=$url\n"
    . "##file-version $schema_version\n";

    if ($config->{individual}) {
        $gvf_header .= "##individual-id " . $config->{individual};
    }
    
    if ($config->{population}) {
        $gvf_header .= "##population " . $config->{population};
    }

    my $fh = $config->{fh};
    print $fh join("\n", $gff_version, $gvf_version, $gff_header, $gvf_header);
    for my $slice (@$slices) {
        print $fh join(' ', '##sequence-region', $slice->seq_region_name, $slice->start, $slice->end), "\n";
    }

}

sub print_gvf_line {
    my ($congig, $gvf_line) = @_;
    my $fh = $config->{fh};  
     
    my $line = join("\t", $gvf_line->{seqid}, $gvf_line->{source}, $gvf_line->{type}, $gvf_line->{start},
            $gvf_line->{end}, $gvf_line->{score}, $gvf_line->{strand}, $gvf_line->{phase});

    my @attributes = ();
    if ($gvf_line->{attributes}) {
        for my $key (keys %{$gvf_line->{attributes}}) {
            my $val = $gvf_line->{attributes}->{$key};
            if (ref $val eq 'ARRAY') {
                push @attributes, $key . '=' . join(',', @$val);
            } else {
                push @attributes, $key . '=' . $val;
            }
        }
        $line .= "\t" . join(';', @attributes);
    }
    print $fh $line, "\n";
}


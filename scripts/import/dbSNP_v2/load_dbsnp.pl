#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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

use strict;
use warnings;

use Data::Dumper;

use Getopt::Long;
use Pod::Usage qw(pod2usage);

use Bio::EnsEMBL::Variation::Utils::QCUtils qw(count_rows get_evidence_attribs check_illegal_characters check_for_ambiguous_alleles remove_ambiguous_alleles);
use Bio::EnsEMBL::Registry;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::Variation::Utils::AncestralAllelesUtils;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(SO_variation_class);
use Bio::EnsEMBL::Variation::Utils::Constants qw(:SO_class_terms);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use File::Basename;
use POSIX qw(strftime);
use JSON;
use FileHandle;

my $config = configure();

print "species = $config->{'species'}\n";
print "registry_file = $config->{'registry'}\n";

my $debug = $config->{'debug'};
my $lines = 0;

# TODO these need to be command line options
my $add_map_weight = 1;
my $max_lines = 5000000;

init_reports($config);

init_db_connections($config);
my $dbh_var = $config->{'dbh_var'};
my $dbh_core= $config->{'dbh_core'};

get_db_info($dbh_var);
get_db_info($dbh_core);

$config->{'source_id'} = get_sources($dbh_var);

# For each load get the fasta index and pass this around
my $fai_index;
if ($config->{'ref_check'}) {
  print "ref checking done\n";
  $fai_index = Bio::DB::HTS::Faidx->new($config->{'fasta_file'});
  print "fai_index done\n";
  die("Unable to get FASTA index") if (!$fai_index);
} else {
  print "ref checking not done\n";
}

my $ancestral_fai_index;
my $ancestral_alleles_utils;
if ($config->{'assign_ancestral_allele'}) {
  print "ancestral allele assignment done\n";
  $ancestral_fai_index = Bio::DB::HTS::Faidx->new($config->{'ancestral_fasta_file'});
  $ancestral_alleles_utils = Bio::EnsEMBL::Variation::Utils::AncestralAllelesUtils->new(-fasta_db => $ancestral_fai_index);
  print "ancestral_fai_index done\n";
  die("Unable to get ancestral FASTA index") if (!$ancestral_fai_index);
} else {
  print "ancestral allele assignment not done\n";
}

my $seq_regions_names = get_seq_region_names($dbh_var);
my $nc_regions = get_nc_regions($dbh_var);
my $nw_regions = get_nw_regions($dbh_var);
my $nt_regions = get_nt_regions($dbh_var);
my $ref_regions = get_ref_regions($dbh_core) if ($add_map_weight);
my $lu_info = get_lu_info($dbh_var);

# Parsing the data file
my ($num_lines) = parse_dbSNP_file($config);

# Close open filehandles used for tracking the
# update of any minor allele changes that are only done for
# GRCh38 imports.
if (($config->{'assembly'} eq 'GRCh38') && ($config->{'add_maf'})) {
  for my $ma_type ('update', 'log') {
    my $fh = $config->{join('-', 'ma', $ma_type, 'fh')};
    $fh->close();
  }
}

report_summary($config, $num_lines);

# Set up the reporting files
sub init_reports {
  my ($config) = @_;

  # Set up an error file
  my $rpt_dir = $config->{'rpt_dir'};
  my @suffixes = ('.json', '.json.gz', '.json.bz2');
  my ($base_filename, $path, $suffix) = fileparse($config->{'input_file'}, @suffixes);
  my $error_filename = join('.', $base_filename, 'errors');
  my $rpt_filename = join('.', $base_filename, 'summary');

  print "error_file = ($error_filename)\n";
  print "rpt_filename = ($rpt_filename)\n";
  print "rpt_dir = ($rpt_dir)\n";

  $config->{'error_file'} = join('/', $rpt_dir, $error_filename);
  $config->{'rpt_file'} = join('/', $rpt_dir, $rpt_filename);
  $config->{'data_file_short'} = join('', $base_filename, $suffix);

  print "rpt_file ($config->{'rpt_file'})\n";

  if (-e $config->{'rpt_file'}) {
    die("$config->{'rpt_file'} already exists. Please rename or delete\n");
  }

  print "error_file ($config->{'error_file'})\n";
  if (-e $config->{'error_file'}) {
    die("$config->{'error_file'} already exists. Please rename or delete\n");
  }

  # For GRCh37 the 1000Genomes minor allele is on the forward strand
  # Some regions between GRCh37 and GRCh38 have been reverse complimented
  # causing miss-matches with the minor allele.
  # Variants with different strand mappings between GRCh37 and GRCh38
  # are identified and the minor allele is reverse complimented.
  # Files for updates and logs for minor allele mismatches.
  if (($config->{'assembly'} eq 'GRCh38') && ($config->{'add_maf'})) {
    print "Flip 1000Genomes minor alleles for strand differences between assembly $config->{'assembly'} and GRCh37\n";
    for my $ma_type ('update', 'log') {
      my $filename = join('-', $base_filename, 'ma', $ma_type) . '.txt';
      if (-e "$rpt_dir/$filename") {
        die("$rpt_dir/$filename already exists. Please rename or delete\n");
      }
      my $fh = FileHandle->new("$rpt_dir/$filename", 'w');
      if (! $fh) {
        die("Unable to open $rpt_dir/$filename");
      } else {
        print "minor allele $ma_type file = $rpt_dir/$filename\n";
        $config->{join('-', 'ma', $ma_type, 'fh')} = $fh;
      }
    }
  } else {
    print "No flip of 1000Genomes minor alleles for assembly $config->{'assembly'}\n";
  }
}

# Sets the variation database connections
# Assumes registry file exists (checked in configure)
sub init_db_connections {
  my ($config) = @_;

  my $reg = 'Bio::EnsEMBL::Registry';
  $reg->no_version_check(1); 
  $reg->load_all($config->{'registry'});

  my $dbVar = $reg->get_DBAdaptor($config->{'species'}, 'variation') || die "Error getting db adaptor variation";
  $config->{'dbh_var'}  = $dbVar->dbc->db_handle;

  $config->{'attribute_adaptor'} = $reg->get_adaptor($config->{'species'}, 'variation', 'attribute') || die "Error getting attribute adaptor";
  
  my $dbCore = $reg->get_DBAdaptor($config->{'species'}, 'core') || die "Error getting db adaptor core";
  $dbCore->dbc->disconnect_if_idle;
  $config->{'core_adaptor'} = $dbCore;
  $config->{'dbh_core'} = $dbCore->dbc->db_handle;
}

sub parse_dbSNP_file {
  my ($config) = @_;
 
  my $inputfile = $config->{'input_file'};
 
  my $refSNP_data;

  debug("\n>>>> Parsing input file <<<<") if ($debug);
  
  my $in;

  if (! -e $inputfile) {
      die("$inputfile does not exist");
  }

  # Note zcat not dying if file does not exist
  if ($inputfile =~ /gz$/) {
    open($in, "zcat $inputfile 2>&1 |") or die("Could not open $inputfile for reading");
  } elsif ($inputfile =~/bz2$/) {
    open($in, "bzcat $inputfile 2>&1 |") or die("Could not open $inputfile for reading");
  } else {
    open($in, '<', $inputfile) or die("Unable to open $inputfile: $!");
  }

  my $batch_id = get_batch_id($config->{'dbh_var'}, $config->{'data_file_short'});
  die("no batch_id for $config->{'data_file_short'}") if (! $batch_id);
 
  $config->{'batch_id'} = $batch_id;
 
  while (my $json_string = <$in>) {

    if ($lines >= $max_lines) {
      print "Max lines ($max_lines) reached so stopping\n";
      last;
    }
    $lines++;

    my $rs_obj = JSON->new->decode($json_string) or throw("ERROR: Failed to read file $inputfile");
    my $refsnp_data = parse_refsnp($config, $rs_obj);
    next if (!$refsnp_data);
    print Dumper($refsnp_data) if (!$config->{'db_load'});

    qc_refsnp($config, $lu_info, $refsnp_data, $fai_index);
    print Dumper($refsnp_data) if (!$config->{'db_load'});
    
    set_SO_variation_class($config, $refsnp_data); 
    import_refsnp($config->{'dbh_var'},$refsnp_data, $config) if ($config->{'db_load'});
  }
  print "Number of lines processed ($lines)\n";
  close($in);
  return $lines;
}

# Parses the JSON structure to obtain the data to stored in the database
sub parse_refsnp {
  my ($config, $rs_json) = @_;
  my $data;
  
  $data->{'refsnp_id'} = $rs_json->{'refsnp_id'};

  # Get the variant_type
  if (exists $rs_json->{'primary_snapshot_data'}) {
    $data->{'variant_type'} = $rs_json->{'primary_snapshot_data'}->{'variant_type'};
  } else {
    $data->{'variant_type'} = undef;
  }   

  # Variation_data
  $data->{'v'}->{'name'} = 'rs' . $rs_json->{'refsnp_id'};
  $data->{'v'}->{'source_id'} = $config->{'source_id'}->{'dbSNP'};
  
  # Variation feature data - location
  $data->{'vfs'} = get_variant_features($rs_json);

  # Merges
  $data->{'merges'} = get_merges($rs_json);
  $data->{'dbsnp2_merges'} = get_dbsnp2_merges($rs_json);

  # HGVS
  $data->{'hgvs'} = get_hgvs($config, $rs_json);

  # Citations
  $data->{'citations'} = get_citations($rs_json);

  # Data
  ($data->{'snp_support'}, $data->{'freq_support'})  = get_support($rs_json);

  # 1000Genomes data
  # Do not add 1000Genomes data
  # $data->{'1000Genomes'} = get_study_frequency($rs_json, '1000Genomes');
  if ($config->{'add_maf'}) {
    $data->{'1000Genomes'} = get_study_frequency($rs_json, '1000Genomes');
  }

  # If: 
  # - the assembly is GRCh38
  # - the 1000Genomes data has a minor_allele
  # check if a flip is needed. 
  # TODO - add a flag if flipping should be done if for GRCh38
  #        no flip is needed
  if ( ($config->{'add_maf'}) &&
       ($config->{'assembly'} eq 'GRCh38') &&
       defined $data->{'1000Genomes'} &&
       $data->{'1000Genomes'}->{'minor_allele'}) {
    my $align_diff = get_align_diff($rs_json);
    if ($align_diff) {
      my $old_minor_allele = $data->{'1000Genomes'}->{'minor_allele'};
      my $new_minor_allele = flip_minor_allele(
                        $old_minor_allele,
                        $data->{'v'}->{'name'},
                        $align_diff,
                        $data->{'vfs'},
                        $config->{'ma-log-fh'});
      if ($new_minor_allele) {
        $data->{'1000Genomes'}->{'org_minor_allele'} = $old_minor_allele;
        $data->{'1000Genomes'}->{'minor_allele'} = $new_minor_allele;
        my $fh = $config->{'ma-update-fh'};
        print $fh join("\t", $data->{'v'}->{'name'},
                         $old_minor_allele,
                         $new_minor_allele), "\n";
      }
    }
  }

  if ($config->{assign_ancestral_allele}) {
    assign_ancestral_alleles($data->{'vfs'}, $ancestral_alleles_utils);
  }

  # To the QC now
  # Check for the allele string matches
  # Check the map weight
  # Set the variant class
  return $data; 

}

# Assigns ancestral allele for each variation feature
# Input:
# vfs              - variation features for the variant
sub assign_ancestral_alleles {
  my ($vfs, $ancestral_alleles_utils) = @_;
  
  for my $vf (@$vfs) { 
    my $seq_name = $seq_regions_names->{$vf->{'seq_region_id'}};
    my $ancestral_allele = $ancestral_alleles_utils->assign($seq_name, $vf->{'seq_region_start'}, $vf->{'seq_region_end'});
    $vf->{'ancestral_allele'} = $ancestral_allele;
  }
}

# Flips a minor allele if there are differences in GRCh37 and GRCh38 alignment
# Input:
# old_minor_allele - allele from GRCh37
# var_name         - name of the variant
# align_diff       - alignment diff info between GRCh37 and GRCh38
# vfs              - variation features for the variant
# fh               - filehandle to log info and errors
#
# Returns:
# new_minor_allele if flip can be done, otherwise undef
sub flip_minor_allele {
  my ($old_minor_allele, $var_name,
      $align_diff, $vfs, $fh) = @_;

  return if (! $old_minor_allele);

  if ($align_diff->{'review_status'} ne 'ALIGN_DIFF') {
    print $fh join("\t", $var_name, 'INFO',
      'review_status: ' . $align_diff->{'review_status'}), "\n";
    return;
  }

  my $grch37_loc = $align_diff->{'grch37_loc'};
  my $grch38_loc = $align_diff->{'grch38_loc'};
  if (($grch37_loc != $grch38_loc) || ($grch37_loc != 1)) {
    print $fh join("\t", $var_name, 'SKIP',
      "Num of GRCh37 loc ($grch37_loc) and num of GRCh38 loc ($grch38_loc) both not equal to 1"), "\n";
    return;
  }

  my $new_minor_allele = $old_minor_allele;
  reverse_comp(\$new_minor_allele);

  # Check the new minor allele in $allele_string
  if (!@$vfs) {
    print $fh join("\t", $var_name, 'ERROR',
               "No variation features"), "\n";
    return;
  }
  my $num_vfs = @$vfs;
  my $found = 0;
  for my $vf (@$vfs) {
    my $vf_found = 0;
    my $allele_string = $vf->{'allele_string'};
    my @alleles = split("/", $allele_string);
    for my $allele (@alleles) {
        if ($new_minor_allele eq $allele) {
          $vf_found = 1;
          last;
        }
    }
    if (! $vf_found) {
      my $loc = join(':',
                     $vf->{'seq_region_id'},
                     $vf->{'seq_region_start'},
                     $vf->{'seq_region_end'});
      print $fh join("\t", $var_name, 'ERROR',
                    "new minor_allele ($new_minor_allele) not found in " .
                    "allele_string ($allele_string) of variation feature " .
                    "($loc)") . "\n";
    }
    $found += $vf_found;
  }
  if (! $found) {
    print $fh join("\t", $var_name, 'ERROR',
              "new minor_allele ($new_minor_allele) not found in any allele_string of variation_feature"), "\n";
    return;
  } elsif ($found != $num_vfs) {
    print $fh join("\t", $var_name, 'WARNING',
                    "new minor_allele ($new_minor_allele) only found in $found/$num_vfs variation_features"), "\n";
  }

  return $new_minor_allele;
}

# Get genomic coords
# Extracting information from {primary_snapshot_data} 
#                             {placements_with_alleles}
#
# Only take NC/NW locations
# TODO remove dependence on global nc_regions, nw_regions
# No QC done. QC done by looping through variants and setting class
sub get_variant_features {
  my ($rs_json) = @_;

  debug("\n>>>> get_variant_features <<<<") if ($debug);

  my @vfs;
  my $vf;
 
  my $psd = $rs_json->{'primary_snapshot_data'};
 
  my $variant_type = $psd->{'variant_type'};
 
  my ($seq_id, $seq_region_id, $seq_region);
  my ($seq_region_start, $seq_region_end,$allele_string);

  my $variation_name = 'rs' . $rs_json->{'refsnp_id'};
  foreach my $location (@{$psd->{'placements_with_allele'}}) {
    $vf = {};
    next if $location->{'seq_id'} !~ /^(NC|NW|NT)/;
    $vf->{'variant_type'} = $variant_type;
    $vf->{'variation_name'} = $variation_name;

    # Look up seq_region_id
    $seq_id = $location->{'seq_id'};
    $vf->{'seq_id'} = $seq_id;
   
    print "seq_id = ($seq_id)\n" if ($debug);
 
    if ($seq_id =~ /^NC/) {
      $seq_region_id = $nc_regions->{$seq_id};
      if (! defined $seq_region_id) {
        # Log and return;
        # Only log those where there is assembly information
        # and for the assembly of interest
        next;
      }
      $vf->{'seq_region_id'} = $seq_region_id;
    } elsif ($seq_id =~ /^NW/) {
      $seq_region = $nw_regions->{$seq_id};
      if (! defined $seq_region) {
        # Log and return;
        # Only log those were there is assembly information
        # and for the assembly of intrest
        next;
       }
      $vf->{'seq_region_id'} = $seq_region->{'seq_region_id'};
      $vf->{'asm_start'} = $seq_region->{'asm_start'};
      $vf->{'asm_end'} = $seq_region->{'asm_end'};
      $vf->{'ori'} = $seq_region->{'ori'};
    } elsif ($seq_id =~ /^NT/) {
      $seq_region = $nt_regions->{$seq_id};
      if (! defined $seq_region) {
        # Log and return;
        # Only log those were there is assembly information
        # and for the assembly of intrest
        next;
      }
      $vf->{'seq_region_id'} = $seq_region->{'seq_region_id'};
      $vf->{'asm_start'} = $seq_region->{'asm_start'};
      $vf->{'asm_end'} = $seq_region->{'asm_end'};
      $vf->{'ori'} = $seq_region->{'ori'};
    }
    
    # Check sequence type 
    my $pa = $location->{'placement_annot'};
    next if ($pa->{'seq_type'} !~ /refseq_chromosome|refseq_contig/);
    
    my $aln_opposite = $pa->{'is_aln_opposite_orientation'} || 0; 
    # Get the allele string
    next if (! @{$location->{'alleles'}});
    $vf->{'alleles'} = $location->{'alleles'};
    my ($allele_string, $position, $length_ref, $allele_errors) = get_allele_info($location->{'alleles'});
    $vf->{'allele_string'} = $allele_string;
    $vf->{'seq_region_start'} = $position + 1;
    $vf->{'seq_region_end'} = $position + $length_ref;
    $vf->{'seq_region_strand'} = 1; 
    
    if (defined $vf->{'asm_start'}) {
      if ($vf->{'ori'} != -1) {
        $vf->{'seq_region_start'} = $position + $vf->{'asm_start'};
        $vf->{'seq_region_end'} = $position - 1 + $length_ref + $vf->{'asm_start'};
      } else {
        $vf->{'seq_region_strand'} = -1;
        # Position is zero-based
        $vf->{'seq_region_end'} = $vf->{'asm_end'} - $position;
        $vf->{'seq_region_start'} = $vf->{'seq_region_end'} - $length_ref + 1;
      }
    }
    $vf->{'position'} = $position;
    $vf->{'aln_opposite'} = $aln_opposite;
    $vf->{'allele_errors'} = $allele_errors;
    push @vfs, $vf;
  }
  return \@vfs;
}

sub get_allele_info {
  my ($alleles_ref) = @_;
  my ($allele_str, $position, $ref_length);
  my %allele_errors;
  my @allele_string;

  my $ref_allele;

  # Loop through all the alleles
  # Is a ref found
  # Is deleted_position_same
  my $ref_found = 0;
  my $del_seq_differ = 0;
  my $pos_differ = 0;

  my $first_del_seq = $alleles_ref->[0]->{'allele'}->{'spdi'}->{'deleted_sequence'};
  my $first_pos     = $alleles_ref->[0]->{'allele'}->{'spdi'}->{'position'};

  for (my $n = 0; $n < @{$alleles_ref}; $n++) {
    # Get the placement annotated allele
    my $paa = $alleles_ref->[$n];
    my $spdi = $paa->{'allele'}->{'spdi'};
    if ($spdi->{'deleted_sequence'} eq $spdi->{'inserted_sequence'}) {
      $ref_allele = $spdi->{'deleted_sequence'} || "-";
      $ref_found = 1;
    } else {
      push @allele_string, $spdi->{'inserted_sequence'} || "-";
    }
    if ($spdi->{'deleted_sequence'} ne $first_del_seq) {
      $del_seq_differ = 1;
      $allele_errors{2} = 1;
    }
    if ($spdi->{'position'} ne $first_pos) {
      $pos_differ = 1;
      $allele_errors{3} = 1;
    }
  }
  if (! $ref_found) {
      $allele_errors{1} = 1;
      # Then set to the first deleted sequence
      # Regardless if there is a del_seq_differ
      # or a pos_differ
      $ref_allele = $first_del_seq || "-";
  }
  unshift @allele_string, $ref_allele;
  my $vf_allele_string;
  if ($del_seq_differ || $pos_differ) {
    $vf_allele_string = 'dbSNP_variant';
  } elsif (@allele_string == 1) {
    $vf_allele_string = 'dbSNP_novariation';
    $allele_errors{4} = 1
  } else {
    $vf_allele_string = join("/", @allele_string);
    if (length($vf_allele_string) > 50000) {
      $vf_allele_string = 'dbSNP_variant';
      $allele_errors{5} = 1;
    } 
  }
  my $length_ref_allele = 0;
  if ($ref_allele && $ref_allele ne '-') {
      $length_ref_allele = length($ref_allele);
  }
  #print "vf_allele_string = ($vf_allele_string)\n";
  #print "ref_allele = ($ref_allele)\n";  
  #print "length = ", $length_ref_allele, "\n";
  #print "position = ($first_pos)\n";
  #print "Allele_errors:", Dumper(\%allele_errors);
  my @allele_errors_a = sort {$a <=> $b} keys %allele_errors;
  return ($vf_allele_string, $first_pos, $length_ref_allele, \@allele_errors_a);
}

sub get_merges {
  my ($rs_json) = @_;

  # Get the dbsnp1_merges
  my $dbsnp1_merges = $rs_json->{'dbsnp1_merges'};
  my @merged_refsnps;

  # Format of merges
  # dbsnp1_merge_event {
  #   merged_rsid (string):
  #   revision (string):
  #   merge_date (string):
  # }
  for my $dbsnp1_merge_event (@$dbsnp1_merges) {
    push @merged_refsnps,  "rs" . $dbsnp1_merge_event->{'merged_rsid'};
  }
  return \@merged_refsnps;
}

sub get_dbsnp2_merges {
  my ($rs_json) = @_;
  my %dbsnp2_merges;
  my $present_obs_movements = $rs_json->{'present_obs_movements'};
  for my $pom (@$present_obs_movements) {
    my $allele_in_cur_release = $pom->{'allele_in_cur_release'};
    if ($allele_in_cur_release->{'deleted_sequence'} eq
      $allele_in_cur_release->{'inserted_sequence'}) {
      next;
    }
    my $prev_rel_rsids = $pom->{'previous_release'}->{'rsids'};
    for my $prev_rel_rsid (@$prev_rel_rsids) {
      # Only store if not equal to refsnp_id
      if ($prev_rel_rsid != $rs_json->{'refsnp_id'}) {
        $dbsnp2_merges{$prev_rel_rsid}++;
      }
    }
  }
  return [ map { 'rs'. $_} keys %dbsnp2_merges];
}

sub get_hgvs {
  my ($config, $rs_json) = @_;

  debug("\n>>>> get_hgvs <<<<") if ($debug);

  my $psd = $rs_json->{'primary_snapshot_data'};

  my %hgvs;
  my $info;
  foreach my $loc (@{$psd->{'placements_with_allele'}}) {
    my $la = $loc->{'placement_annot'};
    if ($la->{'seq_type'} =~ /^(refseq_mrna|refseq_prot)$/) {
      for my $allele (@{$loc->{'alleles'}}) {
        next if ($allele->{'hgvs'} !~ /^(NP|NM)/);
        next if ($allele->{'hgvs'} =~ /^(NP|NM).+=$/);
        if (length($allele->{'hgvs'}) >= 255) {
          $info = join(";",
                 'length_HGVS=' . length($allele->{'hgvs'}),
                 'HGVS_start=' . substr($allele->{'hgvs'}, 0, 50));
          log_errors($config, 'rs' . $rs_json->{'refsnp_id'}, 'HGVS_length_ge_255', $info);
          next;
        }
        $hgvs{$allele->{'hgvs'}}++;
      }
    }
  }
  my @hgvs = keys %hgvs;
  return \@hgvs;
}

#citations (Array[integer]):
# Set of Pubmed IDs (PMIDs) for this RefSNP or its supporting submissions ,
sub get_citations {
  my ($rs_json) = @_;

  debug("\n>>>> get_citations <<<<") if ($debug);

  return $rs_json->{'citations'};
}

sub get_support {
  my ($rs_json) = @_;

  debug("\n>>>> get_support <<<<") if ($debug);

  my %subsnp_support;
  my %freq_support;

  # Extract subsnp and frequency submitters
  foreach my $support (@{$rs_json->{'primary_snapshot_data'}->{'support'}}) {
    if ($support->{'id'}->{'type'} eq 'subsnp') {
      $subsnp_support{$support->{'submitter_handle'}}++;
    } elsif ($support->{'id'}->{'type'} eq 'frequency') {
      $freq_support{$support->{'submitter_handle'}}++;
    }
  }
  return (\%subsnp_support, \%freq_support);
}

# QC the refSNP
sub qc_refsnp {
  my ($config, $lu_info, $rs_data, $fai_index) = @_;
  $rs_data->{'evidence_attribs'} = get_evidence($lu_info, $rs_data); 
  $rs_data->{'map_weight'} = get_map_weight($ref_regions, $rs_data) if ($add_map_weight);
  $rs_data->{'num_par'} = get_vf_par($rs_data, $config->{'assembly'},
                                     $lu_info->{'chrY_seq_region_id'});

  # Warn if map_weight = num_par
  if ($rs_data->{'num_par'} && $rs_data->{'num_par'} == $rs_data->{'map_weight'}) {
      my $info = join(";",
                'map_weight=' . $rs_data->{'map_weight'},
                'num_par=' . $rs_data->{'num_par'});
      log_errors($config, $rs_data->{'v'}->{'name'},
          'num PAR matches map weight',
          $info);
  }
  $rs_data->{'variant_fails'} = get_variant_vf_fails($rs_data, $config->{'ref_check'}, $fai_index, $ref_regions);
  $rs_data->{'display'} = set_display($rs_data);
} 

# TODO remove check the lu for the evidence exists
# Should be part of the initialisation check
sub get_evidence {
  my ($lu_info, $rs_data) = @_;
  my @evidence;

  my $added_freq_evidence = 0;
  if (%{$rs_data->{'freq_support'}}) {
    push @evidence, $lu_info->{'evidence_ids'}->{'Frequency'};
    $added_freq_evidence = 1;
  } 

  if (@{$rs_data->{'citations'}}) {
    push @evidence, $lu_info->{'evidence_ids'}->{'Cited'};
  }

  my $freq_evidence = 0;

  my $subsnp_support = $rs_data->{'snp_support'};

  # Determine 1000Genomes evidence: subsnp support or frequency support
  # submitter_handle has different capitalisation for
  # subsnp suppport (1000GENOMES)
  # frequency support (1000Genomes)
  if (defined $subsnp_support) {
    if (defined $subsnp_support->{'1000GENOMES'}) {
      if (defined $lu_info->{'evidence_ids'}->{'1000Genomes'}) {
        push @evidence, $lu_info->{'evidence_ids'}->{'1000Genomes'};
        $freq_evidence++;
      }
    } elsif (defined $rs_data->{'freq_support'}->{'1000Genomes'}) {
        push @evidence, $lu_info->{'evidence_ids'}->{'1000Genomes'};
    }

    # Determine ESP evidence: subsnp support or frequency support
    # submitter_handle is different for
    # subsnp suppport (NHLBI-ESP)
    # and frequency support (GoESP)
    if (defined $subsnp_support->{'NHLBI-ESP'}) {
      if (defined $lu_info->{'evidence_ids'}->{'ESP'}) {
        push @evidence, $lu_info->{'evidence_ids'}->{'ESP'};
        $freq_evidence++;
      }
    } elsif (defined $rs_data->{'freq_support'}->{'GoESP'}) {
        push @evidence, $lu_info->{'evidence_ids'}->{'ESP'};
    }

    if (defined $subsnp_support->{'EVA_EXAC'}) {
      if (defined $lu_info->{'evidence_ids'}->{'ExAC'}) {
        push @evidence, $lu_info->{'evidence_ids'}->{'ExAC'};
        $freq_evidence++;
      }
    }

    if (defined $subsnp_support->{'TOPMED'}) {
      if (defined $lu_info->{'evidence_ids'}->{'TOPMed'}) {
        push @evidence, $lu_info->{'evidence_ids'}->{'TOPMed'};
        $freq_evidence++;
      }
    }
    if (defined $subsnp_support->{'GNOMAD'}) {
      if (defined $lu_info->{'evidence_ids'}->{'gnomAD'}) {
        push @evidence, $lu_info->{'evidence_ids'}->{'gnomAD'};
        $freq_evidence++;
      }
    }

    if (! $added_freq_evidence && ($freq_evidence)) {
      push @evidence, $lu_info->{'evidence_ids'}->{'Frequency'};
    }
  }

  return \@evidence;
}

sub get_map_weight {
  my ($regions, $rs_data) = @_;
  my $map_weight=0;
  my %seen_vf;

  # Loop through the vfs and count the regions that are top-level
  # Only count the unique seq_region_id:seq_region_start:seq_region_end
  for my $vf (@{$rs_data->{'vfs'}}) {
    if (defined $regions->{$vf->{'seq_region_id'}}) {
      # The import can contain duplicate VFs based on
      # seq_region_id, seq_region_start, seq_region_end
      # These can arise when NTs map to same position as NCs
      # The duplicates are excluded in map_weight calculations
      my $loc = join(':',
                     $vf->{'seq_region_id'},
                     $vf->{'seq_region_start'},
                     $vf->{'seq_region_end'});
      next if (defined $seen_vf{$loc});
      $seen_vf{$loc} = 1;
      $map_weight++;
    }
  }
  return $map_weight;
}

sub get_vf_par {
  my ($rs_data, $assembly, $Y_seq_region_id) = @_;

  my $num_par = 0;

  my ($Y_PAR1_seq_region_start, $Y_PAR1_seq_region_end);
  my ($Y_PAR2_seq_region_start, $Y_PAR2_seq_region_end);

  if ($assembly eq 'GRCh38') {
    $Y_PAR1_seq_region_start = 10001 ;
    $Y_PAR1_seq_region_end = 2781479;
    $Y_PAR2_seq_region_start = 56887903;
    $Y_PAR2_seq_region_end = 57217415;
  } elsif ($assembly eq 'GRCh37') {
    $Y_PAR1_seq_region_start = 10001;
    $Y_PAR1_seq_region_end = 2649520;
    $Y_PAR2_seq_region_start = 59034050;
    $Y_PAR2_seq_region_end = 59363566;
  }

  # Loop the vfs and count the regions that are PAR
  for my $vf (@{$rs_data->{'vfs'}}) {
    $vf->{'par'} = 0;
    if ($vf->{'seq_region_id'} == $Y_seq_region_id) {
      if ($vf->{'seq_region_start'} >= $Y_PAR1_seq_region_start
          && $vf->{'seq_region_end'}   <= $Y_PAR1_seq_region_end
          && $vf->{'seq_region_start'} <= $Y_PAR1_seq_region_end
          && $vf->{'seq_region_end'} >= $Y_PAR1_seq_region_start) {
        $num_par++;
        $vf->{'par'} = 1;
      } elsif ($vf->{'seq_region_start'} >= $Y_PAR2_seq_region_start
          && $vf->{'seq_region_end'} <=   $Y_PAR2_seq_region_end
          && $vf->{'seq_region_start'} <= $Y_PAR2_seq_region_end
          && $vf->{'seq_region_end'} >= $Y_PAR2_seq_region_start) {
        $num_par++;
        $vf->{'par'} = 1;
      } 
    }
  }
  return $num_par;
}

# Based on QCUtils::run_vf_checks
sub get_variant_vf_fails {
  my ($rs_data, $ref_check, $fai_index, $regions) = @_;

  # The variant fails
  my %v_fails;

  # If there are no VF, there are no variant fails based on VF
  if (!@{$rs_data->{'vfs'}}) {
    return {%v_fails};
  }

  # Type 19 - fail variant if >1 mapping seen
  if (defined $rs_data->{'map_weight'} && $rs_data->{'map_weight'} > 1) {
      if (($rs_data->{'map_weight'} - $rs_data->{'num_par'}) > 1) {
        $v_fails{19} = 1;
      }
  }

  # For each variant feature add failure. Keep
  # track of the following so can fail at variant level.
  # Type 13 - Alleles contain non-nucleotide characters
  # Type 14 - Alleles contain ambiguity codes
  # Type 20 - Variant at first base in sequence
  # Type 21 - Reference allele does not match the bases at this genome location
  my %fail_count;
  my @error_codes = (13, 14, 20);
  foreach my $code (@error_codes) {
    $fail_count{$code} = 0;
  }
  my $num_genomic_fail = 0;
  my $num_allele_error = 0;
  my $num_no_variation = 0;

  for my $vf (@{$rs_data->{'vfs'}}) {
    # If had allele_errors that could not be resolved
    # allele string set to dbSNP_, do not QC the variant feature
    if ($vf->{'allele_string'} =~ /^dbSNP_/) {
      $vf->{'fails'} = {};
      $num_allele_error++;
      if ($vf->{'allele_string'} eq 'dbSNP_novariation') {
        $num_no_variation++;
      }
      next;
    }
    $vf->{'fails'} = get_vf_fails($vf, $ref_check, $fai_index);

    foreach my $code (@error_codes) {
      if (defined $vf->{'fails'}->{$code}) {
        $fail_count{$code}++;
      }
      # Is this a genomic location?
      if (defined $regions->{$vf->{'seq_region_id'}}) {
        if (defined $vf->{'fails'}->{21}) {
          $num_genomic_fail++;
        }
      }
    }
  }
  my $num_vf = scalar(@{$rs_data->{'vfs'}});
  foreach my $code (@error_codes) {
    if ($fail_count{$code} == $num_vf) {
      $v_fails{$code} = 1;
    }
  }
  # All VF have allele_errors
  if ($num_vf == $num_allele_error) {
    if ($num_vf == $num_no_variation) {
      $v_fails{4} = 1;
    } else {
      $v_fails{22} = 1;
    }
  }
  if (defined $rs_data->{'map_weight'} &&
      ($rs_data->{'map_weight'} == 1) &&
      ($num_genomic_fail == 1)) {
      $v_fails{21} = 1;
  }

  return {%v_fails};
}

sub get_vf_fails {
  my ($vf, $ref_check, $fai_index) = @_;

  my %vf_fails;

  # Type 20 Variant at first base in sequence
  # Unlikely a variant at the first base of a reference sequence
  # is a good call
  if ($vf->{'seq_region_start'} == 1 ) {
    $vf_fails{20} = 1;
  }

  my $allele_string = $vf->{'allele_string'};
  # Type 13  Alleles contain non-nucleotide characters
  # If illegal characters stop checking this VF
  my $illegal_alleles = check_illegal_characters($allele_string);
  if (defined $illegal_alleles->[0]) {
    $vf_fails{13} = 1;
    return {%vf_fails};
  }

  # Type 14 -  Alleles contain ambiguity codes
  # Currently do not resolve ambiguities.
  # Marked as a vf_failure and do not do reference check
  my $is_ambiguous = check_for_ambiguous_alleles($allele_string);
  if (defined $is_ambiguous) {
    $vf_fails{14} = 1;
    return {%vf_fails};
  }

  # Note: No flipping done as assumed to be on forward strand
  #       for dbSNP
  #
  # Previously Reference match checks and re-ordering only run for variants with 1 genomic location
  # Patches were still checked
  # So check all variant features now

  # Do the reference checks
  # Extract reference sequence to run ref checks
  if ($ref_check) {
      # 21 Reference allele does not match the bases at this genome location
      my $ref_seq = get_reference_base_hts($vf, $fai_index);
      # Keeping a record of map check on variation feature
      # Previously a map failure recorded by variant
      unless ($ref_seq) {
        ## don't check further if obvious coordinate error giving no sequence
        $vf_fails{21} = 1;
        $vf->{'ref_correct'} = 'no';
        return {%vf_fails};
      }
      my @alleles = split/\//, $vf->{allele_string};
      # The empty allele is represented by "-"
      # get_reference_base_hts return "-"
      # No specific checking for "-" needed
      if ($alleles[0] ne $ref_seq) {
        $vf_fails{21} = 1;
        $vf->{'ref_correct'} = 'no';
      }
  }
  return {%vf_fails};
}

sub get_reference_base_hts {
  my ($vf, $fai_index) = @_;

  #print "Getting reference base for a VF\n";
  my $ref_seq; 
  if ( $vf->{'seq_region_end'} + 1 == $vf->{'seq_region_start'}) { 
    # Convention for insertions to reference
    $ref_seq = "-";
  } elsif ($vf->{'seq_region_end'} < $vf->{'seq_region_start'}) {
    # Coordinate error;
    warn "Incorrect coords $vf->{'seq_region_end'} - $vf->{'seq_region_start'} for $vf->{'name'} \n";
    return;
  } else {
    my $seq_name = $seq_regions_names->{$vf->{'seq_region_id'}};
    my $location = $seq_name . ":" . $vf->{'seq_region_start'} . "-" . $vf->{'seq_region_end'};
    $ref_seq = $fai_index->get_sequence_no_length($location);
    # If strand is "-" need to reverse complement it
    if ($vf->{'seq_region_strand'} == -1) {
      reverse_comp(\$ref_seq);
    }
  }
  return $ref_seq;
}

# Set the variants display
# Currently the display for the variant 
# and all its variation_features is the same
sub set_display {
  my ($rs_data) = @_;

  # Assume display = 1;
  my $display = 1;
  
  # If there are citations, then display
  # Even if no mappings?
  if (@{$rs_data->{'citations'}}) {
      return $display;
  }

  # If there are no mappings, then do not display
  if (! (@{$rs_data->{'vfs'}})) {
    $display = 0;
    return $display;
  }

  # If there are failures, then do not display
  if (%{$rs_data->{'variant_fails'}}) {
    $display = 0;
    return $display;
  }
  return $display;
}


# Stores the data in the database
sub import_refsnp {
  my ($dbh, $rs_data, $config) = @_;

  debug("\n>>>> importing record to database <<<<") if ($debug);

  my $variation_id = import_variation($dbh, $rs_data);
  print "Imported $variation_id\n" if ($debug);

  die("no variation id") if ($variation_id < 1);

  import_batch($dbh, $config->{'batch_id'}, $variation_id, $rs_data->{'variant_type'});

  import_variation_feature($dbh, $variation_id, $rs_data, $config, $lu_info->{'failed_set_id'});
  import_merges($dbh, $variation_id,
                $rs_data->{'merges'},
                $config->{'source_id'}->{'Archive dbSNP'});

  import_merges($dbh, $variation_id,
                $rs_data->{'dbsnp2_merges'},
                $config->{'source_id'}->{'Former dbSNP'});

  import_hgvs($dbh, $variation_id,
              $rs_data->{'hgvs'},
              $config->{'source_id'}->{'dbSNP HGVS'});

  import_citations($dbh, $variation_id, $rs_data->{'name'}, $rs_data->{'citations'});

  if (! (@{$rs_data->{'vfs'}})) {
    add_unmapped_variant($dbh, $variation_id);
  }
 
  if (%{$rs_data->{'variant_fails'}}) {
    add_variant_fails($dbh, $variation_id, $rs_data->{'variant_fails'});
  } 
  #add_failed_variation();
  #add_publication();
}

sub import_batch {
  my ($dbh, $batch_id, $variation_id, $variant_type) = @_;

  $dbh->do(qq{INSERT INTO batch_variation
             (batch_id, variation_id, variant_type)
             VALUES
             (?, ?, ?)},
             undef,
             $batch_id, $variation_id, $variant_type);
}

# Add the variation record
# TODO:
# - Add freq and evidence setting to function
# - Precision of 1000Genome frequencey
sub import_variation {
  my ($dbh, $data) = @_;

  my $v = $data->{'v'};
  my $freq = $data->{'1000Genomes'};
 
  # Determine the values for maf, minor_allele, minor_allele_count
  my ($maf, $minor_allele, $minor_allele_count);
  if (defined $freq) {
      $minor_allele = $freq->{'minor_allele'} || '-';
      $maf = format_frequency($freq->{'MAF'});
      $minor_allele_count = $freq->{'minor_allele_count'};
  }

  # Get the evidence
  my $evidence_attribs_str;
  if (@{$data->{'evidence_attribs'}}) {
    $evidence_attribs_str = join(",", @{$data->{'evidence_attribs'}});
  }

  $dbh->do(qq{INSERT INTO variation (name, source_id,
                             minor_allele, minor_allele_freq, minor_allele_count,
                             evidence_attribs, display,
                             class_attrib_id) VALUES
                    (?, ?, 
                     ?, ?, ?,
                     ?, ?,
                     ?)},
                    undef,
                   $v->{'name'}, $v->{'source_id'}, 
                   $minor_allele, $maf, $minor_allele_count,
                   $evidence_attribs_str, $data->{'display'},
                   $v->{'class_attrib_id'});
  my $db_variation_id = $dbh->last_insert_id(undef, undef, qw(variation variation_id)) or die("no insert id for variation");
  if (! $db_variation_id) {
    die("Unable to get variation_id for $v->{'name'}: $!");
  }
  $v->{'variation_id'} = $db_variation_id;
  return $db_variation_id;
}

# Add the variation features
# TODO:
# - Add freq and evidence setting to function
# - Source not be hardcoded
# - Precision of 1000Genomes frequency
# - Remove addition of non-schema columns
sub import_variation_feature {
  my ($dbh, $variation_id, $data, $config, $failed_set_id) = @_;

  die("no variation id") if (! $variation_id);

  my $vfs = $data->{'vfs'};
  my $freq = $data->{'1000Genomes'};
  my $source_id = $data->{'v'}->{'source_id'};

  my ($maf, $minor_allele, $minor_allele_count);
  if (defined $freq) {
      $minor_allele = $freq->{'minor_allele'} || '-';
      $maf = format_frequency($freq->{'MAF'});
      $minor_allele_count = $freq->{'minor_allele_count'};
  }
  
  my $evidence_attribs_str;
  if (@{$data->{'evidence_attribs'}}) {
    $evidence_attribs_str = join(",", @{$data->{'evidence_attribs'}});
  }

  my $map_weight = 0;
  if ($data->{'map_weight'}) {
    $map_weight = $data->{'map_weight'};
  }

  # If there are failed variants, add the failed variant set
  # As for display, set on variant level. To do need to set
  # on variation_feature level
  # variation_set_id defaults to empty string not NULL;
  my $sets = '';
 
  if (%{$data->{'variant_fails'}}) {
    $sets = $failed_set_id;
  }
 
  # TODO standardise allele_error(s) name
  # TODO review storing seq_id 
  #'seq_id' => 'NC_000008.11',

  # 'allele_errors' => '',
  # 'seq_id' => 'NC_000021.9',
  # 'position' => 18862898,
  # 'seq_region_id' => '131543',
  # 'aln_opposite' => 0,
  # 'seq_region_start' => 18862899,
  # 'variant_type' => 'delins',
  # 'variation_name' => 'rs66465154',
  # 'seq_region_end' => 18862899,
  # 'allele_string' => 'G/AA'

  my $sth=$dbh->prepare(qq[INSERT INTO variation_feature
                           (variation_name, map_weight, 
                            seq_region_id, seq_region_start, seq_region_end,
                            seq_region_strand, 
                            variation_id, allele_string, ancestral_allele,
                            source_id, variation_set_id, class_attrib_id,
                            minor_allele, minor_allele_freq, minor_allele_count,
                            evidence_attribs, display
                            )
                          VALUES (
                            ?, ?,
                            ?, ?, ?,
                            ?,
                            ?, ?, ?,
                            ?, ?, ?,
                            ?, ?, ?,
                            ?, ?)]);
  my %seen_vf;

  for my $vf (@$vfs) {
    # The import can contain duplicate VFs based on
    # seq_region_id, seq_region_start, seq_region_end
    # These can arise when NTs map to same position as NCs
    # The duplicates are not imported into the database
    my $loc = join(':',
                    $vf->{'seq_region_id'},
                    $vf->{'seq_region_start'},
                    $vf->{'seq_region_end'});

    if (defined $seen_vf{$loc}) {
      # Keep a record of the location and alleles in the JSON file
      import_placement_allele($dbh, $variation_id, $vf->{'alleles'}, $vf->{'variation_name'});

      # Log only the location not the alleles in the error file
      my $info = join(";",
                'seq_region_id=' . $vf->{'seq_region_id'},
                'seq_region_start=' . $vf->{'seq_region_start'},
                'seq_region_end=' . $vf->{'seq_region_end'},
                'seq_id=' . $vf->{'seq_id'},
                'position=' . $vf->{'position'},
                'strand=' . $vf->{'seq_region_strand'});
      log_errors($config, $vf->{'variation_name'},
                 'duplicate_variation_feature',
                  $info);
      next;
    }
    $seen_vf{$loc} = 1;

    $sth->execute($vf->{'variation_name'}, $map_weight,
                  $vf->{'seq_region_id'}, $vf->{'seq_region_start'}, $vf->{'seq_region_end'},
                  $vf->{'seq_region_strand'},
                  $variation_id, $vf->{'allele_string'}, $vf->{'ancestral_allele'},
                  $source_id, $sets, $vf->{'class_attrib_id'},
                  $minor_allele, $maf, $minor_allele_count,
                  $evidence_attribs_str, $data->{'display'});

    # If there are allele_errors, add these to the failed_variation_feature
    # Need to get the last_insert_id in this case

    if (@{$vf->{'allele_errors'}} ||
           %{$vf->{'fails'}}) {
      my $db_vf_id = $dbh->last_insert_id(undef, undef, qw(variation_feature variation_feature_id)) or die("no insert id for variation_feature");
      if (! $db_vf_id) {
        die("Unable to get variation_feature_id for $vf->{'variation_name'}: $!");
      }
      if (@{$vf->{'allele_errors'}}) {
        import_failed_alleles($dbh, $db_vf_id, $vf->{'allele_errors'}, $vf->{'allele_string'});
      }
      if (%{$vf->{'fails'}}) {
        import_failed_variation_features($dbh, $db_vf_id, $vf->{'fails'});
      }
    }

    # insert the placement allele
    import_placement_allele($dbh, $variation_id, $vf->{'alleles'}, $vf->{'variation_name'});
    
    # Log those with aln_opposite but seq_region_strand +1
    if ($vf->{'aln_opposite'} && ($vf->{'seq_region_strand'} == 1)) {
      my $info = join(";",
                'seq_id=' . $vf->{'seq_id'},
                'position=' . $vf->{'position'},
                'strand=' . $vf->{'seq_region_strand'});
      log_errors($config, $vf->{'variation_name'},
          'aln_opposite_orientation',
          $info);
    }
  }
}

sub import_placement_allele {
  my ($dbh, $variation_id, $alleles, $var_name) = @_;

  debug("import_placement_allele") if ($debug);
 
  die("no variation id") if (! $variation_id);

  # 'alleles' => [
  #                {
  #                  'hgvs' => 'NC_000015.10:g.22999284G=',
  #                  'allele' => {
  #                                'spdi' => {
  #                                  'deleted_sequence' => 'G',
  #                                  'seq_id' => 'NC_000015.10',
  #                                  'position' => 22999283,
  #                                  'inserted_sequence' => 'G'
  #                                }
  #                              }
  #                },
  #                {
  #                  'hgvs' => 'NC_000015.10:g.22999284G>A',
  #                  'allele' => {
  #                               'spdi' => {
  #                                 'deleted_sequence' => 'G',
  #                                 'seq_id' => 'NC_000015.10',
  #                                 'position' => 22999283,
  #                                 'inserted_sequence' => 'A'
  #                                }
  #                              }
  #                 }
  #              ],
  my $sth=$dbh->prepare(qq[INSERT INTO placement_allele
                           (variation_id,
                            seq_id, position, deleted_sequence, inserted_sequence,
                            hgvs)
                          VALUES (
                            ?,
                            ?, ?, ?, ?,
                            ?)]);
  for my $allele (@$alleles) {
    my $spdi = $allele->{'allele'}->{'spdi'};
    #print "seq_id = $spdi->{'seq_id'}\n";
    #print "position = $spdi->{'position'}\n";
    #print "deleted_sequence = $spdi->{'deleted_sequence'}\n";
    #print "insereted sequence = $spdi->{'inserted_sequence'}\n";
    #print "hgvs = $allele->{'hgvs'}\n";

    # The placement_allele table defines deleted_sequence, inserted_sequence
    # and hgvs to varchar(255)
    # Only the first 255 characters of these fields are stored and truncations flagged
    # TODO - alter the schema for placement allele
    #
    check_spdi_lengths($spdi, $variation_id, $var_name);

    my $pa_hgvs = $allele->{'hgvs'};
    if (length($pa_hgvs) > 255) {
      $pa_hgvs = substr($pa_hgvs,0,255);
      my $msg = 'truncated_pa_hgvs';
      my $info = join(";", 'variant_id='. $variation_id,
                           'hgvs=' . $pa_hgvs);
      log_errors($config, $var_name, $msg, $info);
    }

    $sth->execute($variation_id,
                  $spdi->{'seq_id'},
                  $spdi->{'position'},
                  $spdi->{'deleted_sequence'},
                  $spdi->{'inserted_sequence'},
                  $pa_hgvs);
  }
}


sub check_spdi_lengths {
  my ($spdi, $var_id, $var_name) = @_;

  my $msg = 'spdi_seq_gt_255';

  my @fields = ('deleted_sequence', 'inserted_sequence');
  for my $field (@fields) {
    if (length($spdi->{$field}) > 255) {
      $spdi->{$field} = substr($spdi->{$field}, 0,255);
      my $info = join(";", 'variant_id='. $var_id,
                         $field . '=' . $spdi->{$field});
      log_errors($config, $var_name, $msg, $info);
    }
  }
}


# Stores the SPDI allele errors - for local use only
sub import_failed_alleles {
  my ($dbh, $vf_id, $allele_errors, $allele_string) = @_;

  my $sth=$dbh->prepare(qq[INSERT INTO failed_variation_feature_spdi
                           (variation_feature_id, spdi_failed_description_id)
                          VALUES (
                            ?, ?)]);
  my $sth_vf=$dbh->prepare(qq[INSERT INTO failed_variation_feature
                           (variation_feature_id, failed_description_id)
                          VALUES (
                           ?, ?)]);
  for my $ae (@$allele_errors) {
    $sth->execute($vf_id, $ae);
  }
  if ($allele_string eq 'dbSNP_variant') {
    $sth_vf->execute($vf_id, 22);
  } elsif ($allele_string eq 'dbSNP_novariation') {
    $sth_vf->execute($vf_id, 4);
  }
}

sub import_failed_variation_features {
  my ($dbh, $vf_id, $errors) = @_;

  my $sth=$dbh->prepare(qq[INSERT INTO failed_variation_feature
                           (variation_feature_id, failed_description_id)
                          VALUES (
                            ?, ?)]);
  for my $error (sort {$a <=> $b} keys %$errors) {
    $sth->execute($vf_id, $error);
  }
}

# Add the merges into variation_synonym
# The index has changed on variation_synonym, allows
# duplicate synonyms for a given source
sub import_merges {
  my ($dbh, $variation_id, $merges, $source_id) = @_;

  debug("import_merges") if ($debug);
 
  die("no variation id") if (! $variation_id);

  #'merges' => [
  #                      'rs17850737',
  #                      'rs52818902',
  #                      'rs386571803'
  #            ],
  my $sth=$dbh->prepare(qq[INSERT INTO variation_synonym
                           (variation_id,
                            source_id,
                            name) 
                          VALUES (
                            ?,
                            ?,
                            ?)]);
  for my $name (@$merges) {
    $sth->execute($variation_id, $source_id, $name);
  }
}

# Add the HGVS to the database
# The index has changed on variation_synonym, allows
# duplicate synonyms for a given source
sub import_hgvs {
  my ($dbh, $variation_id, $hgvs, $source_id) = @_;
  
  debug("import hgvs") if ($debug);
  
  die("no variation id") if (! $variation_id);

  # [
  #    'NP_000228.1:p.Asn318Ser',
  #    'NM_000237.2:c.953A>G'
  # ]
  #
  my $sth = $dbh->prepare(qq[INSERT INTO variation_synonym
                          (variation_id,
                           source_id,
                           name)
                          VALUES (
                           ?,
                           ?,
                           ?)]);
  for my $name (@$hgvs) {
    $sth->execute($variation_id, $source_id, $name);
  }
}

#
# Add the citations to tmp_variation_citation
# Later this will be moved to publication, variation_citation
sub import_citations {
  my ($dbh, $variation_id, $refsnp_id, $citations) = @_;
  
  debug("import_citations") if ($debug);
  
  die("no variation id") if (! $variation_id);

  my $inserts;
  # [
  #    'NP_000228.1:p.Asn318Ser',
  #    'NM_000237.2:c.953A>G'
  #  ]
 
  # Note doing an INSERT IGNORE just to get past the problem of rs HGVS 
  # This needs to be removed
  my $sth = $dbh->prepare(qq[INSERT INTO tmp_variation_citation
                          (variation_id,
                           pmid )
                          VALUES (
                           ?,?)]);
  for my $citation (@$citations) {
    $sth->execute($variation_id, $citation);
   }
}

sub add_unmapped_variant {
  my ($dbh, $variation_id) = @_;
 
  debug("add_unmapped_variant") if ($debug);

  die("no variation id") if (! $variation_id);

  my $sth = $dbh->prepare(qq[INSERT INTO failed_variation
                             (variation_id, failed_description_id)
                             VALUES (?, ?)
                            ]);
  
  $sth->execute($variation_id, 5);
}

sub add_variant_fails {
  my ($dbh, $variation_id, $fails) = @_;
 
  debug("add_variant_fails") if ($debug);

  die("no variation id") if (! $variation_id);

  my $sth = $dbh->prepare(qq[INSERT INTO failed_variation
                             (variation_id, failed_description_id)
                             VALUES (?, ?)
                            ]);
  for my $fail_id (keys %{$fails}) {
    $sth->execute($variation_id, $fail_id) or die("Error inserting variation fails info");
  }
}

sub get_variation_id {

  # Sub to be removed - so temp declaration of $run_info to pass compile
  my $run_info;
  my $variant_id = $run_info->{'variant_seq'}->{'next_id'};

  $run_info->{'variant_seq'}->{'next_id'}++;

  my $start = $run_info->{'variant_seq'}->{'start'};
  my $end   = $run_info->{'variant_seq'}->{'end'};
  if (($variant_id > $end) || ($variant_id < $start)) {
    die("variant_id ($variant_id) is outside of sequence range ($start - $end)");
  }
  return $variant_id;
}
  
sub report_summary{

  my ($config, $num_lines) = @_;

  my $species = $config->{'species'};
  my $input_file = $config->{'input_file'};
  my $report_file = $config->{'rpt_file'};
  my $batch_id = $config->{'batch_id'};

  debug("\n>>>> report_summary <<<<") if ($debug);

  open my $report, '>', $report_file or die("Failed to open $report_file to write: $!");
  print $report "Report date:\t", strftime('%Y-%m-%d %H:%M:%S', localtime()), "\n";
  print $report "Input file:\t$input_file\n";
  print $report "Lines processed:\t$num_lines\n";
  print $report "Batch id:\t$batch_id\n";

  if (!$config->{'db_load'}) {
    print $report "Configured for no db load\n";
    close($report);
    return;
  }

  my $dbh = $config->{'dbh_var'};

  my $variation_count = count_rows_batch($dbh, 'variation', 'variation_id', $batch_id);
  print $report "Number variation:\t$variation_count\n";
  my $vf_count = count_rows_batch($dbh, 'variation_feature', 'variation_feature_id', $batch_id);
  print $report "Number variation_features:\t$vf_count\n";

  # Get failure rate
  my $failed_variation_count = count_rows_batch($dbh, 'failed_variation', 'variation_id', $batch_id);

  my $var_fail_rate = sprintf("%.4g", 100 * ($failed_variation_count/$variation_count));
  print $report "Variation failure rate:\t$var_fail_rate\% [$failed_variation_count / $variation_count]\n";

  # Get fails by type
  my $sql = qq{
        SELECT fd.description, count(*)
        FROM failed_description fd, failed_variation fv, batch_variation bv
        WHERE fd.failed_description_id = fv.failed_description_id
        AND fv.variation_id = bv.variation_id
        AND bv.batch_id = ?
        GROUP BY fd.description};

  my $sth = $dbh->prepare($sql);
  $sth->execute($batch_id);
  my $tbl_ary_ref = $sth->fetchall_arrayref();
  print $report "\n\nVariation failure reason:\n";
  print $report join("\t", "fail_desc", "num"), "\n";
  foreach my $val (@{$tbl_ary_ref}) {
    print $report join("\t", @{$val}),"\n";
  }

  # Get fails by seq region
  $sql = qq{
        SELECT sr.name, fd.description, count(*)
        FROM failed_description fd, failed_variation fv, batch_variation bv,
             variation_feature vf, seq_region sr
        WHERE fd.failed_description_id = fv.failed_description_id
        AND fv.variation_id = bv.variation_id
        AND bv.batch_id = ?
        AND fv.variation_id = vf.variation_id
        AND vf.seq_region_id = sr.seq_region_id
        GROUP BY sr.name, fd.description};

  $sth = $dbh->prepare($sql);
  $sth->execute($batch_id);
  $tbl_ary_ref = $sth->fetchall_arrayref();
  print $report "\n\nSeq region - Variation failure reason:\n";
  print $report join("\t", 'sr_name', 'fail_desc', 'num'), "\n";
  foreach my $val (@{$tbl_ary_ref}) {
    print $report join("\t", @{$val}),"\n";
  }

  # Breakdown by seq_region
  $sql = qq{
        SELECT sr.name, count(*)
        FROM batch_variation bv,
             variation_feature vf, seq_region sr
        WHERE bv.variation_id = vf.variation_id
        AND vf.seq_region_id = sr.seq_region_id
        AND bv.batch_id = ?
        GROUP BY sr.name};

  $sth = $dbh->prepare($sql);
  $sth->execute($batch_id);
  $tbl_ary_ref = $sth->fetchall_arrayref();
  print $report "\n\nSeq region:\n";
  print $report join("\t", 'sr_name','num'), "\n";
  foreach my $val (@{$tbl_ary_ref}) {
    print $report join("\t", @{$val}),"\n";
  }
  close($report);
}

sub count_rows_batch {
  my ($dbh, $table, $column_name, $batch_id) = @_;

  my $sql = qq{SELECT COUNT(DISTINCT t.${column_name})
               FROM $table t, batch_variation bv
               WHERE t.variation_id = bv.variation_id
               AND batch_id = ?};
  return $dbh->selectrow_array($sql, undef, $batch_id);
}

sub get_seq_region_names {

  my ($dbh) = @_;
  my %seq_region;

  my $sth = $dbh->prepare(qq{SELECT seq_region_id, name FROM seq_region});

  $sth->execute() ||die "Error extracting seq_region_ids";
  my $data = $sth->fetchall_arrayref();

  foreach my $l (@{$data}) {
    $seq_region{$l->[0]} = $l->[1];
  }
  return (\%seq_region);
}

# For a given chr look up the seq_region
# This is currently only used to look up seq_region for Y for
# processing of PAR.
# It assumes that there is only one.
# This should be checked at the start of script.
sub get_seq_region_chr {
  my ($dbh, $chr) = @_;
  my $sth = $dbh->prepare(qq{SELECT seq_region_id FROM seq_region WHERE name = ?});
  $sth->execute($chr) || die "Error getting seq_region_id for chr $chr";
  my ($seq_region_id) = $sth->fetchrow_array();
  $sth->finish();
  return $seq_region_id;
}

# Get command line options
# Get information for run
sub configure {
  # get command-line options
  my $config = {};
  my $args = scalar @ARGV;

  my @variants = ();
  GetOptions(
    $config,

    'help|h',
    'debug',

    'input_file|i=s',
    'rpt_dir=s',
    'fasta_file=s',
    'ancestral_fasta_file=s',

    'species=s',
    'registry|r=s',
    'assembly|a=s',
    'no_db_load',
    'no_ref_check',
    'no_assign_ancestral_allele',
    'add_maf',
    ) or pod2usage(2);
 
  # Print usage message if help requested or no args
  pod2usage(1) if ($config->{'help'} || !$args);

  # Check the required parameters have been provided
  my @required_params = ('input_file', 'rpt_dir', 'registry');
  for my $param (@required_params) {
    if (! defined $config->{$param}) {
      pod2usage({ -message => "Mandatory argument ($param) is missing", 
                  -exitval => 2,
                }
               )
    }
  }

  # Set defaults
  $config->{'species'} ||= 'homo_sapiens';
  $config->{'debug'} ||= 0;
  $config->{'assembly'} ||= 'GRCh38';
  $config->{'db_load'} = 1;
  if (exists $config->{'no_db_load'}) {
    $config->{'db_load'} = 0;
  }

  $config->{'ref_check'} = 1;
  if (exists $config->{'no_ref_check'}) {
    $config->{'ref_check'} = 0;
  } 

  $config->{'assign_ancestral_allele'} = 1;
  if (exists $config->{'no_assign_ancestral_allele'}) {
    $config->{'assign_ancestral_allele'} = 0;
  }

  # Check parameters
  if (! -e $config->{'input_file'}) {
    die("ERROR: Input file does not exist ($config->{'input_file'})\n");
  }

  if (! -e $config->{'registry'}) {
    die("ERROR: Registry file does not exist ($config->{'registry'})\n");
  }

  if (! -d $config->{'rpt_dir'}) {
    die("ERROR: Report directory does not exist ($config->{'rpt_dir'})\n");
  }

  if ($config->{'assembly'} !~ /^(GRCh38|GRCh37)$/) {
      die("ERROR: Assembly is invalid ($config->{'assembly'}). Please specify GRCh38 or GRCh37");
  }

  if ($config->{'ref_check'}) {
    if (! defined $config->{'fasta_file'}) {
      pod2usage({ -message => "Mandatory argument (fasta_file) is missing", 
                  -exitval => 2,
                }
               );
    } elsif (! -e $config->{'fasta_file'}) {
      die("ERROR: FASTA file does not exist ($config->{'fasta_file'})\n");
    }
  }

  if ($config->{'assign_ancestral_allele'}) {
    if (! defined $config->{'ancestral_fasta_file'}) {
      pod2usage({ -message => "Mandatory argument (ancestral_fasta_file) is missing", 
                  -exitval => 2,
                }
               );
    } elsif (! -e $config->{'ancestral_fasta_file'}) {
      die("ERROR: Ancestral FASTA file does not exist ($config->{'ancestral_fasta_file'})\n");
    }
  }

  if (exists $config->{'add_maf'}) {
    $config->{'add_maf'} = 1;
  } else {
    $config->{'add_maf'} = 0;
  }

  return $config;  
}

sub debug {
  my $msg = shift;
  print STDERR "$msg\n";
}

sub get_db_info {
  my ($dbh) = @_;
  
  debug("\n>>>> get_db_info <<<<") if ($debug);
  
  my $sql = "SELECT now(), database(), user()";
  my $sth = $dbh->prepare($sql);
  $sth->execute();
  $sth->dump_results();
}

# Look up has of nc_region ids to enter variation features
sub get_nc_regions {
  my ($dbh) = @_;

  debug("\n>>>> get_nc_regions <<<<") if ($debug);
  
  my %nc_regions;

  my $sth = $dbh->prepare(qq{SELECT seq_region_id, name, srs_synonym FROM tmp_nc_synonym});

  $sth->execute() or die("Error extracting nc_regions");
  
  my $data = $sth->fetchall_arrayref();
  unless(defined $data->[0]->[0]) {
    die ("No nc_regions available");
  } 
    
  foreach my $reg (@{$data}){
    $nc_regions{$reg->[2]} = $reg->[0];
  }
  return (\%nc_regions);
}

# Look up hash of nw_region ids to enter variation features
sub get_nw_regions {
  my ($dbh) = @_;

  debug("\n>>>> get_nw_regions <<<<") if ($debug);
  
  my %nw_regions;

  my $sth = $dbh->prepare(qq{SELECT seq_region_id, srs_synonym, asm_start, asm_end, ori FROM tmp_nw_synonym});

  $sth->execute() or die("Error extracting nw_regions");
  
  my $data = $sth->fetchall_arrayref();
  unless(defined $data->[0]->[0]) {
    die ("No nw_regions available");
  } 
  my $rec;
    
  my ($seq_region_id, $srs_synonym, $asm_start, $asm_end, $ori); 
  foreach my $reg (@{$data}) {
    $rec = {};
    ($seq_region_id, $srs_synonym, $asm_start, $asm_end, $ori) = @$reg; 
    $nw_regions{$srs_synonym} = $rec;
    $rec->{'seq_region_id'} = $seq_region_id;
    $rec->{'asm_start'} = $asm_start;
    $rec->{'asm_end'} = $asm_end;
    $rec->{'ori'} = $ori;
  }
  return (\%nw_regions);
}

# Look up has of nt_region to enter variation features
sub get_nt_regions {
  my ($dbh) = @_;

  debug("\n>>>> get_nt_regions <<<<") if ($debug);
  
  my %nt_regions;

  my $sth = $dbh->prepare(qq{SELECT seq_region_id, srs_synonym, asm_start, asm_end, ori FROM tmp_nt_synonym});

  $sth->execute() or die("Error extracting nt_regions");
  
  my $data = $sth->fetchall_arrayref();
  unless(defined $data->[0]->[0]) {
    die ("No nt_regions available");
  } 
  my $rec;
    
  my ($seq_region_id, $srs_synonym, $asm_start, $asm_end, $ori); 
  foreach my $reg (@{$data}) {
    $rec = {};
    ($seq_region_id, $srs_synonym, $asm_start, $asm_end, $ori) = @$reg; 
    $nt_regions{$srs_synonym} = $rec;
    $rec->{'seq_region_id'} = $seq_region_id;
    $rec->{'asm_start'} = $asm_start;
    $rec->{'asm_end'} = $asm_end;
    $rec->{'ori'} = $ori;
  }
  return (\%nt_regions);
}

sub get_ref_regions {
  my ($dbh) = @_;

  debug("\n>>>> get_ref_regions <<<<") if ($debug);
  
  my %ref_regions;

  my $sth = $dbh->prepare(qq {SELECT sr.seq_region_id 
                              FROM seq_region sr, coord_system cs
                              WHERE sr.coord_system_id = cs.coord_system_id  
                              AND cs.rank = 1
                              AND seq_region_id NOT IN (
                                 SELECT seq_region_id FROM seq_region_attrib WHERE attrib_type_id = 16 )
                             }
                          );
  $sth->execute() or die("Error extracting ref_regions");
 
  my $is_ref = $sth->fetchall_arrayref();
 
  foreach my $srid (@{$is_ref}) {
    $ref_regions{$srid->[0]} = 1;
  }
  return \%ref_regions;
}


sub get_lu_info {
  my ($dbh) = @_;

  debug("\n>>>> get_lu_info <<<<") if ($debug);
  
  my %lu_info;

  $lu_info{'evidence_ids'} = get_evidence_attribs($dbh);
  $lu_info{'failed_set_id'} = find_failed_variation_set_id($dbh);
  $lu_info{'chrY_seq_region_id'} = get_seq_region_chr($dbh, 'Y');
  return \%lu_info;
}

# Taken from VariantQC.pm
# but used dbh instead of adaptor
sub find_failed_variation_set_id {

  my ($dbh) = @_;
  my $sth = $dbh->prepare(qq[ select variation_set_id
                              from variation_set
                              where name = ?
                              ]);
  ## check if present
  $sth->execute('All failed variations')  || die "Failed to extract failed variant set id";
  my $failed_set_id = $sth->fetchall_arrayref();

  die "Failed set not available" unless defined $failed_set_id->[0]->[0];

  return $failed_set_id->[0]->[0];
}

sub get_study_frequency {
  my ($data, $study) = @_;

  return if (! $study);

  my $study_freq;
  my $found = 0;
  my @alleles;

  my $allele_annotations = $data->{'primary_snapshot_data'}->{'allele_annotations'};
  #print "Number of allele_annotations = ", scalar(@$allele_annotations), "\n";

  # Assume allele_annotations exists
  # Is there frequency annotation ?
  return undef if (! exists $allele_annotations->[0]->{'frequency'});

  # Is there annotation for the study of interest.
  my $freqs = $allele_annotations->[0]->{'frequency'};
  for my $freq (@$freqs) {
      #print Dumper($freq), "\n";
      if ($freq->{'study_name'} eq $study) {
        $found = 1;
        push @alleles, $freq;
      }
  }
  return if (! $found);
  for (my $i=1; $i<@$allele_annotations; $i++) {
      my $freqs = $allele_annotations->[$i]->{'frequency'};
      for my $freq (@$freqs) {
        if ($freq->{'study_name'} eq $study) {
          push @alleles, $freq;
        }
      }
  }
  # Check the alleles are consistent
  # Have the same seq_id, position, deleted_sequence. The alleles then become
  # the inserted_sequence
  my $cmp_obs = $alleles[0]->{'observation'};
  my $cmp_seq_id = $cmp_obs->{'seq_id'};
  my $cmp_position = $cmp_obs->{'position'};
  my $cmp_deleted = $cmp_obs->{'deleted_sequence'};
  my %allele_errors;
  for (my $i=0; $i<@alleles; $i++) {
    my $obs = $alleles[$i]->{'observation'};
    if ($obs->{'seq_id'} ne $cmp_seq_id) {
        $allele_errors{1}++;
    }
    if ($obs->{'position'} ne $cmp_position) {
        $allele_errors{2}++;
    }
    if ($obs->{'deleted_sequence'} ne $cmp_deleted) {
      $allele_errors{3}++;
    }
  }
  if (%allele_errors) {
    #print "There are allele errors not going to do MAF, minor_allele\n";
  } else {
    #print "Going to calculate MAF and minor_allele\n";
  }
  $study_freq->{'alleles'} = [@alleles];
  $study_freq->{'minor_allele'} = undef;
  $study_freq->{'MAF'} = undef;
  $study_freq->{'minor_allele_count'} = undef;

  if (@{$study_freq->{'alleles'}} <= 1) {
    $allele_errors{4}++;
  } else {
    my ($maf, $minor_allele, $minor_allele_count) = get_maf($study_freq->{'alleles'}, $data->{'refsnp_id'});
    $study_freq->{'MAF'} = $maf;
    $study_freq->{'minor_allele'} = $minor_allele;
    $study_freq->{'minor_allele_count'} = $minor_allele_count;
  }
  $study_freq->{'allele_errors'} = {%allele_errors};
  return $study_freq;

}

# Get alignment diff between GRCh37 and GRCh38
sub get_align_diff {
  my ($rs_json) = @_;

  my $align_info;

  my $refsnp = $rs_json->{'refsnp_id'};

  my $psd = $rs_json->{'primary_snapshot_data'};
  if (! $psd) {
    return;
  }

  # Number of placements with allele for that primary snapshot data
  my $num_pwa = @{$psd->{'placements_with_allele'}};
  if (! $num_pwa) {
    return;
  }

  my ($grch37_loc, $grch38_loc, $grch37_aln_opp, $grch38_aln_opp) = (0,0,0,0);

  # Loop through each of the placements
  foreach my $loc (@{$psd->{'placements_with_allele'}}) {
    # for each placement get
    # the seq_id, is_ptlp, number of seq_id_trait_by_assemby,is_aln_opposite_orientation
    my $seq_id = $loc->{'seq_id'};

    # Only process the NC
    next if ($seq_id !~ /^NC/);
    my $assembly = $loc->{'placement_annot'}->{'seq_id_traits_by_assembly'}->[0]->{'assembly_name'};
    if ($assembly =~ /GRCh37/) {
      $grch37_loc++;
      $grch37_aln_opp++ if ($loc->{'placement_annot'}->{'is_aln_opposite_orientation'});
    }
    if ($assembly =~ /GRCh38/) {
      $grch38_loc++;
      $grch38_aln_opp++ if ($loc->{'placement_annot'}->{'is_aln_opposite_orientation'});
    }
    my $is_ptlp = $loc->{'is_ptlp'};
  }

  my $review_status;
  if (!$grch37_loc && !$grch38_loc) {
    $review_status = 'NO_MAPPING_GRCh37_GRCh38';
  } elsif (!$grch37_loc && $grch38_loc) {
    $review_status = 'ONLY_GRCh38';
  } elsif ($grch37_loc && !$grch38_loc) {
    $review_status = 'ONLY_GRCh37';
  } elsif ($grch37_aln_opp != $grch38_aln_opp) {
    $review_status = 'ALIGN_DIFF';
  }
  if ($review_status) {
    $align_info->{'grch37_loc'} = $grch37_loc;
    $align_info->{'grch38_loc'} = $grch38_loc;
    $align_info->{'grch37_aln_opp'} = $grch37_aln_opp;
    $align_info->{'grch38_aln_opp'} = $grch38_aln_opp;
    $align_info->{'review_status'} = $review_status;
  }
  return $align_info;
}

sub get_maf {
  my ($alleles_ref, $rsid) = @_;

  #print Dumper($alleles_ref);
  my @sorted_alleles = sort {
                            $b->{'allele_count'} <=> $a->{'allele_count'}
                                     or
                            $a->{'observation'}->{'inserted_sequence'} cmp $b->{'observation'}->{'inserted_sequence'}

                          } @$alleles_ref;
  my $maf = sprintf("%.6f", $sorted_alleles[1]->{'allele_count'}/$sorted_alleles[1]->{'total_count'});
  my $minor_allele = $sorted_alleles[1]->{'observation'}->{'inserted_sequence'};
  my $minor_allele_count = $sorted_alleles[1]->{'allele_count'};

  my $var_name = "rs$rsid";

  if (length($minor_allele) > 50) {
    my $msg = 'truncated_minor_allele';
    my $info = join(";", 'rsid='. $var_name,
                          'minor_allele=' . $minor_allele);
    $minor_allele = substr($minor_allele,0,50);
    log_errors($config, $var_name, $msg, $info);
  }


  #print "Minor Allele Frequence = ($maf)\n";
  #print "Minor allelele = ($minor_allele)\n";
  #print "Minor allele count = ($minor_allele_count)\n";
  #print Dumper(\@sorted_alleles);
  return ($maf, $minor_allele, $minor_allele_count);
}

# Sets the class_attrib_id for variation and variation_features
sub set_SO_variation_class {
  my ($config, $data) = @_;

  my $seq_alt = SO_TERM_SEQUENCE_ALTERATION;
  my $seq_alt_class_id = 
       $config->{'attribute_adaptor'}->attrib_id_for_type_value('SO_term', $seq_alt);
  
  # No variation features set class_attrib_id to sequence alteration
  if (! (@{$data->{'vfs'}})) {
    $data->{'v'}->{'class_attrib_id'} = $seq_alt_class_id;
    return;
  }
  my %class;
  # For each variation_feature set the class_attrib_id
  for my $vf (@{$data->{'vfs'}}) {
    # During mapping checking vf->{'ref_correct'} set to
    # 'no' when there is a mapping failure
    # TODO: Review if vf->{'ref_correct'} should be set to
    # 'yes', 'no', 'unknown'
    my $ref_correct = 1;
    
    # for dbSNP_novariation and dbSNP_variant
    if ($vf->{'allele_string'} =~ /^dbSNP_$/) {
        $vf->{'class_attrib_id'} = $seq_alt_class_id;
      $class{'class_attrib_id'}++;
      next;
    }
    if ((exists $vf->{'ref_correct'}) && $vf->{'ref_correct' eq 'no'})  {
      $ref_correct = 0;
    }
    my $so_term  = SO_variation_class($vf->{'allele_string'}, $ref_correct);
    my $class_id = $config->{'attribute_adaptor'}->attrib_id_for_type_value('SO_term', $so_term);
    $vf->{'class_attrib_id'} = $class_id;
    $class{$class_id}++;
  }
  # Set the variant class the the first variation_feature set 
  $data->{'v'}->{'class_attrib_id'} = $data->{'vfs'}->[0]->{'class_attrib_id'};
  
  if (scalar(keys %class) > 1) {
    my $msg = "multiple class_attrib_id";
    my $info;
    $info = join(";",
           'variant_class_attrib_id=' . $data->{'v'}->{'class_attrib_id'},
           'num_class_attrib_id=' . scalar(keys %class),
           map {join(":" , "class_attrib_id=$_", "num=$class{$_}")} sort {$a <=> $b} keys %class);
    log_errors($config, $data->{'v'}->{'name'},
               $msg,
               $info);
  }
}

sub log_errors {
  my ($config, $rsid, $msg, $info) = @_;
  $info ||= '';
  $rsid ||= '';
  
  my $error_file = $config->{error_file};
  open(my $fh, '>>', $config->{'error_file'}) or die "Failed to open error_file: $config->{'error_file'}: $!";
  print $fh join("\t",
    strftime('%Y-%m-%d %H:%M:%S', localtime()),
    $rsid, $msg, $info), "\n";
  $fh->close;
}

# Based on Bio::EnsEMBL::Variation::VCFVariation.pm
sub format_frequency {
  my ($freq_val) = @_;
  return sprintf("%.4g", $freq_val || 0);
}

sub get_batch_id {
  my ($dbh, $filename, $parent_file) = @_;

  $dbh->do(qq{INSERT INTO batch
              (filename, parent_filename)
              VALUES
               (?, ?)},
                undef,
                $filename,
                $parent_file);
  
  my $batch_id = $dbh->last_insert_id(undef, undef, qw(batch batch_id)) or die("no batch_id for batch");
  if (! $batch_id) {
    die("Unable to get batch_id $filename: $!");
  }
  return $batch_id;
}

sub get_sources {
  my ($dbh) = @_;

  my @source_names = ('dbSNP', 'Archive dbSNP', 'dbSNP HGVS', 'Former dbSNP');

  my $sources = {};
  my $source_id;
  my $sth = $dbh->prepare(qq{select source_id from source where name = ?});
  for my $source (@source_names) {
    $sth->execute($source);
    ($source_id) =  $sth->fetchrow_array();
    $sth->finish();
    if (defined $source_id) {
      $sources->{$source} = $source_id;
    } else {
      die("Unable to find source_id for $source");
    }
  }
  return $sources;
}

__END__

=head1 NAME

load_dbsnp.pl

=head1 DESCRIPTION

Imports variations from a dbSNP JSON file into an Ensembl variation DB

=head1 SYNOPSIS

load_dbsnp.pl [arguments]

=head1 EXAMPLE COMMAND LINE

 load_dbsnpl.pl --input_file test.json \
                --rpt_dir  /user1/load_reports \
                --fasta_file /user1/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
                --registry  /user1/ensembl.registry

=head1 OPTIONS

=over 4

=item B<--help>

Displays this documentation

=item B<--input_file FILE>

JSON file provided by dbSNP

=item B<--rpt_dir DIR>

Directory to store summary report and error logs

=item B<--assembly ASSEMBLY>

Specify assembly to use. GRCh38 (default) or GRCh37

=item B<--ancestral_fasta_file FILE>

Path to FASTA file containing ancestral sequence

=item B<--fasta_file FILE>

Path to FASTA file containing reference sequence

=item B<--registry FILE>

Read database connections from this registry file

=item B<--species NAME>

Species to use [default: 'homo_sapiens']

=item B<--no_db_load>

No database load. Only parses the file. Used for testing.

=item B<--no_ref_check>

No reference checking

=item B<--no_assign_ancestral_allele>

No ancestral allele assignment

=item B<--no_maf>

Not include MAF information

=item B<--debug>

Debug mode

=back

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <https://www.ensembl.org/Help/Contact>.

=cut

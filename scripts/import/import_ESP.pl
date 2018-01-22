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
use strict;
use warnings;

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::Utils::VEP qw/parse_line get_all_consequences/;
use Bio::EnsEMBL::Variation::Utils::Sequence qw/SO_variation_class/;
use Bio::EnsEMBL::Utils::Sequence qw/reverse_comp/;
use FileHandle;
use Getopt::Long;
use Data::Dumper;
=doc
  - Parse variation and variation_feature from VCF file
  - Find co-located variants from dbSNP, if variation class is the same: link ESP variant to co-located variant
  - If alleles from co-located variant cannot be merged create new variant and variation_feature objects
=end
=cut

my $config = {};
GetOptions(
  $config,
  'tmp_dir=s',
  'vcf_files_dir=s',
  'registry=s',
  'host=s',
  'dbname=s',
  'user=s',
  'pass=s',
  'port=i',
  'test',
  'import_vcf',
  'assign_evidence',
  'create_set',
  'release_version=i',
  'esp_version=i',
  'assembly=i',
  'help!',
) or die "Error: Failed to parse command line arguments\n";

usage() if ($config->{help});

foreach my $param (qw/vcf_files_dir tmp_dir assembly registry release_version esp_version/) {
  die ("Parameter (--$param) is required") unless (defined($config->{$param}));
}

if ($config->{test}) {
  print STDERR "Run in test mode\n";
}

setup_db_connections($config);
set_ESP_related_parameters($config);
clear_old_data($config);

my $tmp_dir = $config->{tmp_dir};
my $assembly = $config->{assembly};
my $release_version = $config->{release_version};
my $fh_variation_ids = FileHandle->new("$tmp_dir/variation_ids_$assembly\_$release_version.txt", 'w');
main($config);
$fh_variation_ids->close();

if ($config->{assign_evidence}) {
  prepare_update_evidence($config);
  cleanup_old_evidence($config);
  update_new_evidence($config);
}

if ($config->{create_set}) {
  create_set($config)
}

sub setup_db_connections {
  my $config = shift;
  my $registry = 'Bio::EnsEMBL::Registry';
  my $species = 'human';
  my ($dbh, $dbc, $vdba);
  if ($config->{registry}) {
    $registry->load_all($config->{registry});
#    $registry->set_disconnect_when_inactive();
#    $registry->set_reconnect_when_lost();
    $vdba = $registry->get_DBAdaptor($species, 'variation');
  } else {
    $vdba = new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
        -host => $config->{host},
        -user => $config->{user},
        -pass => $config->{pass},
        -port => $config->{port},
        -dbname => $config->{dbname},
        ) or die("Could not get VDBA");
  }

  $dbh = $vdba->dbc->db_handle;
  $dbc = $vdba->dbc;
  $config->{dbh} = $dbh;
  $config->{dbc} = $dbc;
  $config->{allele_adaptor} = $vdba->get_AlleleAdaptor;
  $config->{attribute_adaptor} = $vdba->get_AttributeAdaptor;
  $config->{population_adaptor} = $vdba->get_PopulationAdaptor;
  $config->{population_gt_adaptor} = $vdba->get_PopulationGenotypeAdaptor;
  $config->{source_adaptor} = $vdba->get_SourceAdaptor;
  $config->{variation_adaptor} = $vdba->get_VariationAdaptor;
  $config->{variation_adaptor}->db->include_failed_variations(1);
  $config->{variation_set_adaptor} = $vdba->get_VariationSetAdaptor;
  $config->{vfa} = $vdba->get_VariationFeatureAdaptor;
  $config->{vfa}->db->include_failed_variations(1);
  $config->{seq_region_ids} = get_seq_region_ids($config);   
}

sub set_ESP_related_parameters {
  my $config = shift;
  my $pa = $config->{population_adaptor};
  my $aa = $pa->fetch_by_name('ESP6500:African_American');

  if (!$aa) {
    $aa = Bio::EnsEMBL::Variation::Population->new(
        -name => 'ESP6500:African_American',
        -display => 'UNDISPLAYABLE',
        -display_group_id => 3,
    );
    $pa->store($aa);
    $config->{AA}->{population_id} = $aa->{dbID};
  } else {
    $config->{AA}->{population_id} = $aa->dbID();
  }

  my $ea = $pa->fetch_by_name('ESP6500:European_American');
  if (!$ea) {
    $ea = Bio::EnsEMBL::Variation::Population->new(
      -name => 'ESP6500:European_American',
      -display => 'UNDISPLAYABLE',
      -display_group_id => 3,
    );
    $pa->store($ea);
    $config->{EA}->{population_id} = $ea->{dbID};
  } else {
    $config->{EA}->{population_id} = $ea->dbID();
  }

  $config->{AA}->{population} = $aa;
  $config->{EA}->{population} = $ea;

  my $source_name = 'ESP';
  my $source_id;
  my $source_adaptor = $config->{source_adaptor};
  my $version = $config->{esp_version}; 
  my $source = $source_adaptor->fetch_by_name($source_name);
  if (!$source) {
    $source = Bio::EnsEMBL::Variation::Source->new(
      -name => 'ESP',
      -description => 'Data from NHLBI ESP version v.0.0.30. The goal of the NHLBI GO Exome Sequencing Project is to discover novel genes and mechanisms contributing to heart, lung and blood disorders by sequencing the protein coding regions of the human genome.',
      -version => $version,
      -url => 'http://evs.gs.washington.edu/EVS/',
      -data_types => 'variation',         
    );
    $source_adaptor->store($source);
    $config->{source_id} = $source_id->{dbID};
  } else {
    # version update
    my $source_id = $source->dbID;
    my $dbh = $config->{dbh};
    $dbh->do(qq/UPDATE source SET version=$version WHERE source_id=$source_id;/) or die $dbh->errstr;
    $config->{source_id} = $source_id;
  } 
  $config->{source} = $source_name;
}

sub clear_old_data {
  my $config = shift;
  my $dbh = $config->{dbh};
  my $source_id = $config->{source_id};
  $dbh->do(qq/DELETE FROM variation WHERE source_id=$source_id;/) or die  $dbh->errstr;
  $dbh->do(qq/DELETE FROM variation_feature WHERE source_id=$source_id;/) or die  $dbh->errstr;

  my $EA_population_id = $config->{EA}->{population_id};
  my $AA_population_id = $config->{AA}->{population_id};
  $dbh->do(qq/DELETE FROM allele WHERE population_id=$EA_population_id OR population_id=$AA_population_id;/) or die  $dbh->errstr;
  $dbh->do(qq/DELETE FROM population_genotype WHERE population_id=$EA_population_id OR population_id=$AA_population_id;/) or die  $dbh->errstr;
}


sub main {
  my $config = shift;
  my $vcf_files_dir = $config->{vcf_files_dir};
  my @files = ();
  opendir (DIR, $vcf_files_dir) or die $!;
  while (my $vcf_file = readdir(DIR)) {
    if ($vcf_file =~ m/\.vcf$/) {
      push @files, $vcf_file;
    }
  }   
  closedir(DIR);
  foreach my $vcf_file (@files) {
    import_vcf($config, "$vcf_files_dir/$vcf_file");
  }
}

sub import_vcf {
  my $config = shift;
  my $vcf_file = shift;

  my %headers = ();

  my $fh = FileHandle->new($vcf_file, 'r');

  while (<$fh>) {
    chomp;
    next if /^##/;
    my @split = split /\t/;
    my $data = {};
    $data->{line} = $_;

    if (/^#/) {
      %headers = %{parse_header($config, \@split)};
    } else {
      $data->{$_} = $split[$headers{$_}] for keys %headers;
      if ($data->{ALT} eq '.') {
        $config->{skipped}->{non_variant}++;
        next;
      }
      # parse info column
      my %info;
      foreach my $chunk (split /\;/, $data->{INFO}) {
        my ($key, $value) = split /\=/, $chunk;
        $info{$key} = $value;
      }
      $data->{info} = \%info;

      ($data->{tmp_vf}) = @{parse_line($config, $data->{line})};
      $data->{tmp_vf}->{seq_region_id} = $config->{seq_region_ids}->{$data->{tmp_vf}->{chr}};

      my $tmp_vf_start  = $data->{tmp_vf}->{start};
      my $tmp_vf_end    = $data->{tmp_vf}->{end};
      my $allele_string = $data->{tmp_vf}->{allele_string};

      my $esp_tmp_id = 'TMP_ESP_' . $data->{'#CHROM'} . '_' . $tmp_vf_start . '_' . $tmp_vf_end;

      my $id = $data->{ID};

      $data->{tmp_id} = $esp_tmp_id;
      unless ($id =~ /^rs/) {
        $data->{ID} = $esp_tmp_id;
      }

      if ($config->{assembly} == 38) {      
        my $grch38_position = $data->{info}->{GRCh38_POSITION};
        if ($grch38_position eq '-1') {
          next;
        }
        my ($grch38_chrom, $grch38_pos) = split(':', $grch38_position);  
        unless ($grch38_chrom =~ /^\d+$|^X$|^Y$|^MT$/) {
          next;
        }      
        my @line = split("\t", $data->{line});
        shift @line;
        shift @line;
        my $line_38 = join("\t", $grch38_chrom, $grch38_pos, @line);
        ($data->{tmp_vf}) = @{parse_line($config, $line_38)};
        $data->{tmp_vf}->{seq_region_id} = $config->{seq_region_ids}->{$data->{tmp_vf}->{chr}};

      } # end changes for new choordinates

      $data->{variation} = variation($config, $data);
      $data->{vf}        = variation_feature($config, $data);
      
      if ($data->{vf}) {
        if (!$config->{test}) {
          population_genotype($config, $data);
        }
      }
    }
  }
  $fh->close();
}

sub variation {
  my $config = shift;
  my $data = shift;

  my $variation_name = $data->{ID};

  my $variation = $config->{variation_adaptor}->fetch_by_name($variation_name);

  my $vf = $data->{tmp_vf}; 
  my $so_term = SO_variation_class($vf->allele_string, 1);
  my $class_id = $config->{attribute_adaptor}->attrib_id_for_type_value('SO_term', $so_term);

  if (!defined($variation)) {
    $variation = Bio::EnsEMBL::Variation::Variation->new_fast({
        name => $variation_name,
        _source_id => $config->{source_id},
        is_somatic => 0,
    });
    $variation->{class_attrib_id} = $class_id;
  }

  $vf->{variation} = $variation;
  return $variation;
}

sub variation_feature {
  my $config = shift;
  my $data = shift; 

  my $dbh = $config->{dbh};
  my $vfa = $config->{vfa};

  my $vf = $data->{tmp_vf};
  my $new_variation = $vf->{variation};

  my @new_alleles = split /\//, $vf->allele_string;  # ESP alleles

  my $var_in_db = defined($data->{variation}->{dbID}) ? 1 : 0;

  my $existing_vfs = $vfa->_fetch_all_by_coords($config->{seq_region_ids}->{$vf->{chr}}, $vf->{start}, $vf->{end}, $config->{somatic});

  my @dbSNP_existing_vfs = grep {$_->source_name eq 'dbSNP'} @$existing_vfs;

  my $count_match = 0;
  my @match = (); 
  my $new_name = $new_variation->{name};
  my $new_allele_string = $vf->allele_string; # ESP alleles
  my $so_term = SO_variation_class($vf->allele_string, 1);

  $data->{ESP_allele_string} = $vf->allele_string; 

  foreach my $existing_vf (@dbSNP_existing_vfs) {
    my $existing_so_term       = SO_variation_class($existing_vf->allele_string, 1);
    my $existing_name          = $existing_vf->variation_name;
    my $existing_source        = $existing_vf->source_name;
    my $existing_allele_string = $existing_vf->allele_string; 
    my $existing_strand        = $existing_vf->seq_region_strand;

    if ($so_term eq $existing_so_term) {
      my $update_allele_string = 0;
      my $merged_allele_string = '';        

      if ($new_allele_string ne $existing_allele_string) {
        $update_allele_string = 1; 
        if ($existing_strand == -1) {
          # reverse_comp ESP allele string
          my @rev_allele_strings = split('/', $new_allele_string);
          foreach my $allele (@rev_allele_strings) {
            reverse_comp(\$allele); 
          }
          $new_allele_string = join('/', @rev_allele_strings);
          $data->{ESP_rev_allele_string} = $new_allele_string;
        }

        my @new_alleles = split('/', $new_allele_string);
        my $new_ref = shift @new_alleles; 
        my @existing_alleles = split('/', $existing_allele_string);
        my $existing_ref = shift @existing_alleles;

        if ($new_ref eq $existing_ref) {
          my $merged_hash = {};
          foreach my $allele (@new_alleles, @existing_alleles) {
            $merged_hash->{$allele} = 1;
          }
          my @sorted_alleles = sort {$a cmp $b} keys %$merged_hash;
          my @merged_alleles = ();
          push @merged_alleles, $new_ref;
          push @merged_alleles, @sorted_alleles;
          $merged_allele_string = join('/', @merged_alleles);
          if ($merged_allele_string ne $existing_allele_string) {
            print STDERR "Merge alleles for $existing_name from $existing_allele_string to $merged_allele_string\n";
            if (!$config->{test}) {
              my $sth = $dbh->prepare(qq{
                UPDATE variation_feature
                SET allele_string = ?
                WHERE variation_feature_id = ?;
              });
              $sth->execute($merged_allele_string, $existing_vf->dbID);
            }
            $existing_vf->{allele_string} = $merged_allele_string;
          }
          # link existing variant to ESP
          print $fh_variation_ids $existing_vf->get_Variation_dbID(), "\n"; 
          $count_match++;
          $vf = $existing_vf;
        }
      } # Allele strings match: 
      print $fh_variation_ids $existing_vf->get_Variation_dbID(), "\n"; 
      $count_match++;
      $vf = $existing_vf;
    }
  }

  if ($count_match == 0) {
    # create new variation feature
    if ($var_in_db) {
      # create new variation object
      my $variation_name = $new_variation->{name};
      if ($variation_name =~ /^rs/) {
        $variation_name = $data->{tmp_id};
        # check that tmp_id is not already stored in database
        my $is_new = 0;
        while (!$is_new) {
          my $new_var = $config->{variation_adaptor}->fetch_by_name($variation_name); 
          if ($new_var) {
            my @components = split('_', $variation_name); # TMP_ESP_Chrom_start_end
            if (scalar @components == 5) {
              $variation_name = "$variation_name\_2";
            } else {
              my $count = pop @components;
              $count++;
              $variation_name = join('_', @components) . '_' . $count;
            }
          } else {
            $is_new = 1;
          }
        }
      } else {
        my @components = split('_', $variation_name); # TMP_ESP_Chrom_start_end
        if (scalar @components == 5) {
          $variation_name = "$variation_name\_2";
        } else {
          my $count = pop @components;
          $count++;
          $variation_name = join('_', @components) . '_' . $count;
        }
      }
      my $so_term = SO_variation_class($vf->allele_string, 1);
      my $class_id = $config->{attribute_adaptor}->attrib_id_for_type_value('SO_term', $so_term);

      $new_variation = Bio::EnsEMBL::Variation::Variation->new_fast({
          name => $variation_name,
          _source_id => $config->{source_id},
          is_somatic => 0,
      });
      $new_variation->{class_attrib_id} = $class_id;
    } 
    if (!$config->{test}) {
      $config->{variation_adaptor}->store($new_variation);
    } else {
      print STDERR Dumper($new_variation), "\n";
    }
    $vf->{variation} = $new_variation;
    $vf->{_source_id} = $config->{source_id};
    $vf->{is_somatic} = 0;
    my $so_term = SO_variation_class($vf->{allele_string}, 1);
    $vf->{class_attrib_id} = $config->{attribute_adaptor}->attrib_id_for_type_value('SO_term', $so_term);
    $vf->{variation_name} = $vf->{variation}->name;
    
    if (!$config->{test}) {
      $vf->{variation_id} = $vf->{variation}->{dbID};
      $vfa->store($vf);
    } else {
      print STDERR Dumper($vf), "\n";
    }
    $data->{variation} = $new_variation;
  } else {
    $data->{variation} = $vf->variation;
  }

  $data->{vf} = $vf;
  return $vf;
}

sub population_genotype {
  my $config = shift;
  my $data   = shift;

  my $allele_adaptor        = $config->{allele_adaptor};
  my $population_gt_adaptor = $config->{population_gt_adaptor};

  my $vf = $data->{vf};
  my $variation_name = $data->{ID};
  if ($variation_name =~ /\;/) {
    print STDERR "WARN multiple IDS for $variation_name\n";
    next;
  }

  my @alt_alleles = split /\//, $data->{ESP_allele_string};
  my $is_rev = $data->{ESP_rev_allele_string};
  if ($is_rev) {
    print STDERR "IS REV $variation_name\n";
  }

  my $ref = shift @alt_alleles;

# observed genotypes
  my @gts = split(',', $data->{info}->{GTS});
  my @sorted_gts = ();
  foreach my $gt (@gts) {
    my @tmp_gt = ();
    if ($gt =~ /([ACGT])([ACGT])/) {
      foreach my $allele ($1, $2) {
        push @tmp_gt, $allele;
      }
    } elsif ($gt =~ /^[ACGT]$/) { # X, Y
      push @tmp_gt, $gt;
    } elsif ($gt eq 'R') { # X, Y
      push @tmp_gt, $ref;
    } elsif ($gt =~ /^[A]\d$/) { # X, Y
      my @allele_idx = split('', $gt);
      my $idx = $allele_idx[1];
      my $alt_allele = $alt_alleles[$idx - 1];
      push @tmp_gt, $alt_allele;
    } elsif ($gt =~ /^([A]\d|[R])([A]\d|[R])$/) {
      foreach my $allele ($1, $2) {
        if ($allele eq 'R') {
          push @tmp_gt, $ref;
        } else {
          my @allele_idx = split('', $allele);
          my $idx = $allele_idx[1];
          my $alt_allele = $alt_alleles[$idx - 1];
          push @tmp_gt, $alt_allele;
        }
      }
    } else {
      print STDERR "Unmatched GT for $variation_name\n";
    }
    my $sorted_gt = join('|', sort @tmp_gt);
    push @sorted_gts, $sorted_gt;
  }

  foreach my $population_id (qw/AA EA/) {
    my @allele_objs = ();
    my @population_gt_objs = ();

    my $hash = {};
    my @ac = split(',', $data->{info}->{$population_id . '_AC'});
    my $ref_count = pop @ac;
    $hash->{$ref} = $ref_count;
    my $idx = 0;
    foreach my $count (@ac) {
      $hash->{$alt_alleles[$idx]} = $count;
      $idx++;
    }
    my $total = 0;
    $total += $_ for values %$hash;
    foreach my $allele (keys %$hash){
      if ($hash->{$allele}) {
        my $mod_allele = $allele;
        if ($is_rev) {
          reverse_comp(\$mod_allele);
        }
        my $allele_obj = Bio::EnsEMBL::Variation::Allele->new_fast({
            allele     => $mod_allele,
            count      => $hash->{$allele},
            frequency  => $hash->{$allele} / $total,
            population => $config->{$population_id}->{population},
            variation  => $data->{variation},
            });
        push @allele_objs, $allele_obj;
      }
    }
    if (defined($config->{test})) {
      foreach my $allele (@allele_objs) {
        print STDERR 'A  ', $allele->population->name, ' ', $allele->variation->name, ' ', $allele->variation->dbID, ' ', $allele->allele, ' ', $allele->count, ' ', $allele->frequency, "\n";
      }
    } else {
      foreach my $allele (@allele_objs) {
        print STDERR 'A  ', $allele->population->name, ' ', $allele->variation->name, ' ', $allele->variation->dbID, ' ', $allele->allele, ' ', $allele->count, ' ', $allele->frequency, "\n";
      }
      $allele_adaptor->store_multiple(\@allele_objs);
    }
# observed genotypes
    my @gtc = split(',', $data->{info}->{$population_id . '_GTC'});
    my $total_gtc = 0;
    $total_gtc += $_ for @gtc;
    $idx = 0;
    foreach my $pop_gt (@sorted_gts) {
      if ($gtc[$idx]) {
        my @gt_alleles = split /\|/, $pop_gt;
        if ($is_rev) {
          foreach my $a (@gt_alleles) {
            reverse_comp(\$a);
          }
        }
        my @sorted_gts_2 = sort {$a cmp $b} @gt_alleles;
        my $pop_gt_obj = Bio::EnsEMBL::Variation::PopulationGenotype->new_fast({
            genotype   => \@sorted_gts_2,
            population => $config->{$population_id}->{population},
            variation  => $data->{variation},
            frequency  => $gtc[$idx] / $total_gtc,
            count      => $gtc[$idx],
            });
        push @population_gt_objs, $pop_gt_obj;
      }
      $idx++;
    }
    if (defined($config->{test})) {
      foreach my $gt (@population_gt_objs) {
        print STDERR 'GT ', $gt->population->name, ' ', $gt->variation->name, ' ', $gt->variation->dbID, ' ', $gt->genotype_string, ' ', $gt->count, ' ', $gt->frequency, "\n";
      }
    } else {
      foreach my $gt (@population_gt_objs) {
        print STDERR 'GT ', $gt->population->name, ' ', $gt->variation->name, ' ', $gt->variation->dbID, ' ', $gt->genotype_string, ' ', $gt->count, ' ', $gt->frequency, "\n";
      }

      $population_gt_adaptor->store_multiple(\@population_gt_objs);
    }
  }
}

sub prepare_update_evidence {
  my $config = shift;
  my $dbh = $config->{dbh};

  my $afr_population_id = $config->{AA}->{population_id};
  my $eur_population_id = $config->{EA}->{population_id};
  my $stmt;

  foreach my $stmt (
    q/DROP TABLE IF EXISTS ESP_variation_ids_old/,
    q/DROP TABLE IF EXISTS ESP_variation_ids_new/,
    q/CREATE TABLE `ESP_variation_ids_old` (`variation_id` int(10) unsigned NOT NULL, PRIMARY KEY (`variation_id`))/,
    q/CREATE TABLE `ESP_variation_ids_new` (`variation_id` int(10) unsigned NOT NULL, PRIMARY KEY (`variation_id`))/,
    q/INSERT INTO ESP_variation_ids_old SELECT variation_id FROM variation WHERE evidence_attribs LIKE '%372%';/) {
    print STDERR "Run $stmt\n";
    $dbh->do($stmt) or die $dbh->errstr;
  }

  $stmt = q/
    INSERT INTO ESP_variation_ids_new
    SELECT distinct variation_id
    FROM population_genotype
    WHERE population_id = ? OR population_id = ?;/;

  $dbh->do($stmt, undef, $afr_population_id, $eur_population_id) or die $dbh->errstr;

}

sub cleanup_old_evidence {
  my $config = shift;
  my $tmp_dir = $config->{tmp_dir};
  my $version = $config->{release_version};

  my $dbh = $config->{dbh};

  # get all variation ids from ESP_variation_ids_old that are not in ESP_variation_ids_new
  my $fh = FileHandle->new("$tmp_dir/fetch_old_evidence_values_$version.txt", 'w');
  my $sth = $dbh->prepare(qq{
    SELECT v.variation_id, v.evidence_attribs
    FROM ESP_variation_ids_old old LEFT JOIN
    ESP_variation_ids_new new ON
    old.variation_id = new.variation_id
    LEFT JOIN variation v ON
    old.variation_id = v.variation_id
    WHERE new.variation_id IS NULL;
  }, {mysql_use_result => 1});

  $sth->execute();
  my ($variation_id, $evidence);
  $sth->bind_columns(\($variation_id, $evidence));
  while ($sth->fetch) {
    $evidence ||= '';
    print $fh $variation_id, "\t", $evidence, "\n";
  }
  $sth->finish();
  $fh->close();

  $fh = FileHandle->new("$tmp_dir/fetch_old_evidence_values_$version.txt", 'r');

  while (<$fh>) {
    chomp;
    my ($id, $evdn) = split("\t", $_);
    my @values = ();
    if ($evdn =~ /372/) {
      my @old_evdns = split(',', $evdn);
      foreach my $value (@old_evdns) {
        unless ($value eq '372') {
          push @values, $value;
        }
      }
    }
    my $new_evdns = join(',', @values);
    if ($new_evdns eq '') {
      $new_evdns = undef;
    }
    $dbh->do(q{UPDATE variation_feature SET evidence_attribs = ? WHERE variation_id = ?;}, undef, $new_evdns, $id) or die $dbh->errstr;
    $dbh->do(q{UPDATE variation SET evidence_attribs = ? WHERE variation_id = ?;}, undef, $new_evdns, $id) or die $dbh->errstr;
    $new_evdns ||= '';
    print STDERR "update variation set evidence_attribs = $new_evdns where variation_id = $id\n";
  }
  $fh->close();
}

sub update_new_evidence {
  my $config = shift;
  my $version = $config->{release_version};
  my $tmp_dir = $config->{tmp_dir};

  my $dbh = $config->{dbh};

  my $fh = FileHandle->new("$tmp_dir/fetch_new_evidence_values_$version.txt", 'w');
  my $sth = $dbh->prepare(qq{
    SELECT v.variation_id, v.evidence_attribs
    FROM ESP_variation_ids_new new LEFT JOIN
    ESP_variation_ids_old old ON
    new.variation_id = old.variation_id
    LEFT JOIN variation v ON
    new.variation_id = v.variation_id
    WHERE old.variation_id IS NULL;
  }, {mysql_use_result => 1});

  $sth->execute();
  my ($variation_id, $evidence);
  $sth->bind_columns(\($variation_id, $evidence));
  while ($sth->fetch) {
    $evidence ||= '';
    print $fh $variation_id, "\t", $evidence, "\n";
  }
  $sth->finish();
  $fh->close();
  $fh = FileHandle->new("$tmp_dir/fetch_new_evidence_values_$version.txt", 'r');

  while (<$fh>) {
    chomp;
    my ($id, $evdn) = split("\t", $_);
    if ($evdn eq '') {
      $evdn = '372';
    } else {
      unless ($evdn =~ /372/) {
        my @values = ();
        push @values, $evdn;
        push @values, '372';
        $evdn = join(',', @values);
      }
    }
    $dbh->do(q{UPDATE variation_feature SET evidence_attribs = ? WHERE variation_id = ?;}, undef, $evdn, $id) or die $dbh->errstr;
    $dbh->do(q{UPDATE variation SET evidence_attribs = ? WHERE variation_id = ?;}, undef, $evdn, $id) or die $dbh->errstr;
    print STDERR "update variation set evidence_attribs = $evdn where variation_id = $id\n";
  }
  $fh->close();
}

sub create_set {
  my $config = shift;
  my $dbh = $config->{dbh};

  my $vsa = $config->{variation_set_adaptor};
  my $set = $vsa->fetch_by_name('ESP_6500');
  my $set_id = '';

  if (!$set) {
    $dbh->do(qq{INSERT INTO variation_set(name, description, short_name_attrib_id) values('ESP_6500', 'Variants from the NHLBI Exome Sequencing Project (investigating heart, lung and blood disorders)', 344)});
    my $sth = $dbh->prepare(qq{SELECT variation_set_id FROM variation_set WHERE name = 'ESP_6500'});
    $sth->execute();
    $sth->bind_columns(\$set_id);
    $sth->fetch;
    $sth->finish;
  } else {
    $set_id = $set->dbID;
  }
  die "Variation set id for ESP_6500 is not defined." unless defined($set_id);

  my $afr_population_id = $config->{AA}->{population_id};
  my $eur_population_id = $config->{EA}->{population_id};

  my $stmt = q/CREATE TABLE `ESP_variation_set` (`variation_id` int(10) unsigned NOT NULL, PRIMARY KEY (`variation_id`))/;
  $dbh->do($stmt) or die $dbh->errstr;

  $stmt = q/DELETE from variation_set_variation WHERE variation_set_id = ?;/;
  $dbh->do($stmt, undef, $set_id) or die $dbh->errstr;

  $stmt = q/
    INSERT INTO ESP_variation_set
    SELECT distinct variation_id
    FROM population_genotype
    WHERE population_id = ? OR population_id = ?;/;
  $dbh->do($stmt, undef, $afr_population_id, $eur_population_id) or die $dbh->errstr;

  $stmt = qq/INSERT INTO variation_set_variation(variation_id, variation_set_id) SELECT variation_id, $set_id FROM ESP_variation_set;/;
  $dbh->do($stmt) or die $dbh->errstr;

}

sub parse_header {
  my $config = shift;
  my $split_ref = shift;
  my @split = @$split_ref;
  my %headers;
  $headers{$split[$_]} = $_ for (0..$#split);
  return \%headers;
}

sub get_seq_region_ids {
  my $config = shift;
  my $dbh = $config->{dbh};
  my ($seq_region_id, $chr_name, %seq_region_ids);

  my $sth = $dbh->prepare(qq{SELECT seq_region_id, name FROM seq_region});
  $sth->execute;
  $sth->bind_columns(\($seq_region_id, $chr_name)); 
  $seq_region_ids{$chr_name} = $seq_region_id while $sth->fetch;
  $sth->finish;

  return \%seq_region_ids;
}


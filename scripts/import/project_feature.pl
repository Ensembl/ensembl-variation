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

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use ImportUtils qw(load);

use DBI;
use FileHandle;
use File::Path qw(make_path);
use Getopt::Long;

$| = 1;
my $args = scalar @ARGV;
my $config = {};
GetOptions(
  $config,
  'help',                 # display help message

  'oldasm_name=s',        # name of old assembly e.g. GRCh37 
  'newasm_name=s',        # name of new assembly e.g. GRCh38

  'working_dir=s',                # working dir to store tmp results (it creates a dir for dumping features on old assembly and dir to store projection results) dir_oldasm, dir_newasm
  'feature_type|ft=s',            # run projection for variation_feature (vf) or structural_variation_feature (svf)
  'feature_table_name_oldasm|ftn_oldasm=s',     # feature table name for which to compute projections
  'feature_table_name_newasm|ftn_newasm=s',    

  'use_new_vdb_only',             
  'include_asm_exceptions',       # Include features on patches and alt loci. Else only project features on ref sequence (e.g. human 1..22,X,Y,MT)
  'load_failed_projections',
  'flip',                         # flip on the fly. This is useful if you don't want to run Variant QC for the projected features.
  'qc_ref_allele',

  # args need to be defined if feature table from old asm is not copied to new variation database
  'vdbname_oldasm|vdb_oldasm=s', # required unless use_new_vdb_only is set
  'vuser_oldasm=s',                # default ensro
  'vhost_oldasm=s',                # required unless use_neq_vdb_olny is set
  'vport_oldasm=i',                # default 3306

  'vdbname_newasm|vdb_newasm=s',   # required        
  'vuser_newasm=s',                # default ensadmin
  'vpass_newasm|p=s',              # required
  'vhost_newasm=s',                # required
  'vport_newasm=i',                # default 3306

  'cdbname_oldasm|cdb_oldasm=s',   # required
  'cuser_oldasm=s',                # default ensro
  'chost_oldasm=s',                # required
  'cport_oldasm=i',                # default 3306

  'cdbname_newasm|cdb_newasm=s',   # required
  'cuser_newasm=s',                # default ensro
  'chost_newasm=s',                # required
  'cport_newasm=i',                # default 3306

) or die "Error: Failed to parse command-line args. Try --help for usage instructions\n";

if ($config->{help} || !$args) {
  &usage;
  exit(0);
}

# vdbname_oldasm is only required if variation_feature table from vdbname_oldasm is projected
unless ($config->{vdbname_oldasm} && $config->{vhost_oldasm}) {
  unless ($config->{vdbname_newasm} && $config->{vhost_newasm}) {
    die 'The script requires information about the variation database that contains the table which needs to be projected. Try --help for usage instructions';
  } else {
    $config->{use_new_vdb_only} = 1;
  }
}


# Test argument settings are correct -------------------------
my $args_complete = 1;
foreach my $arg (qw/newasm_name oldasm_name working_dir feature_type
    feature_table_name_oldasm feature_table_name_newasm 
    vdbname_newasm vhost_newasm vpass_newasm 
    cdbname_oldasm chost_oldasm
    cdbname_newasm chost_newasm/) {
  unless ($config->{$arg}) {
    warn "--$arg must be defined\n";
    $args_complete = 0;
  }
}
unless ($args_complete) {
  die 'Missing arguments. Try --help for usage instructions';
}

unless ($config->{feature_type} eq 'vf' || $config->{feature_type} eq 'svf') {
  die '--feature_type must be vf or svf. Try --help for usage instructions';
} 

unless ($config->{use_new_vdb_only}) {
  die '--vdbname_oldasm must be defined. Try --help for usage instructions' unless ($config->{vdbname_oldasm});
  die '--vhost_oldasm must be defined. Try --help for usage instructions' unless ($config->{vhost_oldasm});
}

foreach my $arg (qw/vport_oldasm vport_newasm cport_newasm cport_oldasm/) {
  $config->{$arg} ||= 3306;
}

foreach my $arg (qw/vuser_oldasm cuser_newasm cuser_oldasm/) {
  $config->{$arg} ||= 'ensro';
}

$config->{vuser_newasm} ||= 'ensadmin';

# check that working dir exists and is empty, create all directories needed for running the script
my $working_dir = $config->{working_dir};
die "$working_dir doesn't exist. Please create $working_dir before running the script" unless (-d $working_dir);
foreach my $dir (qw/tmp_dir dir_oldasm dir_newasm/) {
  die "$dir already exists. Please delete $dir before running the script" if (-d "$working_dir/$dir");
  make_path("$working_dir/$dir") or die "Failed to create dir $working_dir/$dir $!";
  $config->{$dir} = "$working_dir/$dir";
}

# Start projection ------------------------------------------------
init_db_connections($config);
init_synonyms($config);
dump_features($config);
project_features($config);
load_features($config, 'projection');
if ($config->{load_failed_projections}) {
  load_features($config, 'no_projection');
}
sub init_db_connections {
  my $config = shift;

  unless ($config->{use_new_vdb_only}) {
    my $vdba_oldasm = new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
        -host   => $config->{vhost_oldasm},
        -user   => $config->{vuser_oldasm},
        -port   => $config->{vport_oldasm},
        -dbname => $config->{vdbname_oldasm},
        ); 
    $config->{vdba_oldasm} = $vdba_oldasm;
  }
  my $vdba_newasm = new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
      -host   => $config->{vhost_newasm},
      -user   => $config->{vuser_newasm},
      -pass   => $config->{vpass_newasm},
      -port   => $config->{vport_newasm},
      -dbname => $config->{vdbname_newasm},
      ); 
  $config->{vdba_newasm} = $vdba_newasm;

  my $cdba_oldasm = new Bio::EnsEMBL::DBSQL::DBAdaptor(
      -host   => $config->{chost_oldasm},
      -user   => $config->{cuser_oldasm},
      -port   => $config->{cport_oldasm},
      -dbname => $config->{cdbname_oldasm},
      );
  $config->{cdba_oldasm} = $cdba_oldasm;

  my $cdba_newasm = new Bio::EnsEMBL::DBSQL::DBAdaptor(
      -host   => $config->{chost_newasm},
      -user   => $config->{cuser_newasm},
      -port   => $config->{cport_newasm},
      -dbname => $config->{cdbname_newasm},
      );
  $config->{cdba_newasm} = $cdba_newasm;
}

# get synonyms for alt_loci in new assembly
sub init_synonyms {
  my $config = shift;
  my $seq_names_synonyms = {};
  if ($config->{include_asm_exceptions}) {
    my $cdba_newasm = $config->{cdba_newasm};
    my $sa_new = $cdba_newasm->get_SliceAdaptor; 
    my $aefa   = $cdba_newasm->get_AssemblyExceptionFeatureAdaptor;
    my $slices = $sa_new->fetch_all('chromosome', undef, 1, 1);

    foreach my $slice (@$slices) {
      my $seq_region_name = $slice->seq_region_name;
      if ($seq_region_name =~ /^\d+$|^X$|^Y$|^MT$/) { 
        my $assembly_exception_features = $aefa->fetch_all_by_Slice($slice);
        foreach my $feature (@$assembly_exception_features) {
          my $alt_slice = $feature->alternate_slice();
          foreach my $synonym (@{$alt_slice->get_all_synonyms()}) {
            $seq_names_synonyms->{$alt_slice->seq_region_name} = $synonym->name;
          } 
        }  
      }
    }
  }
  $config->{synonyms} = $seq_names_synonyms;
}

sub dump_features {
  my $config = shift;

  my $feature_table = $config->{feature_table_name_oldasm};

  my ($dbname, $vdba);
  if ($config->{use_new_vdb_only}) {
    $dbname = $config->{vdbname_newasm};
    $vdba   = $config->{vdba_newasm};
  } else {
    $dbname = $config->{vdbname_oldasm};
    $vdba   = $config->{vdba_oldasm};
  }

  my $dbh          = $vdba->dbc->db_handle;
  my $sth = $dbh->prepare(qq{
      SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS 
      WHERE TABLE_SCHEMA = '$dbname'
      AND TABLE_NAME = '$feature_table';
      });
  $sth->execute();

# QC that all necessary columns are there: e.g. seq_region_id, ...
  my @column_names = ();
  while (my @name = $sth->fetchrow_array) {
    push @column_names, $name[0];
  }
  $sth->finish();
  @column_names = sort @column_names;
  my $column_names_concat = join(',', @column_names); 
  $config->{sorted_column_names} = $column_names_concat;
  $sth->finish();

  my $dir_oldasm = $config->{dir_oldasm};

  my $cdba_oldasm = $config->{cdba_oldasm}; 
  my $sa_old = $cdba_oldasm->get_SliceAdaptor;

  my $slices;
  if ($config->{include_asm_exceptions}) {
    $slices = $sa_old->fetch_all('chromosome', undef, 1); 
  } else {
    $slices = $sa_old->fetch_all('chromosome', undef, 0); 
  }

  $sth = $dbh->prepare(qq{
    SELECT $column_names_concat FROM $feature_table WHERE seq_region_id = ?;
  }, {mysql_use_result => 1}); 

  foreach my $slice (@$slices) {
    my $seq_region_name = $slice->seq_region_name;
    next if ($seq_region_name =~ /PATCH/);
    my $seq_region_id = $slice->get_seq_region_id; 

    my $fh = FileHandle->new("$dir_oldasm/$seq_region_name.txt", 'w');
    $sth->execute($seq_region_id);

    while (my $row = $sth->fetchrow_arrayref) {
      my @values = map { defined $_ ? $_ : '\N' } @$row;
      my @pairs = ();
      for my $i (0..$#column_names) {
        push @pairs, "$column_names[$i]=$values[$i]";
      } 
      print $fh join("\t", @pairs), "\n";
    }
    $sth->finish();
    $fh->close(); 
  }
}

sub project_features {
  my $config = shift;
  my $dir_oldasm = $config->{dir_oldasm};    

  opendir(DIR, $dir_oldasm) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /(.+)\.txt$/) {
      my $seq_region_name = $1;     
      project_features_in_seq_region($config, $seq_region_name);
    }
  }
  closedir(DIR); 
}

sub project_features_in_seq_region {
  my $config = shift;
  my $seq_region_name = shift;

  my $feature_type = $config->{feature_type}; 
  my $synonyms     = $config->{synonyms};

  my $cdba_newasm = $config->{cdba_newasm};
  my $sa_newasm   = $cdba_newasm->get_SliceAdaptor;

  my $oldasm_name = $config->{oldasm_name};
  my $dir_oldasm  = $config->{dir_oldasm};
  my $dir_newasm  = $config->{dir_newasm};

  my $slice_newdb_oldasm = $sa_newasm->fetch_by_region(
    'chromosome',
    $seq_region_name,
    undef,
    undef,
    undef,
    $oldasm_name,
  );

  my $out = {};
  foreach my $report_type (qw/projection no_projection not_useful_projection/) {
    my $fh = FileHandle->new("$dir_newasm/$report_type\_$seq_region_name.txt", 'w');
    $out->{$report_type} = $fh;
  }
  if ($config->{flip}) {
    my $fh = FileHandle->new("$dir_newasm/flip_$seq_region_name.txt", 'w');
    $out->{flip} = $fh;
  }

  my $fh = FileHandle->new("$dir_oldasm/$seq_region_name.txt", 'r');

  while (<$fh>) {
    chomp;
    my $data = read_line($_);

    if ($feature_type eq 'vf') {
      my $map_weight = $data->{map_weight};
      next if ($map_weight > 1);
      my $start  = $data->{seq_region_start};
      my $end    = $data->{seq_region_end};
      my $strand = $data->{seq_region_strand};    
      my $is_insertion = 0;
      if ($start > $end) {
        $is_insertion = 1;
        ($start, $end) = ($end, $start);
      }
      my $feature = project_feature($config, $slice_newdb_oldasm, $start, $end, $strand);

      if ($feature) {
        my $new_start           = $feature->seq_region_start;
        my $new_end             = $feature->seq_region_end;
        my $new_strand          = $feature->seq_region_strand;
        my $new_seq_region_id   = $feature->slice->get_seq_region_id;  
        my $new_seq_region_name = $feature->slice->seq_region_name;

        # QC project from chrom to chrom, or alt_loci to alt_loci
        if ( $seq_region_name eq $new_seq_region_name || ($synonyms->{$new_seq_region_name} eq $seq_region_name) ) {
          if ($is_insertion) {
            ($new_start, $new_end) = ($new_end, $new_start);
          }
          if ($new_strand == -1 && $config->{flip}) {
            my $vf_id = $data->{variation_feature_id};
            my $var_name = $data->{variation_name};           

            my $report_flip_fh = $out->{flip};

            my $prev_allele_string = $data->{allele_string};
            reverse_comp(\($data->{allele_string}));        
            my $new_allele_string = $data->{allele_string};
            $new_strand = 1;

            print $report_flip_fh join("\t", ("VF_ID=$vf_id", "VAR_NAME=$var_name", "PREV_ALLELE_STRING=$prev_allele_string", "NEW_ALLELE_STRING=$new_allele_string", "PREV_STRAND=$strand", "NEW_STRAND=$new_strand")), "\n";
          }
          if ($config->{qc_ref_allele}) {
          }
          # update data
          $data->{seq_region_start}  = $new_start;
          $data->{seq_region_end}    = $new_end;
          $data->{seq_region_strand} = $new_strand;
          $data->{seq_region_id}     = $new_seq_region_id;

          write_output($out->{projection}, $data); 
        } else { # not a useful_projection, could be that feature from alt loci is projected to ref seq?
          write_output($out->{not_useful_projection}, $data); 
        }
      } else { # report no projection for feature
        write_output($out->{no_projection}, $data); 
      }
    } else { # project svf
      my $strand = $data->{seq_region_strand};    
      my $new_locations = {};
      my $failed_projection = 0;
      my $new_seq_region_id;
      foreach my $location (qw/outer_start seq_region_start inner_start inner_end seq_region_end outer_end/) {
        next if ($data->{$location} eq '\N');
        my $feature = project_feature($config, $slice_newdb_oldasm, $data->{$location}, $data->{$location}, $strand);
        if ($feature) {
          my $new_start  = $feature->seq_region_start;
          my $new_strand = $feature->seq_region_strand;
          $new_locations->{$location} = $new_start;
          $new_seq_region_id = $feature->slice->get_seq_region_id;
        } else {
          write_output($out->{no_projection}, $data); 
          $failed_projection = 1;
          last;
        }
      }
      # update data
      unless ($failed_projection) {
        foreach my $location (qw/outer_start seq_region_start inner_start inner_end seq_region_end outer_end/) {
          $data->{$location} = $new_locations->{$location}; 
        }
        $data->{seq_region_id} = $new_seq_region_id;
        write_output($out->{projection}, $data); 
      }
    }
  }
  $fh->close();

  foreach my $report_type (keys %$out) {
    my $fh = $out->{$report_type};
    $fh->close();
  } 
}

sub project_feature {
  my ($config, $slice, $start, $end, $strand) = @_; 

  my $newasm_name = $config->{newasm_name}; 
  my $feature = new Bio::EnsEMBL::SimpleFeature(
      -start  => $start,
      -end    => $end,
      -strand => $strand,
      -slice  => $slice,
      );
  my $projected_feature = $feature->transform('chromosome', $newasm_name);
  return $projected_feature;
}

sub read_line {
  my $line = shift;

  my @key_values = split("\t", $line);

  my $mapping = {};
  foreach my $key_value (@key_values) {
    my ($table_name, $value) = split('=', $key_value, 2);
    $mapping->{$table_name} = $value;
  }

  return $mapping;
}

sub write_output {
  my ($fh, $data) = @_;

  my @output = ();
  foreach my $column_name (sort keys %$data) {
    my $column_value = ($data->{$column_name}) ? $data->{$column_name} : '';
    push @output, $column_value;
  } 
  print $fh join("\t", @output), "\n";
}

sub load_features {
  my $config = shift;
  my $type = shift; # projection or no_projection 
    my $vdba_newasm = $config->{vdba_newasm};
  my $dbc         = $vdba_newasm->dbc;

  # create or truncate table
  my $result_table = $config->{feature_table_name_newasm};
  my $feature_table = $config->{feature_table_name_oldasm};
  if ($type eq 'no_projection') {
    $result_table = $result_table . '_failed';
  }

  $dbc->do(qq{ DROP TABLE IF EXISTS $result_table});
  $dbc->do(qq{ CREATE TABLE $result_table like $feature_table });
  my $column_names = $config->{sorted_column_names};
  my $dir_newasm = $config->{dir_newasm};
  opendir(DIR, $dir_newasm) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /^$type\_(.+)\.txt$/) {
      my $seq_region_name = $1;     
      run_cmd("cp $dir_newasm/$type\_$seq_region_name.txt $dir_newasm/load\_$type\_$seq_region_name.txt");
      my $TMP_DIR = $dir_newasm;
      my $tmp_file = "load\_$type\_$seq_region_name.txt";
      $ImportUtils::TMP_DIR = $TMP_DIR;
      $ImportUtils::TMP_FILE = $tmp_file;
      load($dbc, ($result_table, $column_names));
    }
  }
  closedir(DIR); 
}

sub run_cmd {
  my $cmd = shift;
  if (my $return_value = system($cmd)) {
    $return_value >>= 8;
    die "system($cmd) failed: $return_value";
  }
}

sub usage {
  my $usage =<<END;
Usage:
  perl project_features.pl [arguments]
Description:
    The script computes projections for all variation_features or structural_variation_features from the appropriate table.
    A workflow would look like this:
    - load data in variation_feature table, eg variation_feature_hgmd_37
    - create empty variation_feature table, eg variation_feature_hgmd_38
    - run script with all necessary parameters

    The script populates variation_feature_hgmd_38. Features that couldn't be projected are stored in
    working_dir/dir_newasm/no_projection_seq_region_name.txt.

Options:

    --help                                       # Display this message and quit

    --oldasm_name                                # Name of old assembly e.g. GRCh37
    --newasm_name                                # Name of new assembly e.g. GRCh38

    --working_dir                                # working dir to store tmp results (it creates a dir for dumping features on old assembly and dir to store projection results: dir_oldasm, dir_newasm)
    --ft | --feature_type                        # run projection for variation_feature (vf) or structural_variation_feature (svf): accepted values are vf and svf
    --ftn_oldasm | --feature_table_name_oldasm   # feature table name with feature to project
    --ftn_newasm | --feature_table_name_newasm   # feature table name to store projections
    --use_new_vdb_only                           # if this is set, feature table with features to project must be stored in vdb_newas
    --include_asm_exceptions                     # include features on patches and alt loci. Else only project features on ref sequence (e.g. human 1..22,X,Y,MT)
    --load_failed_projections                    # load failed projections into variation database 
    --flip                                       # flip on the fly. This is useful if you don't want to run Variant QC for the projected features.
    --qc_ref_allele

    --cdb_oldasm | --cdbname_oldasm              # required
    --cuser_oldasm                               # default ensro
    --chost_oldasm                               # required
    --cport_oldasm                               # default 3306

    --cdb_newasm | --cdbname_newasm              # required
    --cuser_newasm                               # default ensro
    --chost_newasm                               # required
    --cport_newasm                               # default 3306

    # Connection params for vdb_oldasm are needed, if feature table from old asm is not copied to new variation database
    --vdb_oldasm | --vdbname_oldasm              # required unless feature table that needs to be projected is in vdbname_newasm
    --vuser_oldasm                               # default ensro
    --vhost_oldasm                               # required unless feature table that needs to be projected is in vdbname_newasm
    --vport_oldasm                               # default 3306

    --vdb_newasm | --vdbname_newasm              # required
    --vuser_newasm                               # default ensadmin
    --p | --vpass_newasm                         # required
    --vhost_newasm                               # required
    --vport_newasm                               # default 3306
END

  print $usage;
}



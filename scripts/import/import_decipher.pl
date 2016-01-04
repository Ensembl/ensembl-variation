#!/usr/bin/env perl

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug create_and_load load);
use LWP::Simple;

our ($species, $input_file, $input_dir, $source_name, $TMP_DIR, $TMP_FILE, $mapping, $num_gaps,
     $target_assembly, $cs_version_number, $size_diff, $version, $registry_file, $replace_data, $debug);

GetOptions('species=s'         => \$species,
           'source_name=s'     => \$source_name,
           'input_file=s'      => \$input_file,
           'input_dir=s'       => \$input_dir,
           'tmpdir=s'          => \$ImportUtils::TMP_DIR,
           'tmpfile=s'         => \$ImportUtils::TMP_FILE,
           'mapping'           => \$mapping,
           'gaps=i'            => \$num_gaps,
           'target_assembly=s' => \$target_assembly,
           'size_diff=i'       => \$size_diff,
           'version=i'         => \$version,
           'registry=s'        => \$registry_file,
           'replace!'          => \$replace_data,
           'debug!'            => \$debug,
          );
$registry_file ||= $Bin . "/ensembl.registry";
$source_name   ||= 'DECIPHER';
$num_gaps = 1 unless defined $num_gaps;

usage('-species argument is required')    if (!$species);
usage('-input_file or -input_dir argument is required') if (!$input_file and !$input_dir);
usage('-version argument is required')    if (!$version);

my @files;
my $f_count = 1;
my $f_nb    = 1;
if ($input_dir) {
  opendir(DIR, $input_dir) or die $!;
  @files = readdir(DIR); 
  @files = sort(@files);
  close(DIR);
  $f_nb = (scalar(@files)-2);
}
else {
  @files = ($input_file);
}

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

my $add         = '';
my $study_table = "study$add";
my $sv_table    = "structural_variation$add";
my $svf_table   = "structural_variation_feature$add";
my $sva_table   = "structural_variation_association$add";
my $sv_failed   = "failed_structural_variation$add";
my $set_table   = "variation_set_structural_variation$add";
my $svs_table   = "structural_variation_sample$add";
my $pf_table    = "phenotype_feature$add";
my $pfa_table   = "phenotype_feature_attrib$add";
my $temp_table  = "temp_cnv$add";

my $tmp_sv_col  = 'tmp_class_name';
  
# connect to databases
Bio::EnsEMBL::Registry->load_all( $registry_file );
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation_private');
my $dbCore = $cdb->dbc->db_handle;
my $dbVar = $vdb->dbc->db_handle;

my $csa = Bio::EnsEMBL::Registry->get_adaptor($species, "core", "coordsystem");
our $default_cs = $csa->fetch_by_name("chromosome");


# set the target assembly
$target_assembly ||= $default_cs->version;

# get the default CS version
$cs_version_number = $target_assembly;
$cs_version_number =~ s/\D//g;

# Inheritance type
my $inheritance_attrib_type = 'inheritance_type';
my $stmt = qq{ SELECT attrib_type_id FROM attrib_type WHERE code='$inheritance_attrib_type'};
my $inheritance_attrib_type_id = ($dbVar->selectall_arrayref($stmt))->[0][0];


# run the mapping sub-routine if the data needs mapping
my (%num_mapped, %num_not_mapped, %samples, %subjects, %study_done);
my $no_mapping_needed = 0;
my $skipped = 0;
my $failed = [];
my $somatic_study;
my $fname;

my $source_id = source();



##########
#  Main  #
##########

pre_processing();

foreach my $in_file (@files) {
  next if ($in_file !~ /\.txt$/);
  
  my $msg ="File: $in_file ($f_count/$f_nb)";
  print "$msg\n";
  debug("\n##############################\n$msg\n##############################");
  
  $fname = undef;
  if ($input_dir) {
    $fname = "$input_dir/$in_file";
  }
  else {
    $fname = $in_file;
  }
  
  %num_mapped = ();
  %num_not_mapped = ();
  %samples = ();
  $no_mapping_needed = 0;
  $skipped = 0;
  $failed = [];
  $somatic_study = 0;
  

  # Parsing
  parse_input($fname);
  load_file_data();
	
	# Removing old data
  remove_data() if (defined($replace_data) && defined($source_id));
	
  # Insertion
  structural_variation();
  failed_structural_variation();
  structural_variation_feature();
  #structural_variation_association();
  #somatic_study_processing() if ($somatic_study == 1);
  #structural_variation_set() if ($var_set_id);
  structural_variation_sample();
  phenotype_feature();
  drop_tmp_table() if (!defined($debug));
  debug(localtime()." Done!\n");
  $f_count ++;
}

## Post processings ##

# Post processing for mouse annotation (delete duplicated entries in structural_variation_sample)
#post_processing_annotation() if ($species =~ /mouse|mus/i);
post_processing_feature();
#post_processing_sample();
#post_processing_phenotype();

# Finishing methods
meta_coord();
#verifications(); # URLs
cleanup() if (!defined($debug));

debug(localtime()." All done!");



#############
#  Methods  #
#############

sub load_file_data{
  debug(localtime()." Loading file into temporary table");
  
  # drop any table that may be there
  $dbVar->do("DROP TABLE IF EXISTS $temp_table;");
  
  create_and_load(
  $dbVar, "$temp_table", "id *", "type", "chr *", "outer_start i", "start i", "inner_start i",
  "inner_end i", "end i", "outer_end i", "strand i", "inheritance", "phenotype *", "sample *", "is_failed i");    
  
  
  # fix nulls
  if (($dbVar->selectrow_arrayref(qq{SELECT count(*) FROM $temp_table WHERE (outer_start=0 AND start=0 AND inner_start=0) OR (inner_end=0 AND end=0 AND outer_end=0);}))->[0] != 0) {
    warn "Structural variants with start and/or end coordinates equal to 0 found in the file!";
    $dbVar->do(qq{UPDATE $temp_table SET outer_start=outer_end, start=end, inner_start=inner_end WHERE outer_start=0 AND start=0 AND inner_start=0});
    $dbVar->do(qq{UPDATE $temp_table SET outer_end=outer_start, end=start, inner_end=inner_start WHERE outer_end=0 AND end=0 AND inner_end=0});
  }
  
  foreach my $coord('outer_start', 'inner_start', 'inner_end', 'outer_end', 'start', 'end') {
    $dbVar->do(qq{UPDATE $temp_table SET $coord = NULL WHERE $coord = 0;});
  }
  
  # Case with insertions            
  $dbVar->do(qq{UPDATE $temp_table SET start=outer_start, end=inner_start, 
                outer_end=inner_start, inner_start=NULL, inner_end=NULL
                WHERE outer_start=outer_end AND inner_start=inner_end;});                                         
}


sub source{
  debug(localtime()." Inserting into source table");
  
	my $name = $source_name;
	my $url  = 'http://decipher.sanger.ac.uk';
	my $desc = 'Database of Chromosomal Imbalance and Phenotype in Humans Using Ensembl Resources';
  # Check if the DGVa source already exists, else it create the entry
  if ($dbVar->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name='$name';})) {
    $dbVar->do(qq{UPDATE IGNORE source SET description='$desc',url='$url',version=$version where name='$name';});
  }
  else {
    $dbVar->do(qq{INSERT INTO source (name,description,url,version) VALUES ('$name','$desc','$url',$version);});
  }
  my @source_id = @{$dbVar->selectrow_arrayref(qq{SELECT source_id FROM source WHERE name='$name';})};
  return $source_id[0];
}


# Structural variations & supporting structural variations
sub structural_variation {
  
  debug(localtime()." Inserting into $sv_table table");
  
  my $stmt = qq{
    INSERT IGNORE INTO
    $sv_table (
      variation_name,
      source_id,
      class_attrib_id,
      $tmp_sv_col
    )
    SELECT
      DISTINCT
      t.id,
      $source_id,
      a1.attrib_id,
      t.type
    FROM
      $temp_table t 
      LEFT JOIN attrib a1 ON (t.type=a1.value) 
  };
  $dbVar->do($stmt);
}

  
# Failed structural variations
sub failed_structural_variation   {  
  debug(localtime()." Inserting into $sv_failed table for SVs");
  
  my $stmt = qq{
    INSERT IGNORE INTO
    $sv_failed (
      structural_variation_id,
      failed_description_id
    )
    SELECT
      sv.structural_variation_id,
      17
    FROM
      $sv_table sv,
      $temp_table t
    WHERE
      sv.variation_name=t.id AND
      (t.is_failed!=0 OR t.chr NOT IN (SELECT name FROM seq_region))
  };
  $dbVar->do($stmt);
}
  
  
# Structural variation features & supporting structural variation features
sub structural_variation_feature {
  debug(localtime()." Inserting into $svf_table table");
  
  my $stmt = qq{
    INSERT IGNORE INTO $svf_table (
      seq_region_id,
      outer_start,
      seq_region_start,
      inner_start,
      inner_end,
      seq_region_end,
      outer_end,
      seq_region_strand,
      structural_variation_id,
      variation_name,
      class_attrib_id,
      source_id
    )
    SELECT
      DISTINCT
      q.seq_region_id,
      t.outer_start,
      t.start,
      t.inner_start,
      t.inner_end,
      t.end,
      t.outer_end,
      t.strand,
      sv.structural_variation_id,
      t.id,
      sv.class_attrib_id,
      $source_id
    FROM
      seq_region q,
      $temp_table t,
      $sv_table sv
    WHERE
      q.name = t.chr AND
      t.id=sv.variation_name AND
      t.is_failed=0 AND
      sv.source_id=$source_id
  };
  $dbVar->do($stmt);
}


sub structural_variation_sample {
  my $stmt;
  
  debug(localtime()." Inserting into $svs_table table");

  # Create individual entries
  $stmt = qq{ SELECT DISTINCT sample FROM $temp_table WHERE sample NOT IN (SELECT DISTINCT name from individual)};
  my $rows_inds = $dbVar->selectall_arrayref($stmt);
  foreach my $row (@$rows_inds) {
    my $sample = $row->[0];
    my $gender = 'Unknown';
    next if ($sample eq  '');

    my $itype_val = 3;
    
    $dbVar->do(qq{ INSERT IGNORE INTO individual (name,description,gender,individual_type_id) VALUES ('$sample','Subject from DECIPHER','$gender',$itype_val)});
  }
  
  # Create sample entries
  $stmt = qq{ SELECT DISTINCT sample FROM $temp_table WHERE sample NOT IN (SELECT DISTINCT name from sample)};
  my $rows_samples = $dbVar->selectall_arrayref($stmt);
  foreach my $row (@$rows_samples) {
    my $sample = $row->[0];
    next if ($sample eq  '');

    $dbVar->do(qq{ INSERT IGNORE INTO sample (name,description,individual_id) SELECT '$sample','Sample from DECIPHER',individual_id FROM individual WHERE name='$sample'});
    }
  
  $stmt = qq{
    INSERT IGNORE INTO
    $svs_table (
      structural_variation_id,
      sample_id
    )
    SELECT DISTINCT
      sv.structural_variation_id,
      s.sample_id
    FROM
      $sv_table sv,
      $temp_table t
      LEFT JOIN sample s ON (s.name=t.sample)
    WHERE
      sv.variation_name=t.id AND
      t.sample is not null
  };
  $dbVar->do($stmt);
}


sub phenotype_feature {
  my $stmt;
  
  # Check if there is some phenotype entries
  $stmt = qq{ SELECT count(phenotype) FROM $temp_table WHERE phenotype is not null };
  if (($dbVar->selectall_arrayref($stmt))->[0][0] == 0) {
    debug(localtime()." Inserting into $pf_table table skipped: no phenotype data for this study");
    return;
  }
  
  debug(localtime()." Inserting into $pf_table table");
  
  # Create phenotype entries
  $stmt = qq{ SELECT DISTINCT phenotype FROM $temp_table 
              WHERE phenotype NOT IN (SELECT description FROM phenotype)
            };
  my $rows = $dbVar->selectall_arrayref($stmt);
  while (my $row = shift(@$rows)) {
    my $phenotype = $row->[0];
    next if ($phenotype eq '');
    
    $phenotype =~ s/'/\\'/g;
    $dbVar->do(qq{ INSERT INTO phenotype (description) VALUES ('$phenotype')});  
  }
    
  $stmt = qq{
    INSERT IGNORE INTO 
    $pf_table (
      phenotype_id,
      object_id,
      source_id,
      type,
      seq_region_id,
      seq_region_start,
      seq_region_end,
      seq_region_strand
    )
    SELECT DISTINCT
      p.phenotype_id,
      t.id,
      $source_id,
      'StructuralVariation',
      q.seq_region_id,
      t.start,
      t.end,
      t.strand
    FROM
      $temp_table t,
			seq_region q,
      phenotype p
    WHERE 
			q.name = t.chr AND
      p.description=t.phenotype
  };
  $dbVar->do($stmt);
	
	# Phenotype feature attrib
	$stmt = qq{
    INSERT IGNORE INTO 
    $pfa_table (
      phenotype_feature_id,
      attrib_type_id,
      value
    )
    SELECT DISTINCT
      pf.phenotype_feature_id,
      $inheritance_attrib_type_id,
      t.inheritance
    FROM
      $temp_table t, 
      $pf_table pf,
			phenotype p
    WHERE 
		  pf.source_id=$source_id AND
      pf.object_id=t.id AND
			pf.phenotype_id=p.phenotype_id AND
			p.description=t.phenotype
  };
  $dbVar->do($stmt);
}  


sub drop_tmp_table {
  #$dbVar->do(qq{DROP TABLE $temp_table;});
}




#### Parsing methods ####


sub parse_input {
  my $infile = shift;
  
  debug(localtime()." Parse input file: $infile");
  debug(localtime()." Start mapping features to $target_assembly if needed") if $mapping;
  
  open IN, $infile or die "Could not read from input file $infile\n";
  
  # get a slice adaptor
  my $sa = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice');
  
  open OUT, ">$TMP_DIR/$TMP_FILE" or die "Could not write to output file $TMP_DIR/$TMP_FILE\n";
  
  my $header;
  my $assembly;
  
  while(<IN>) {
    next if ($_ =~ /^#/);
		
    my $current_line = $_;
    chomp ($current_line);
    my $info = parse_line($current_line);
    
    my $is_failed = 0;
    
    
    #foreach my $info (@$infos) {
    
      ###### Mapping work ######  
      if ($mapping && $assembly !~ /$cs_version_number/) {
        
        # get a slice on the old assembly
        my $slice;
        my $start_c = $info->{start};
        my $end_c   = $info->{end};
  
        if (!$info->{start}) {
          $start_c = ($info->{outer_start}) ? $info->{outer_start} : $info->{inner_start};
        }
        if (!$info->{end}) {
          $end_c = ($info->{outer_end}) ? $info->{outer_end} : $info->{inner_end};
        }
  
        eval { $slice = $sa->fetch_by_region('chromosome', $info->{chr}, $start_c, $end_c, 1, $assembly); };
  
        # check got the slice OK
        if(!defined($slice)) {
           warn("Structural variant '".$info->{ID}."' (study ".$header->{study}."): Unable to map from assembly $assembly or unable to retrieve slice ".$info->{chr}."\:$start_c\-$end_c");
          $skipped++;
          $num_not_mapped{$assembly}++;
          $is_failed = 2;
        }
        else {
          my $count = 0;
          my ($min_start, $max_end) = (999999999, -999999999);
          my $to_chr;
    
          # project the slice to the latest assembly
          # this may produce several "segments"
          foreach my $segment(@{$slice->project('chromosome', $target_assembly)}) {
    
              # get the slice that this segment is on
            my $to_slice = $segment->to_Slice();
    
            # check this segment is on the same chrom
            if((defined $to_chr) && ($to_chr ne $to_slice->seq_region_name)) {
              $count = 0;
              warn("Segments of ".$info->{ID}." map to different chromosomes ($to_chr and ", $to_slice->seq_region_name, ")");
              $is_failed = 2;
            }
    
            # get the chromosome name
            $to_chr = $to_slice->seq_region_name;
    
            # get the slice start and end
            my $to_start = $to_slice->start;
            my $to_end = $to_slice->end;
    
            # now calculate the coords of the segment on the slice
            $to_start = ($to_start + $segment->from_start) - 1;
            $to_end = ($to_end + $segment->from_start) - 1;
    
            # store the min/max
            $min_start = $to_start if $to_start < $min_start;
            $max_end = $to_end if $to_end > $max_end;
    
            # print this segment
            #print "$id\t$chr\t$start\t$end\tsize ",($end - $start + 1), "\t\-\>\t", $to_chr, "\t", $to_start, "\t", $to_end, " size ", ($to_end - $to_start + 1), "\n";
          
            $count++;
          }
      
          #my $diff = abs(100 - (100 * (($end - $start + 1) / ($max_end - $min_start + 1))));
  
          #print "Before ",($end - $start + 1), " After ", ($max_end - $min_start + 1), " Diff $diff count $count\n";
  
          # allow gaps?
          if((($count - 1) <= $num_gaps && $count && !$is_failed)) {# && $diff < $size_diff) {
          
            $num_mapped{$assembly}++ if (!$info->{is_ssv});
        
            $info->{chr}   = $to_chr;
            $info->{start} = $min_start;
            $info->{end}   = $max_end;
        
            # calculate inner/outer coords
            $info->{outer_start} = $min_start - ($start_c - $info->{outer_start}) if $info->{outer_start} >= 1;
            $info->{inner_start} = $min_start + ($info->{inner_start} - $start_c) if $info->{inner_start} >= 1;
            $info->{inner_end}   = $max_end - ($end_c - $info->{inner_end}) if $info->{inner_end} >= 1;
            $info->{outer_end}   = $max_end + ($info->{outer_end} - $end_c) if $info->{outer_end} >= 1;

          }
          else {
            if (!$info->{is_ssv}) {
              warn ("Structural variant '".$info->{ID}."' (study ".$header->{study}.") has location '$assembly:".$info->{chr}.":$start_c\-$end_c' , which could not be re-mapped to $target_assembly. This variant will be labelled as failed");
              $num_not_mapped{$assembly}++;
            }  
            $is_failed = 1;
          }
        }
      }
    
      $info->{is_failed} = $is_failed;
               
      my $data = generate_data_row($info);           
    
      print OUT (join "\t", @{$data})."\n";
    #}
  }
  close IN;
  close OUT;
  
  if ($mapping && $assembly !~ /$cs_version_number/) {
    debug(localtime()." Finished SV mapping:\n\t\tSuccess: ".(join " ", %num_mapped)." not required $no_mapping_needed\n\t\tFailed: ".(join " ", %num_not_mapped)." (In no existing chromosome: $skipped)");
  }
}


sub parse_line {
  my $line     = shift;
  my @data     = split(/\t/, $line);
  my $info;
  
	$info->{subject}     = $data[0];
  $info->{chr}         = $data[3];
  $info->{SO_term}     = ($data[4] > 0) ? 'duplication' : 'deletion';
  $info->{start}       = $data[1];
  $info->{end}         = $data[2];
  $info->{strand}      = 1;
	$info->{inheritance} = $data[5];
	$info->{phenotype}   = $data[6];
	$info->{ID}          = "DEC".$info->{subject}."_".$info->{chr}."_".$info->{start}."_".$info->{end};

  return $info;
}


sub decode_text {
  my $text = shift;
  
  $text  =~ s/%3B/;/g;    
  $text  =~ s/%3D/=/g;
  $text  =~ s/%25/%/g;
  $text  =~ s/%26/&/g;
  $text  =~ s/%2C/,/g;
  
  return $text;
}


#### Pre processing ####
sub pre_processing {

  debug(localtime()." Prepare the structural variation tables by adding temporary keys and columns");
  debug(localtime()."\t - Add temporary unique keys in $sv_table, $svf_table and $svs_table tables");
  
  # Prepare the structural_variation table
  if ($dbVar->do(qq{SHOW KEYS FROM $sv_table WHERE Key_name='name_key';}) < 1){
    $dbVar->do(qq{ALTER TABLE $sv_table ADD CONSTRAINT UNIQUE KEY `name_key` (`variation_name`)});
  }
  
  # Prepare the structural_variation_feature table
  if ($dbVar->do(qq{SHOW KEYS FROM $svf_table WHERE Key_name='name_coord_key';}) < 1){
    $dbVar->do(qq{ALTER TABLE $svf_table ADD CONSTRAINT  UNIQUE KEY
                  `name_coord_key` (`variation_name`,`seq_region_id`,`seq_region_start`,`seq_region_end`)
                 });
  }
  if ($dbVar->do(qq{SHOW KEYS FROM $svf_table WHERE Key_name='name_key';}) < 1){
    $dbVar->do(qq{ALTER TABLE $svf_table ADD KEY `name_key` (`variation_name`)});
  }
  
  # Prepare the structural_variation_sample table
  if ($dbVar->do(qq{SHOW KEYS FROM $svs_table WHERE Key_name='sv_id_key';}) < 1){
    $dbVar->do(qq{ALTER TABLE $svs_table ADD CONSTRAINT  UNIQUE KEY `sv_id_key` (`structural_variation_id`,`sample_id`)});
  }
  
  # Prepare the individual table
  debug(localtime()."\t - Add temporary unique keys in the individual table");
  if ($dbVar->do(qq{SHOW KEYS FROM individual WHERE Key_name='name_key';}) < 1){
    $dbVar->do(qq{ALTER TABLE individual ADD KEY `name_key` (`name`)});
  }
  
  # Prepare the sample table
  debug(localtime()."\t - Add temporary unique keys in the sample table");
  if ($dbVar->do(qq{SHOW KEYS FROM sample WHERE Key_name='name_key';}) < 1){
    $dbVar->do(qq{ALTER TABLE sample ADD KEY `name_key` (`name`)});
  }
  
  debug(localtime()."\t - Add temporary columns in the $sv_table table");
  
  if ($dbVar->do(qq{show columns from $sv_table like '$tmp_sv_col';}) != 1){
    $dbVar->do(qq{ALTER TABLE $sv_table ADD COLUMN $tmp_sv_col varchar(255);});
  }
}


#### Post processing ####

sub post_processing_annotation {
  debug(localtime()." Post processing of the table $svs_table");
  
  my $stmt = qq{ DELETE FROM $svs_table WHERE sample_id IS NULL };
  $dbVar->do($stmt);
}


sub post_processing_feature {
  debug(localtime()." Post processing of the table $svf_table");
  
  my $stmt_s = qq{ UPDATE $svf_table SET outer_start=NULL, inner_start=NULL 
                   WHERE outer_start=seq_region_start AND inner_start=seq_region_start
                 };
  $dbVar->do($stmt_s);                         

  my $stmt_e = qq{ UPDATE $svf_table SET outer_end=NULL, inner_end=NULL 
                   WHERE outer_end=seq_region_end AND inner_end=seq_region_end
                 };
  $dbVar->do($stmt_e);
}


#### Finishing methods ####

sub meta_coord{
  
  debug(localtime()." Adding entries to meta_coord table");
  
  # Structural variation feature
  my $max_length_ref = $dbVar->selectall_arrayref(qq{
    SELECT MAX(seq_region_end - seq_region_start + 1) FROM $svf_table
  });
  my $max_length = $max_length_ref->[0][0];
  if (!defined($max_length)) { $max_length = 0; }
  
  my $cs = $default_cs->dbID;
  
  $dbVar->do(qq{DELETE FROM meta_coord where table_name = '$svf_table';});
  $dbVar->do(qq{INSERT INTO meta_coord(table_name, coord_system_id, max_length) VALUES ('$svf_table', $cs, $max_length);});
}


sub cleanup {
  debug(localtime()." Cleanup temporary columns and keys");

  # Drop constraints in structural_variation_feature;
  $dbVar->do(qq{ALTER TABLE $svf_table DROP KEY name_coord_key});
  $dbVar->do(qq{ALTER TABLE $svf_table DROP KEY name_key});
  debug(localtime()."\t - Table $svf_table: cleaned");

  # Drop a unique constraint in structural_variation_sample;
  $dbVar->do(qq{ALTER TABLE $svs_table DROP KEY sv_id_key});
  debug(localtime()."\t - Table $svs_table: cleaned");
  
  # Drop a unique constraint in individual;
  $dbVar->do(qq{ALTER TABLE individual DROP KEY name_key});
  debug(localtime()."\t - Table individual: cleaned");
  
  # Drop a unique constraint in sample;
  $dbVar->do(qq{ALTER TABLE sample DROP KEY name_key});
  debug(localtime()."\t - Table sample: cleaned");
  
  # structural_variation table
  
  # Drop a unique constraint in structural_variation
  $dbVar->do(qq{ALTER TABLE $sv_table DROP KEY name_key});
  
  my $sv_flag = 0;
  # Column tmp_class_name" in structural_variation
  my $sth1 = $dbVar->prepare(qq{ SELECT count(*) FROM $sv_table WHERE source_id=$source_id AND class_attrib_id=0});
  $sth1->execute();
  my $sv_count = ($sth1->fetchrow_array)[0];
  $sth1->finish;
  if ($sv_count != 0) {
    print STDERR "\tThe table $sv_table has $sv_count variants with no class_attrib_id defined!\nPlease, look at the column '$tmp_sv_col'\n";
    $sv_flag = 1;
  } 
  else {
    $dbVar->do(qq{ALTER TABLE $sv_table DROP COLUMN $tmp_sv_col});
  }
  
  debug(localtime()."\t - Table $sv_table: cleaned") if ($sv_flag == 0);

}



#### Other methods ####

sub remove_data {
  
  debug(localtime()." Remove the structural variation data of this study");

  # structural_variation_sample
  $dbVar->do(qq{ DELETE from $svs_table WHERE structural_variation_id IN 
                 (SELECT structural_variation_id FROM $sv_table WHERE source_id=$source_id);
            });
  # structural_variation_association          
  $dbVar->do(qq{ DELETE from $sva_table WHERE structural_variation_id IN 
                 (SELECT structural_variation_id FROM $sv_table WHERE source_id=$source_id)
            });  
  # structural_variation_feature                  
  $dbVar->do(qq{ DELETE from $svf_table WHERE structural_variation_id IN 
                 (SELECT structural_variation_id FROM $sv_table WHERE source_id=$source_id)
            });
  # variation_set_structural_variation          
  $dbVar->do(qq{ DELETE from $set_table WHERE structural_variation_id IN 
                 (SELECT structural_variation_id FROM $sv_table WHERE source_id=$source_id)
            });
						
	# phenotype_feature_attrib                  
  $dbVar->do(qq{ DELETE from $pfa_table WHERE phenotype_feature_id IN 
                 (SELECT phenotype_feature_id FROM $pf_table WHERE source_id=$source_id AND 
                 (type='StructuralVariation' OR type='SupportingStructuralVariation'))
            }); 					          
  # phenotype_feature                  
  $dbVar->do(qq{ DELETE from $pf_table WHERE object_id IN 
                 (SELECT variation_name FROM $sv_table WHERE source_id=$source_id) AND 
                 (type='StructuralVariation' OR type='SupportingStructuralVariation')
            });          
  # failed_structural_variation
  $dbVar->do(qq{ DELETE from $sv_failed WHERE structural_variation_id IN 
                 (SELECT structural_variation_id FROM $sv_table WHERE source_id=$source_id)
            });
  # structural_variation
  $dbVar->do(qq{ DELETE from $sv_table WHERE source_id=$source_id });                  
}


sub generate_data_row {
  my $info = shift;
 
  $info->{phenotype} = decode_text($info->{phenotype}); 
 
  my @row = ($info->{ID},
             $info->{SO_term},
             $info->{chr}, 
             $info->{outer_start}, 
             $info->{start}, 
             $info->{inner_start}, 
             $info->{inner_end}, 
             $info->{end}, 
             $info->{outer_end},
             $info->{strand}, 
             $info->{inheritance},
						 $info->{phenotype},
						 $info->{subject},
						 $info->{is_failed}
            );
  return \@row;
}


sub usage {
  die shift, qq{

Options -
  -species         : species name (required)
  -target_assembly : assembly version to map to (optional)
  -tmpdir          : (optional)
  -tmpfile         : (optional)
  -input_file      : file containing DGVa data dump (required if no input_dir)
  -input_dir       : directory containing DGVa data dump (required if no input_file)
  -mapping         : if set, the data will be mapped to $target_assembly using the Ensembl API
  -gaps            : number of gaps allowed in mapping (defaults to 1) (optional)
  -size_diff       : % difference allowed in size after mapping (optional)
  -version         : version number of the data (required)
  -registry        : registry file (optional)
  -replace         : flag to remove the existing study data from the database before import them (optional)
  -debug           : flag to keep the $temp_table table (optional)
  };
}

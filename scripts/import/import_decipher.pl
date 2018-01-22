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
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug create_and_load load);
use LWP::Simple;
use JSON;

our ($species, $input_file, $source_name, $TMP_DIR, $TMP_FILE,
     $version, $registry_file, $debug);

GetOptions('species=s'         => \$species,
           'source_name=s'     => \$source_name,
           'input_file=s'      => \$input_file,
           'tmpdir=s'          => \$ImportUtils::TMP_DIR,
           'tmpfile=s'         => \$ImportUtils::TMP_FILE,
           'version=i'         => \$version,
           'registry=s'        => \$registry_file,
           'debug!'            => \$debug,
          );
$registry_file ||= $Bin . "/ensembl.registry";
$source_name   ||= 'DECIPHER';

usage('-species argument is required')    if (!$species);
usage('-input_file argument is required') if (!$input_file);
usage('-version argument is required')    if (!$version);

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
  
my $long_variant_type = 'CNV';
  
# connect to databases
Bio::EnsEMBL::Registry->load_all( $registry_file );
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation_private');
my $dbCore = $cdb->dbc->db_handle;
my $dbVar = $vdb->dbc->db_handle;

my $csa = Bio::EnsEMBL::Registry->get_adaptor($species, "core", "coordsystem");
our $default_cs = $csa->fetch_by_name("chromosome");

# Clinical significance
my $cs_attr_type_sth = $dbVar->prepare(qq{ SELECT attrib_type_id FROM attrib_type WHERE code='clinvar_clin_sig'});
$cs_attr_type_sth->execute;
my $cs_attr_type_id = ($cs_attr_type_sth->fetchrow_array)[0];
$cs_attr_type_sth->finish;
die "No attribute type found for the entry 'clinvar_clin_sig'" if (!defined($cs_attr_type_id));
my $clin_sign_types = get_attrib_list($cs_attr_type_id);

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

%num_mapped = ();
%num_not_mapped = ();
%samples = ();
$no_mapping_needed = 0;
$skipped = 0;
$failed = [];
$somatic_study = 0;
  

# Parsing
parse_input($input_file);
load_file_data();
  
# Insertion
structural_variation();
failed_structural_variation();
structural_variation_feature();
structural_variation_sample();
phenotype_feature();
drop_tmp_table() if (!defined($debug));
debug(localtime()." Done!\n");

## Post processings ##

# Post processing (delete duplicated entries in structural_variation_sample)
post_processing_feature();

# Finishing methods
meta_coord();
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
  "inner_end i", "end i", "outer_end i", "strand i", "inheritance", "phenotype *", "sample *", "is_failed i", "clin_sign");    
  
  
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
      study_id,
      class_attrib_id,
      clinical_significance,
      $tmp_sv_col
    )
    SELECT
      DISTINCT
      t.id,
      $source_id,
      $source_id,
      a1.attrib_id,
      t.clin_sign,
      t.type
    FROM
      $temp_table t 
      LEFT JOIN attrib a1 ON (t.type=a1.value) 
  };
  $dbVar->do($stmt);
  
  my $pp_stmt =  qq{ UPDATE $sv_table SET clinical_significance=null where clinical_significance='' };
  $dbVar->do($pp_stmt);
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
      source_id,
      study_id
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
      $source_id,
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
    $dbVar->do(qq{ INSERT IGNORE INTO phenotype (description) VALUES ('$phenotype')});  
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
  
  open IN, $infile or die "Could not read from input file $infile\n";
  
  # get a slice adaptor
  my $sa = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice');
  
  open OUT, ">$TMP_DIR/$TMP_FILE" or die "Could not write to output file $TMP_DIR/$TMP_FILE\n";
  
  my $header;
  my $assembly;
  
  while(<IN>) {
    next if ($_ =~ /^track/);

    my $current_line = $_;
    chomp ($current_line);
    
    my $info = parse_line($current_line);
    
    if ($info->{'phenotype'}) {
            
      my @phenotypes = keys(%{$info->{'phenotype'}});
 
      foreach my $phe (@phenotypes) {

        $info->{'phenotype'} = $phe;
        my $data = generate_data_row($info);
        print OUT (join "\t", @{$data})."\n";
      }
    }
    else {
      my $data = generate_data_row($info);
      print OUT (join "\t", @{$data})."\n";
    }
  }
  close IN;
  close OUT;
}


sub parse_line {
  my $line     = shift;
  my @data     = split(/\t/, $line);
  my $info;
  
  $info->{chr}         = $data[0];
  $info->{start}       = $data[1]+1;
  $info->{end}         = $data[2];
  $info->{subject}     = $data[3];
  $info->{strand}      = ($data[5] eq '+') ? 1 : -1;
  
  $info->{chr} =~ s/chr//i;
  $info->{chr} = 'MT' if ($info->{chr} eq 'M');
  
  my $json_text = decode_json($data[13]);
  
  return if ($json_text->{'variant_type'} ne $long_variant_type);
  
  $info->{SO_term}     = ($json_text->{'mean_ratio'} > 0) ? 'duplication' : 'deletion';
  $info->{inheritance} = $json_text->{'inheritance'};
  $info->{clin_sign}   = lc($json_text->{'pathogenicity'}) if ($json_text->{'pathogenicity'});
  $info->{ID}          = "DEC".$info->{subject}."_".$info->{chr}."_".$info->{start}."_".$info->{end};
  
  foreach my $phe (@{$json_text->{phenotypes}}) {
    my $p_name = $phe->{name};
    $info->{phenotype}{$p_name} = 1;
  }

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
             $info->{is_failed},
             $info->{clin_sign}
            );
  return \@row;
}


sub get_attrib_list {
  my $type_id = shift;
  my %attrib_list;
  my $attr_list_sth = $dbVar->prepare(qq{ SELECT attrib_id,value FROM attrib WHERE attrib_type_id=$type_id});
  $attr_list_sth->execute;
  while(my @row = $attr_list_sth->fetchrow_array) {
   $attrib_list{$row[1]} = $row[0];
  }
  return \%attrib_list;
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
  -version         : version number of the data (required)
  -registry        : registry file (optional)
  -debug           : flag to keep the $temp_table table (optional)
  };
}

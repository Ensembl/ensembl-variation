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
use DBI;
use Socket;
use Bio::EnsEMBL::Registry;
use Getopt::Long;
use POSIX qw(strftime);
use Cwd qw(cwd);

$| = 1;
socketpair(CHILD, PARENT, AF_UNIX, SOCK_STREAM, PF_UNSPEC) or die "ERROR: Failed to open socketpair: $!";
PARENT->autoflush(1);

my $registry = 'Bio::EnsEMBL::Registry';
my $config = configure();

my $reg_file = $config->{registry}; # the reg_file needs to be added to the config.
$registry->load_all($reg_file); # loading the registry is important 
my $db_adaptor = $registry->get_DBAdaptor("homo_sapiens", "variation");
my $dbh = $db_adaptor->dbc;
my $dbname = $dbh->dbname;
debug($config, "Connected to the new database $dbname");

my $old_dbh;
my $old_reg_file = $config->{old_registry} if (defined ($config->{old_registry}));
my $release =  $config->{release} if (defined ($config->{release}));

my $old_var_file = "old_variation.txt";
my $sorted_old_var = "sorted_old_variation.txt";
my $new_var = "new_variation.txt";
my $old_vf_file = "old_variation_feature.txt";
my $sorted_old_var_feat = "sorted_old_variation_feature.txt";
my $new_var_feat = "new_variation_feature.txt";
my $sorted_new_vf = "sorted_new_variation_feature.txt";
my $old_pheno_feat_file = "old_pheno_feature.txt";
my $old_pheno_feat_attrib_file = "old_pheno_feature_attrib.txt";
my $sorted_old_pheno_feat = "sorted_old_pheno.txt";
my $new_pf_file = "new_pf.txt";
my $sorted_new_pf = "sorted_new_pf.txt";
my $sorted_old_pfa = "sorted_old_pfa.txt";
my $sorted_new_pfa = "sorted_new_pfa.txt";
my $old_tv_file = "old_tv.txt";
my $sorted_old_tv_file = "sorted_old_tv.txt";
my $new_tv_file = "new_tv.txt";
my $old_mtmp_file = "old_mtmp.txt";
my $sorted_old_mtmp_file = "sorted_old_mtmp.txt";
my $new_mtmp_file = "new_mtmp.txt";


if ($old_reg_file) {
  $registry->load_all($old_reg_file);
  my $db_adaptor = $registry->get_DBAdaptor("homo_sapiens", "variation");
  # Connect to DBI
  $old_dbh = $db_adaptor->dbc;
  my $old_dbname = $old_dbh->dbname;
  debug($config, "Connected to the old database $old_dbname");
} else {
  my $old_release = $release - 1;
  my $old_dbname = $dbh->dbname =~ s/_${release}_/_${old_release}_/gr;
  my $old_host = $dbh->host;
  my $old_port = $dbh->port;
  #Connect to DBI
  $old_dbh = DBI->connect("DBI:mysql:database=${old_dbname};host=${old_host};port=${old_port}", "ensro", "");
  debug($config, "Connected to the old database $old_dbname");
}

my $TMP_DIR = $config->{tmp};
if ($config->{test}) {
  test();
} else {
  main();
}

sub test {

  debug($config, "Running in test mode"); 
  
  my $old_dbname = $old_dbh->dbname;

  debug($config, "Dumping HGMD data from tables variation and variation_feature from the $old_dbname");
  dump_vdata_into_file($old_dbh);

  debug($config, "Dumping old phenotype data with HGMD as source from $old_dbname");
  dump_pheno($old_dbh);

  debug($config, "Sorting files based on the variation_id column");
  $TMP_DIR = $TMP_DIR . "/";
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_var} ${TMP_DIR}${old_var_file}");
  system("sort -k 6 -o ${TMP_DIR}${sorted_old_var_feat} ${TMP_DIR}${old_vf_file}");


  if ($config->{transcript}) {
    debug($config, "Sorting transcript variation files based on the variation_feature_id column");
    system("sort -k 2 -o ${TMP_DIR}${sorted_old_tv_file} ${TMP_DIR}${old_tv_file}");
    system("sort -k 2 -o ${TMP_DIR}${sorted_old_mtmp_file} ${TMP_DIR}${old_mtmp_file}");
  }

  debug($config, "Sorting files based on the phenotype_feature_id column");
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_pheno_feat} ${TMP_DIR}${old_pheno_feat_file}");
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_pfa} ${TMP_DIR}${old_pheno_feat_attrib_file}");

  debug($config, "Assigns new variation_ids and variation_feature_ids using the existing maximum variation_id and variation_feature_id");
  manipulate_var_ids($dbh);

  debug($config, "Manipulating Phenotype ids in each file");
  manipulate_pheno_ids($dbh);

  debug($config, "Sorting the new variation_feature file by seq regions. Test mode no insertion");
  system("sort -k2,2n -k3,3n -k4,4n -o ${TMP_DIR}${sorted_new_vf} ${TMP_DIR}${new_var_feat}");


  debug($config, "Sorting the new phenotype_feature file by seq regions. Test mode no insertion"); 
  system("sort -k7,7n -k8,8n -k9,9n -o ${TMP_DIR}${sorted_new_pf} ${TMP_DIR}${new_pf_file}");
}

sub main { 
  
   my $old_dbname = $old_dbh->dbname;

  debug($config, "Dumping HGMD data from tables variation and variation_feature from the $old_dbname");
  dump_vdata_into_file($old_dbh);

  debug($config, "Dumping old phenotype data with HGMD as source from $old_dbname");
  dump_pheno($old_dbh);

  debug($config, "Sorting files based on the variation_id column");
  $TMP_DIR = $TMP_DIR . "/";
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_var} ${TMP_DIR}${old_var_file}");
  system("sort -k 6 -o ${TMP_DIR}${sorted_old_var_feat} ${TMP_DIR}${old_vf_file}");
  
  if ($config->{transcript}) {
    debug($config, "Sorting transcript variation files based on the variation_feature_id column");
    system("sort -k 2 -o ${TMP_DIR}${sorted_old_tv_file} ${TMP_DIR}${old_tv_file}");
    system("sort -k 2 -o ${TMP_DIR}${sorted_old_mtmp_file} ${TMP_DIR}${old_mtmp_file}");
  }

  debug($config, "Sorting files based on the phenotype_feature_id column");
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_pheno_feat} ${TMP_DIR}${old_pheno_feat_file}");
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_pfa} ${TMP_DIR}${old_pheno_feat_attrib_file}");

  debug($config, "Assigns new variation_ids and variation_feature_ids using the maximum variation_id and variation_feature_id");
  manipulate_var_ids($dbh);

  debug($config, "Assigns new phenotype_feature id using the existing maximum phenotype_feature_id");
  manipulate_pheno_ids($dbh);

  debug($config, "Inserting HGMD data into variation table");
  insert_into_variation_table($dbh, $new_var);

  debug($config, "Sorting the new variation_feature file by seq regions");
  system("sort -k2,2n -k3,3n -k4,4n -o ${TMP_DIR}${sorted_new_vf} ${TMP_DIR}${new_var_feat}");
  debug($config, "Inserting HGMD data into variation_feature table");
  insert_variation_feature($dbh, $sorted_new_vf);

  debug($config, "Inserting HGMD data into variation_set_variation table");
  load_all_variation_sets($dbh, $new_var);

  debug($config, "Sorting the new phenotype_feature file by seq regions"); 
  system("sort -k7,7n -k8,8n -k9,9n -o ${TMP_DIR}${sorted_new_pf} ${TMP_DIR}${new_pf_file}");
  debug($config, "Inserting into phenotype_feature table");
  insert_pheno_feature($dbh, $sorted_new_pf);

  debug($config, "Inserting into phenotype_feature_attrib table");
  insert_pheno_feature_attrib($dbh, $sorted_new_pfa);

  debug($config, "Removing files");
  system("rm *.txt");

}

sub dump_vdata_into_file {
  my $old_dbhvar = shift;

  my $dump_var = $old_dbhvar->prepare(qq [ SELECT variation_id, source_id, name, class_attrib_id, somatic, evidence_attribs, display FROM variation WHERE source_id = 8 ]);
  my $dump_var_feat = $old_dbhvar->prepare(qq[ SELECT variation_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, variation_id, allele_string, ancestral_allele, variation_name, map_weight, source_id, consequence_types, variation_set_id, class_attrib_id, evidence_attribs from variation_feature where source_id = 8  ]);
 
  open(my $var, ">>$TMP_DIR/$old_var_file") or die ("Cannot open $TMP_DIR/$old_var_file: $!");
  open(my $var_fea, ">>$TMP_DIR/$old_vf_file") or die ("Cannot open $TMP_DIR/$old_vf_file: $!");

  $dump_var->execute();
  $dump_var_feat->execute();

  while ( my $aref = $dump_var->fetchrow_arrayref() ) {
    my @a =  map {defined($_) ? $_ : '\N'} @$aref;
    print $var join("\t", @a), "\n";
  }
  
  while ( my $varref = $dump_var_feat->fetchrow_arrayref() ) {
    my @var_feat = map {defined($_) ? $_ : '\N'} @$varref;
    print $var_fea join("\t", @var_feat), "\n";
  }
  
  $dump_var->finish();
  $dump_var_feat->finish();

  close $var;
  close $var_fea;

  if ($config->{transcript}) {
    my $dump_tv = $old_dbhvar->prepare(qq[ SELECT transcript_variation_id, variation_feature_id, feature_stable_id, allele_string, somatic, consequence_types, cds_start, cds_end, cdna_start, cdna_end, translation_start, translation_end, distance_to_transcript, display from transcript_variation where variation_feature_id in  (SELECT variation_feature_id from variation_feature where source_id = 8) ]);
    open(my $tv, ">>$TMP_DIR/$old_tv_file") or die ("Cannot open $TMP_DIR/$old_tv_file: $!");
    $dump_tv->execute();

    while ( my $tvref = $dump_tv->fetchrow_arrayref() ){
      my @tv = map {defined($_) ? $_: '\N' } @$tvref;
      print $tv join("\t", @tv), "\n";
    }

    $dump_tv->finish();

    my $dump_mtmp = $old_dbhvar->prepare(qq[ SELECT  variation_feature_id, feature_stable_id, allele_string, consequence_types, cds_start, cds_end, cdna_start, cdna_end, translation_start, translation_end, distance_to_transcript from MTMP_transcript_variation where variation_feature_id in  (SELECT variation_feature_id from variation_feature where source_id = 8) ]);
    open (my $mtmp, ">>$TMP_DIR/$old_mtmp_file") or die ("Cannot open $TMP_DIR/$old_mtmp_file: $!");
    $dump_mtmp->execute();

    while ( my $mtmpref = $dump_mtmp->fetchrow_arrayref() ){
      my @mtmp = map {defined($_) ? $_: '\N'} @$mtmpref;
      print $mtmp join("\t", @mtmp), "\n";
    }

    $dump_mtmp->finish();
  }
 
}

sub manipulate_var_ids {
  my $dbhvar = shift; 

  my $select_max_var = $dbhvar->prepare(qq[SELECT MAX(variation_id) from variation]);
  $select_max_var->execute();
  my $var = $select_max_var->fetchall_arrayref();
  my $max_var = $var->[0]->[0];
  $max_var = $max_var + 1;
  
  my $select_max_feat = $dbhvar->prepare(qq[SELECT MAX(variation_feature_id) from variation_feature]);
  $select_max_feat->execute();
  my $vf = $select_max_feat->fetchall_arrayref();
  my $max_var_feat = $vf->[0]->[0];
  $max_var_feat = $max_var_feat + 1;
  
  $select_max_var->finish();
  $select_max_feat->finish();
 
  open(my $fh, '<', "$TMP_DIR/$sorted_old_var") or die "Cannot open file: $!";
  open(my $new_var_file, '>', "$TMP_DIR/$new_var") or die "Cannot open file: $!";
  while ( my $line = <$fh> ) {
    chomp $line;
    my @columns = split("\t", $line);

    $columns[0] = $max_var++;
    print $new_var_file join("\t", @columns), "\n";
  }
  close $fh;
  close $new_var_file;
  
  $max_var = $var->[0]->[0];
  $max_var = $max_var + 1;

  open(my $old_var, '<', "$TMP_DIR/$sorted_old_var_feat") or die "Cannot open file: $!";
  open(my $new_vf, '>', "$TMP_DIR/$new_var_feat") or die "Cannot open file: $!";

  while (my $var_feat = <$old_var> ) {
    chomp $var_feat;
    my @columns = split ("\t", $var_feat);
    
    $columns[0] = $max_var_feat++;
    $columns[5] = $max_var++;
    print $new_vf join ("\t", @columns), "\n";
  }

  close $old_var;
  close $new_vf;
  
  if ($config->{transcript}) {

    my $select_max_tv = $dbhvar->prepare(qq[SELECT MAX(transcript_variation_id) from transcript_variation]);
    $select_max_tv->execute();
    my $tv = $select_max_tv->fetchall_arrayref();
    my $max_tv = $tv->[0]->[0];
    $max_tv = $max_tv + 1;

    open($old_var, '<', "$TMP_DIR/$sorted_old_var_feat") or die "Cannot open file: $!";
    open($new_vf, '<', "$TMP_DIR/$new_var_feat") or die "Cannot open file: $!";
    
    my @old_vid;
    my @new_vid;
    while (<$old_var>) {
      chomp;
      my $key = (split)[0];
      push @old_vid, $key;
    }
    while (<$new_vf>) {
      chomp;
      my $new_vid = (split)[0];
      push @new_vid, $new_vid;
    }
    close $old_var;
    close $new_vf;

    my %v_feat = map { $old_vid[$_] => $new_vid[$_] } 0..$#old_vid;
    
    open (my $old_tv, "<", "$TMP_DIR/$sorted_old_tv_file") or die "Cannot open file: $!";
    open (my $old_mtmp, "<", "$TMP_DIR/$sorted_old_mtmp_file") or die "Cannot open file: $!";
    open (my $new_tv_feat, ">", "$TMP_DIR/$new_tv_file") or die "Cannot open file: $!";
    open (my $new_mtmp_feat, ">", "$TMP_DIR/$new_mtmp_file") or die "Cannot open file: $!";
    
    while (my $tv_data = <$old_tv>){
      chomp $tv_data;
      my @columns = split ("\t", $tv_data);
      
      my $change_id = $columns[1];
      if (exists $v_feat{$change_id}) {
        $columns[0] = $max_tv++;
        $columns[1] = $v_feat{$change_id};
        print $new_tv_feat join ("\t", @columns), "\n";
      } 
    }
    close $old_tv;
    close $new_tv_feat;

    while (my $mtmp_data = <$old_mtmp>) { 
      chomp $mtmp_data;
      my @columns = split ("\t", $mtmp_data);
      
      my $mtmp_id = $columns[0];
      if (exists $v_feat{$mtmp_id}) {
        $columns[0] = $v_feat{$mtmp_id};
        print $new_mtmp_feat join ("\t", @columns), "\n";
      } 
    }
    close $old_mtmp;
    close $new_tv_feat;
   
    
  }

}

sub manipulate_pheno_ids {
  my $dbhvar = shift;

  my $select_max_feat = $dbhvar->prepare(qq[SELECT MAX(phenotype_feature_id) from phenotype_feature]);
  $select_max_feat->execute();
  my $pf = $select_max_feat->fetchall_arrayref();
  my $max_pf = $pf->[0]->[0];
  $max_pf = $max_pf + 1;
  
  $select_max_feat->finish();

  open(my $old_pf, '<', "$TMP_DIR/$sorted_old_pheno_feat") or die "Cannot open file: $!";
  open(my $new_pf, '>', "$TMP_DIR/$new_pf_file") or die "Cannot open file: $!";

  while ( my $pf = <$old_pf> ) {
    chomp $pf;
    my @columns = split ("\t", $pf);
    
    $columns[0] = $max_pf++;
    print $new_pf join("\t", @columns), "\n";
  }

  open(my $old_pfa, "<", "$TMP_DIR/$sorted_old_pfa" ) or die "Cannot open file: $!";
  open(my $new_pfa, ">", "$TMP_DIR/$sorted_new_pfa" ) or die "Cannot open file: $!";

  while (my $pfa = <$old_pfa> ) {
    chomp $pfa;
    my @columns = split ("\t", $pfa);

    $columns[0] = $max_pf++;
    print $new_pfa join("\t", @columns), "\n";
  } 

  close $old_pf;
  close $new_pf;
  close $new_pfa;
  close $old_pfa;  
}

sub insert_into_variation_table {
  my $dbhvar = shift;
  my $load_file = shift;

  my $sql = $dbhvar->prepare(qq[INSERT INTO variation (variation_id, source_id, name, class_attrib_id, somatic, evidence_attribs, display) VALUES (?,?,?,?,?,?,?)]);
  
  local *FH;
  open FH, "<", "$load_file" or die "Can not open $load_file $!";
  while (<FH>) {
    chomp;
    my $variation_id = (split)[0];
    my $source_id = (split)[1];
    my $name = (split)[2];
    my $class_attrib_id = (split)[3];
    my $somatic = (split)[4];
    my $evidence_attribs = (split)[5];
    my $display = (split)[6];

    $sql->execute($variation_id, $source_id, $name, $class_attrib_id, $somatic, $evidence_attribs, $display);
  }

  close FH;
  $sql->finish();

}

sub insert_variation_feature {
  my $dbhvar = shift;
  my $load_file = shift; 

  my $sql = $dbhvar->prepare(qq[INSERT INTO variation_feature (variation_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, variation_id, allele_string, ancestral_allele, variation_name, map_weight, source_id, consequence_types, variation_set_id, class_attrib_id, evidence_attribs) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)]);

  local *FH;
  open FH, "<", "$load_file" or die "Can not open $load_file $!";
  while (<FH>) {
    chomp;
    my $var_feat_id = (split)[0];
    my $seq_region_id = (split)[1];
    my $seq_region_start = (split)[2];
    my $seq_region_end = (split)[3];
    my $seq_region_strand = (split)[4];
    my $variation_id = (split)[5];
    my $allele_string = (split)[6];
    my $ancestral_allele = (split)[7];
    my $variation_name = (split)[8];
    my $map_weight = (split)[9];
    my $source_id = (split)[10];
    my $consequence_types = (split)[11];
    my $variation_set_id = (split)[12];
    my $class_attrib_id = (split)[13];
    my $evidence_attribs = (split)[14];

    $sql->execute($var_feat_id, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $variation_id, $allele_string, $ancestral_allele, $variation_name, $map_weight, $source_id, $consequence_types, $variation_set_id, $class_attrib_id, $evidence_attribs);

  }

  close FH;
  $sql->finish();
}

sub load_all_variation_sets {
  my $dbhvar = shift;
  my $load_file = shift;

  my $sql = qq{INSERT INTO variation_set_variation (variation_id, variation_set_id ) VALUES (?, ?)};
  my $sth = $dbhvar->prepare($sql);

  open FH, "<", "$load_file" or die "Can not open $load_file: $!";
  while (<FH>) {
    chomp;
    my $var_id = (split)[0];
    my $var_set_id = 29;

    $sth->execute($var_id, $var_set_id);
  }
  close FH;
  $sth->finish();

}

sub insert_pheno_feature {
  my $dbhvar = shift; 
  my $load_file = shift;

  my $insert_pheno = $dbhvar->prepare(qq[ INSERT INTO phenotype_feature (phenotype_feature_id, phenotype_id, source_id, type, object_id, is_significant, seq_region_id, seq_region_start, seq_region_end, seq_region_strand) VALUES (?,?,?,?,?,?,?,?,?,?) ] );

  local *FH;
  open FH, "<$load_file" or die "Can not open $load_file $!";
  while (<FH>) {
    chomp;
    my $phenotype_feature_id = (split)[0];
    my $phenotype_id = (split)[1];
    my $source_id = (split)[2];
    my $type = (split)[3];
    my $object_id = (split)[4];
    my $is_significant = (split)[5];
    my $seq_region_id = (split)[6];
    my $seq_region_start = (split)[7];
    my $seq_region_end = (split)[8];
    my $seq_region_strand = (split)[9];

    $insert_pheno->execute($phenotype_feature_id, $phenotype_id, $source_id, $type, $object_id, $is_significant, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand);
  }
  close FH;
  $insert_pheno->finish();

}

sub insert_pheno_feature_attrib {
  my $dbhvar = shift;
  my $load_file = shift;

  my $insert_pfa = $dbhvar->prepare(qq{INSERT INTO phenotype_feature_attrib (phenotype_feature_id, attrib_type_id, value) VALUES (?,?,?)});

  local *FH;
  open FH, "<$load_file" or die "Can not open $load_file $!";
  while (<FH>) {
    chomp;
    my $phenotype_feature_id = (split)[0];
    my $attrib_type_id = (split)[1];
    my $value = (split)[2];

    $insert_pfa->execute($phenotype_feature_id, $attrib_type_id, $value);
  }
  close FH;
  $insert_pfa->finish();
  
}

sub dump_pheno {
  my $old_dbhvar = shift;

  my $dump_pheno = $old_dbhvar->prepare(qq[SELECT phenotype_feature_id, phenotype_id, source_id, type, object_id, is_significant, seq_region_id, seq_region_start, seq_region_end, seq_region_strand from phenotype_feature where source_id = 8 ]);
  my $dump_pheno_attrib = $old_dbhvar->prepare(qq[SELECT * from phenotype_feature_attrib where phenotype_feature_id in (SELECT phenotype_feature_id from phenotype_feature where source_id = 8)]);

  open(my $pf, ">>$TMP_DIR/$old_pheno_feat_file") or die ("Cannot open $TMP_DIR/$old_pheno_feat_file: $!");
  open(my $pfa, ">>$TMP_DIR/$old_pheno_feat_attrib_file") or die ("Cannot open $TMP_DIR/$old_pheno_feat_attrib_file: $!");

  $dump_pheno->execute();
  $dump_pheno_attrib->execute();

  while (my $a = $dump_pheno->fetchrow_arrayref()) {
    my @p_f =  map {defined($_) ? $_ : '\N'} @$a;
    print $pf join("\t", @p_f), "\n";
  }

  while (my $apf = $dump_pheno_attrib->fetchrow_arrayref()) {
    my @pfa = map {defined($_) ? $_ : '\N'} @$apf;
    print $pfa join("\t", @pfa), "\n";
  }

  $dump_pheno->finish();
  $dump_pheno_attrib->finish();

  close $pf;
  close $pfa;

}

sub usage {

  die "\n\tUsage: import_oldhgmd.pl -registry [registry file] -release [release number] \tOptional:  -tmp [temp folder] or gets set based on current directory 
   -test to use in test mode, does not insert into database, only dumps files and sorts them
   -transcript to use only if to load into transcript and MTMP_transcript_variation
  -old_registry [old database registry file] use only if release is not defined
  --release or --old_registry needs to be defined
  NB - if using release, old database needs to be on the same server as new release. \n";
}

sub configure {
  my $config = {};
  my $args = scalar @ARGV;

  GetOptions (
    $config,
    "registry=s",
    "old_registry=s",
    "release=s",
    "tmp=s",
    "test|t",
    "transcript|tr",
    "help|h",
  ) or die "ERROR: Failed to parse command line arguments - check the documentation \n";

    if ($config->{help} || !$args){
    usage();
  }
  if (!defined($config->{tmp})) {
    $config->{tmp} = cwd;
  }
  if ( !$config->{release} && !$config->{old_registry} ) {
    usage();
  }
  return $config;

}

sub debug {
  my $config = shift;

  my $text = (@_ ? (join "", @_) : "No message");
  my $date_time = strftime("%Y-%m-%d %H:%M:%S", localtime());

  print $date_time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}
#!/usr/bin/env perl
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
if $config->{test} {
  test();
} else {
  main();
}
sub test {
  debug($config, "Running in test mode"); 

  debug($config, "Dumping old variation tables with HGMD as source from old database");
  dump_vdata_into_file($old_dbh);

  debug($config, "Dumping old phenotype with HGMD as source from old database");
  dump_pheno($old_dbh);

  debug($config, "Sorting files based on the variation_id column");
  $TMP_DIR = $TMP_DIR . "/";
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_var} ${TMP_DIR}${old_var_file}");
  system("sort -k 6 -o ${TMP_DIR}${sorted_old_var_feat} ${TMP_DIR}${old_var_far_file}");

  debug($config, "Sorting files based on the phenotype_feature_id column");
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_pheno_feat} ${TMP_DIR}${$old_pheno_feat_file}");
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_pfa} ${TMP_DIR}${old_pheno_feat_attrib_file}");
  
  debug($config, "Manipulating variation ids in each file");
  manipulate_var_ids($dbh);

  debug($config, "Manipulating Phenotype ids in each file");
  manipulate_pheno_ids($dbh);

  debug($config, "Inserting into variation_feature table but first sort, this is test mode so no insertion");
  system("sort -k 2,2 -k3,3 -k4,4 -o ${TMP_DIR}${sorted_new_vf} ${TMP_DIR}${new_var_feat}");

  debug($config, "Inserting into phenotype_feature but first sort, test mode so no insertion");
  system("sort -k7,7 -k8,8 -k9,9 -o ${TMP_DIR}${$sorted_new_pf} ${TMP_DIR}${$new_pf_file}");

}

sub main { 
  debug($config, "Dumping old variation tables with HGMD as source from old database");
  dump_vdata_into_file($old_dbh);

  debug($config, "Dumping old phenotype with HGMD as source from old database");
  dump_pheno($old_dbh);

  debug($config, "Sorting files based on the variation_id column");
  $TMP_DIR = $TMP_DIR . "/";
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_var} ${TMP_DIR}${old_var_file}");
  system("sort -k 6 -o ${TMP_DIR}${sorted_old_var_feat} ${TMP_DIR}${old_var_far_file}");

  debug($config, "Sorting files based on the phenotype_feature_id column");
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_pheno_feat} ${TMP_DIR}${$old_pheno_feat_file}");
  system("sort -k 1 -o ${TMP_DIR}${sorted_old_pfa} ${TMP_DIR}${old_pheno_feat_attrib_file}");

  debug($config, "Manipulating variation ids in each file");
  manipulate_var_ids($dbh);

  debug($config, "Manipulating Phenotype ids in each file");
  manipulate_pheno_ids($dbh);

  debug($config, "Inserting into variation table");
  insert_into_variation_table($dbh, $new_var);

  debug($config, "Inserting into variation_feature table but first sort");
  system("sort -k 2,2 -k3,3 -k4,4 -o ${TMP_DIR}${sorted_new_vf} ${TMP_DIR}${new_var_feat}");
  insert_variation_feature($dbh, $sorted_new_vf);

  debug($config, "Inserting into variation_set_variation table");
  load_all_variation_sets($dbh, $new_var);

  debug($config, "Inserting into phenotype_feature but first sort");
  system("sort -k7,7 -k8,8 -k9,9 -o ${TMP_DIR}${$sorted_new_pf} ${TMP_DIR}${$new_pf_file}");
  insert_pheno_feature($dbh, $sorted_new_pf);

  debug($config, "Inserting into phenotype_feature attrib");
  insert_pheno_feature_attrib($dbh, $sorted_old_pfa);

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
  open(my $new_var_feat, '>', "$TMP_DIR/$new_var_feat") or die "Cannot open file: $!";

  while (my $var_feat = <$old_var> ) {
    chomp $var_feat;
    my @columns = split ("\t", $var_feat);
    
    $columns[0] = $max_var_feat++;
    $columns[5] = $max_var++;
    print $new_var_feat join ("\t", @columns), "\n";
  }

  close $old_var;
  close $new_var_file;
  
  $select_max_var->finish();
  $select_max_feat->finish();
}

sub manipulate_pheno_ids {
  my $dbhvar = shift;

  my $select_max_feat = $dbhvar->prepare(qq[SELECT MAX(variation_feature_id) from phenotype_feature]);
  $select_max_feat->execute();
  my $pf = $select_max_feat->fetchall_arrayref();
  my $max_pf = $pf->[0]->[0];
  $max_pf = $max_pf + 1;
  
  $select_max_feat->finish();

  open(my $old_pf, '<', "$TMP_DIR/$sorted_old_pheno_feat") or die "Cannot open file: $!";
  open(my $new_pf, '>', "$TMP_DIR/$new_pf_file ") or die "Cannot open file: $!";

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

sub insert_into_variation_table  {
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

  my $insert_pheno = $dbhvar->prepare(qq[INSERT INTO phenotype_feature (phenotype_feature_id, phenotype_id, source_id, type, object_id, is_significant, seq_region_id, seq_region_start, seq_region_end, seq_region_strand) VALUES (?,?,?,?,?,?,?,?,?,?)]);

  local *FH;
  open FH, ">$load_file" or die "Can not open $load_file $!";
  while (<FH>) {
    chomp;
    #my $phenotype_feature_id = (split)[0];
    #my $phenotype_id = (split)[1];
    #my $source_id = (split)[2];
    #my $type = (split)[3];
    #my $object_id = (split)[4];
    #my $is_significant = (split)[5];
    #my $seq_region_id = (split)[6];
    #my $seq_region_start = (split)[7];
    #my $seq_region_end = (split)[8];
    #my $seq_region_strand = (split)[9];

    $insert_pheno->execute((split)[0], (split)[1], (split)[2], (split)[3], (split)[4], (split)[5], (split)[6], (split)[7], (split)[8], (split)[9]);
  }
  close FH;
  $insert_pheno->finish();

}

sub insert_pheno_feature_attrib {
  my $dbhvar = shift;
  my $load_file = shift;

  my $insert_pfa = $dbhvar->prepare(qq{INSERT INTO phenotype_feature_attrib (phenotype_feature_id, attrib_type_id, value) VALUES (?,?,?)});

  local *FH;
  open FH, ">$load_file" or die "Can not open $load_file $!";
  while (<FH>) {
    chomp;
    $insert_pfa->execute((split)[0], (split)[1], (split)[2]);
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
  -old_registry [old database registry file] use only if release is not defined
  --release or --old_registry needs to be defined
  NB - if using release, old database needs to be on the same server. \n";
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
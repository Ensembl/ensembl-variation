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
my $sorted_old_var = "old_new_variation.txt";
my $old_var_far_file = "old_variation_feature.txt";
my $sorted_old_var = "old_new_variation_feature.txt";
#my $old_var_syn_file = "old_variation_synonym.txt";
my $old_var_set_file = "old_variation_set.txt";
my $sorted_old_var_set = "old_new_variation_set.txt";
#my $failed_var_file = "old_failed_variation.txt"; - does not exist for 37 and 38
my $old_pheno_feat_file = "old_pheno_feature.txt";
my $old_pheno_feat_attrib_file = "old_pheno_feature_attrib.txt";


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

debug($config, "Creating temporary tables for insertion");
creating_temp_table($dbh);

debug($config, "Selecting minimum and maximum variation feature associated with HGMD from the old database");
my $sth = $old_dbh->prepare(qq[ SELECT MIN(variation_id), MAX(variation_id) FROM variation_feature where source_id = 8 ] );
$sth->execute();
my $vf = $sth->fetchall_arrayref();
my $min_id = $vf->[0]->[0];
my $chunk = 1000000;
my $max_id = $vf->[0]->[1];


my $TMP_DIR = $config->{tmp};

debug($config, "Dumping old variation tables with HGMD as source from old database");
dump_vdata_into_file($old_dbh);


debug($config, "Dumping old phenotype features with HGMD from old database");
dump_pheno($old_dbh);

debug($config, "Sorting files based on the variation_id column");
$TMP_DIR = $TMP_DIR . "/";
system("sort -k 1 -o ${TMP_DIR}${sorted_old_var} ${TMP_DIR}${old_var_file}");
system("sort -k 6 -o ${TMP_DIR}${sorted_old_var} ${TMP_DIR}${old_var_far_file}");
system("sort -k 1 -o ${TMP_DIR}${sorted_old_var_set} ${TMP_DIR}${old_var_set_file}");

sub creating_temp_table {
  my $dbhvar = shift; 

  #my $failed_var_sth = $dbhvar->prepare(qq[ CREATE TABLE IF NOT EXISTS TMP_failed_variation LIKE failed_variation]);
  my $var_sth = $dbhvar->prepare(qq[ CREATE TABLE IF NOT EXISTS TMP_variation LIKE variation]);
  #my $var_syn_sth = $dbhvar->prepare(qq[CREATE TABLE IF NOT EXISTS TMP_variation_synonym like variation_synonym]); does not exist in the table 
  my $var_set_sth = $dbhvar->prepare(qq[ CREATE TABLE IF NOT EXISTS TMP_variation_set_variation like variation_set_variation]);
  my $var_feat_sth = $dbhvar->prepare(qq[CREATE TABLE IF NOT EXISTS TMP_variation_feature like variation_feature]);
  my $pheno_feat_sth = $dbhvar->prepare(qq[CREATE TABLE IF NOT EXISTS TMP_phenotype_feature LIKE phenotype_feature]);
  my $pheno_feat_attrib_sth = $dbhvar->prepare(qq[CREATE TABLE IF NOT EXISTS TMP_phenotype_feature_attrib like phenotype_feature_attrib]);
  
  #$failed_var_sth->execute();
  $var_sth->execute();
  #$var_syn_sth->execute();
  $var_set_sth->execute();
  $var_feat_sth->execute();
  $pheno_feat_sth->execute();
  $pheno_feat_attrib_sth->execute();
}

sub dump_vdata_into_file {
  my $old_dbhvar = shift;

  my $dump_var = $old_dbhvar->prepare(qq [ SELECT variation_id, source_id, name, class_attrib_id, somatic, evidence_attribs, display FROM variation WHERE source_id = 8 ]);
  my $dump_var_feat = $old_dbhvar->prepare(qq[ SELECT variation_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, variation_id, allele_string, ancestral_allele, variation_name, map_weight, source_id, consequence_types, variation_set_id, class_attrib_id, evidence_attribs from variation_feature where source_id = 8  ]);
  #my $dump_failed_var = $old_dbhvar->prepare(qq [SELECT * from failed_variation where variation_id in (SELECT variation_id from variation where source_id = 8) ]);
  my $dump_var_set = $old_dbhvar->prepare(qq [SELECT * from variation_set_variation where  variation_id in (SELECT variation_id from variation where source_id = 8)  ]);
 
  open(my $var, ">>$TMP_DIR/$old_var_file") or die ("Cannot open $TMP_DIR/$old_var_file: $!");
  open(my $var_fea, ">>$TMP_DIR/$old_var_far_file") or die ("Cannot open $TMP_DIR/$old_var_far_file: $!");
  #open(my $failed_var, ">>$TMP_DIR/$failed_var_file ") or die ("Cannot open $TMP_DIR/$failed_var_file : $!");
  open(my $var_set, ">>$TMP_DIR/$old_var_set_file") or die ("Cannot open $TMP_DIR/$old_var_set_file: $!");

  $dump_var->execute();
  $dump_var_feat->execute();
  #$dump_failed_var->execute();
  $dump_var_set->execute();

  while ( my $aref = $dump_var->fetchrow_arrayref() ) {
    my @a =  map {defined($_) ? $_ : '\N'} @$aref;
    print $var join("\t", @a), "\n";
  }
  
  while ( my $varref = $dump_var_feat->fetchrow_arrayref() ) {
    my @var_feat = map {defined($_) ? $_ : '\N'} @$varref;
    print $var_fea join("\t", @var_feat), "\n";
  }

  #while ( my $fail_feat = $dump_failed_var->fetchrow_arrayref() ) {
   # my @fail_var = map {defined($_) ? $_ : '\N'} @$fail_feat;
    #print $failed_var join("\t", @fail_var), "\n";
  #}
  
  while ( my $v_set = $dump_var_set->fetchrow_arrayref() ) {
    my @set = map {defined($_) ? $_ : '\N'} @$v_set;
    print $var_set join("\t", @set), "\n";
  }
  
  $dump_var->finish();
  $dump_var_feat->finish();
  #$dump_failed_var->finish();
  $dump_var_set->finish();

  close $var;
  close $var_fea;
  #close $failed_var;
  close $var_set;


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


sub insert_into_table {
  my $dbhvar = shift; 
  

}

sub usage {

  die "\n\tUsage: update_new_sets.pl -registry [registry file] -release [release number] \tOptional:  -tmp [temp folder] or gets set based on current directory 
  -old_registry [old database registry file] use only if release is not defined
  --release or --old_registry needs to be defined \n";
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
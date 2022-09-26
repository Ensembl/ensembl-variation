#!/usr/bin/perl

## This script should get rs_names + variation_set_ids from last release
## and "lift" all variation_id / variation_set_id to new dbSNP import.

use strict;
use DBI;

## Connect to Databases

my $dbname = "";
my $old_dbname = "";

my $host = "mysql-ens-var-prod-3";
my $port = 4606;
my $user = "ensadmin";
my $pass = "";

# Set TMP_DIR + TMP_FILE
my $TMP_DIR = "<PATH-TO>/tmp";
my $TMP_FILE = "variation_feature.txt";

my $tmp_vset_table = "variation_set_variation_e109";

# define variants and move through the list
my $dbh = DBI->connect("DBI:mysql:database=${dbname};host=${host};port=${port}",
                       $user, $pass);

my $old_dbh = DBI->connect("DBI:mysql:database=${old_dbname};host=${host};port=${port}",
                       $user, $pass);


## Stable variation_set.variation_set_id's values
my @stable = qw(
    3, 4, 5, 6, 7, 8, 9, 10, 
    11, 12, 13, 14, 15, 16, 17, 
    27, 40, 43, 44, 45, 47, 48, 
    49, 50, 51, 52, 53, 54, 55, 56, 57
  );

### Get old rs_names / variation_set_ids
## It will create a tmp file with variation_feature table values

sub dumpPreparedSQL {
  my $dbVar = shift;
  my $tmp_num = shift;
  my $chunk = shift;
  my $size = shift;

  my $start = $chunk * $tmp_num;
  my $end = $chunk + $start;
  $end = $end < $size ? $end : $size;


  my $sql = qq{SELECT variation_name, variation_set_id FROM variation_feature WHERE variation_feature_id > $start AND variation_feature_id <= $end};

  my $sth = $dbVar->prepare( $sql );

  local *FH;
  open ( FH, ">>$TMP_DIR/$TMP_FILE" )
      or die( "Cannot open $TMP_DIR/$TMP_FILE: $!" );

  $sth->execute();

  while ( my $aref = $sth->fetchrow_arrayref() ) {
    my @a = map {defined($_) ? $_ : '\N'} @$aref;
    print FH join("\t", @a), "\n";

  }

  close FH;
  $sth->finish();

  return $end, $tmp_num;
}

# Get Variation_feature size
my $sql = qq{SELECT count(*) FROM variation_feature};
my $sth = $old_dbh->prepare( $sql );
$sth->execute();
my $size = $sth->fetchall_arrayref()->[0]->[0];
my $chunk = 1000000; # Buffer_size number to avoid memory issues

for my $tmp_num (map { $_ } 0 .. $size/$chunk) {
  dumpPreparedSQL($old_dbh, $tmp_num, $chunk, $size);
}

system("awk '{if (\$2) print \$0;}' $TMP_DIR/$TMP_FILE > $TMP_DIR/$TMP_FILE.not_empty");

### Dump new variation_set_variation / Update new variation_feature.variation_set_id

my $var_ext_sth = $dbh->prepare(qq[ SELECT variation_id FROM variation WHERE name = ? limit 1]);
my $syn_ext_sth = $dbh->prepare(qq[ SELECT variation_id FROM variation_synonym WHERE name= ? limit 1]);
my $vfid_ext_sth = $dbh->prepare(qq[ SELECT variation_feature_id, variation_set_id from variation_feature WHERE variation_id = ? limit 1]);

my $vsv_create_sth = $dbh->prepare(qq[ CREATE TABLE IF NOT EXISTS $tmp_vset_table LIKE variation_set_variation ]);
$vsv_create_sth->execute();

my $vsv_ins_sth = $dbh->prepare(qq[ INSERT IGNORE INTO $tmp_vset_table (variation_id, variation_set_id) VALUES (?,?)]);
my $vf_upd_sth = $dbh->prepare(qq[ UPDATE variation_feature SET variation_set_id = ? WHERE variation_feature_id = ?]);

open my $list, "$TMP_DIR/$TMP_FILE.not_empty" || die "Failed to open var list $ARGV[0] : $!\n"; 
while(<$list>){
  next unless /rs/; 
  chomp;

  my $rs  = (split)[0];
  my $sets = (split)[1];

  $var_ext_sth->execute($rs)||die;
  my $id =   $var_ext_sth->fetchall_arrayref();

  unless (defined $id->[0]->[0]){
    $syn_ext_sth->execute($rs)||die;
    $id =   $syn_ext_sth->fetchall_arrayref();
  }

  unless (defined $id->[0]->[0]){
    print "Skipping $rs\t$sets\n";
    next;
  }

  ## Insert into variation_feature
  $vfid_ext_sth->execute($id->[0]->[0])||die;
  my $vf = $vfid_ext_sth->fetchall_arrayref();
  my $vf_id = $vf->[0]->[0];
  my $vset_id = $vf->[0]->[1];

  if ($vset_id ne ""){
    warn "Skipping vf_id $vf_id with already filled sets $vset_id\n";
    next;
  };

  warn "adding sets $sets to variation_feature_id $vf_id\n";
  $vf_upd_sth->execute($sets, $vf_id);

  ## Insert into variation_set_variation
  my @sets = split/\,/, $sets;
  foreach my $set (@sets){
    next unless (grep { /$set/ } @stable);
    warn "adding $rs  $id->[0]->[0] to set $set\n";
    $vsv_ins_sth->execute( $id->[0]->[0], $set ) ||die "Error adding $set to $rs: $!\n";
  }
}

# Remove temp files
system("rm $TMP_DIR/$TMP_FILE $TMP_DIR/$TMP_FILE.not_empty");


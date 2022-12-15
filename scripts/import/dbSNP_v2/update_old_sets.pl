#!/usr/bin/perl

## This script should get rs_names + variation_set_ids from last release
## and "lift" all variation_id / variation_set_id to new dbSNP import.

use strict;
use DBI;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

## Connect to Databases

my ($reg_file, $old_reg_file, $release, $tmp, $help);

GetOptions ("registry=s"        => \$reg_file,
            "old_registry:s"        => \$old_reg_file,
            "release=s"      =>  \$release,
            "tmp=s"         => \$tmp,
            "help|h"           => \$help,
);

usage() unless defined $reg_file && defined $release && defined $tmp;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($reg_file);

my $db_adaptor = $registry->get_DBAdaptor("homo_sapiens", "variation");
my $dbh = $db_adaptor->dbc;

my $TMP_DIR = $tmp;
my $TMP_FILE = "variation_feature.txt";

my $tmp_vset_table = "variation_set_variation_e${release}";

my $old_dbh;

# define variants and move through the list
if ($old_reg_file) {
  $registry->load_all($old_reg_file);
  $db_adaptor = $registry->get_DBAdaptor("homo_sapiens", "variation");
  # Connect to DBI
  $old_dbh = $db_adaptor->dbc;
} else {
  my $old_release = $release - 1;
  my $old_dbname = $dbh->dbname =~ s/_${release}_/_${old_release}_/gr;
  my $old_host = $dbh->host;
  my $old_port = $dbh->port;
  # Connect to DBI
  $old_dbh = DBI->connect("DBI:mysql:database=${old_dbname};host=${old_host};port=${old_port}", "ensro", "");
}

## Stable variation_set.variation_set_id's values
my @stable = qw(
    3, 4, 5, 6, 7, 8, 9, 10, 
    11, 12, 13, 14, 15, 16, 17, 
    27, 40, 43, 44, 45, 47, 48, 
    49, 50, 51, 52, 53, 54, 55, 56, 57
  );

# Get Variation_feature size
my $sql = qq{ SELECT MIN(variation_feature_id), MAX(variation_feature_id) FROM variation_feature };
my $sth = $old_dbh->prepare( $sql );
$sth->execute();
my $vf = $sth->fetchall_arrayref();
my $min_id = $vf->[0]->[0];
my $max_id = $vf->[0]->[1];
my $chunk = 1000000; # Buffer_size number to avoid memory issues

for my $tmp_num (map { $_ } $min_id/$chunk .. $max_id/$chunk) {
  dumpPreparedSQL($old_dbh, $tmp_num, $chunk, $max_id);
}

system("awk '{if (\$2) print \$0;}' $TMP_DIR/$TMP_FILE | sort -u > $TMP_DIR/$TMP_FILE.not_empty");

### Dump new variation_set_variation / Update new variation_feature.variation_set_id

my $var_ext_sth = $dbh->prepare(qq[ SELECT variation_id FROM variation WHERE name = ? limit 1]);
my $vfid_ext_sth = $dbh->prepare(qq[ SELECT variation_set_id from variation_feature WHERE variation_id = ? limit 1]);

my $vsv_create_sth = $dbh->prepare(qq[ CREATE TABLE IF NOT EXISTS $tmp_vset_table LIKE variation_set_variation ]);
my $vsv_dis_sth = $dbh->prepare(qq[ ALTER TABLE $tmp_vset_table DISABLE KEYS ]);
$vsv_create_sth->execute();
$vsv_dis_sth->execute();


local *FH;
open ( FH, ">>$TMP_DIR/$TMP_FILE.dump" )
or die( "Cannot open $TMP_DIR/$TMP_FILE.dump: $!" );


open my $list, "$TMP_DIR/$TMP_FILE.not_empty" || die "Failed to open var list $ARGV[0] : $!\n"; 
while(<$list>){
  next unless /^rs/; 
  chomp;

  my $rs  = (split)[0];
  my $sets = (split)[1];

  $var_ext_sth->execute($rs) or die;
  my $id =   $var_ext_sth->fetchall_arrayref();

  next unless (defined $id->[0]->[0]);

  ## Check if variation_feature already has variation_set filled
  $vfid_ext_sth->execute($id->[0]->[0]) or die;
  my $vf = $vfid_ext_sth->fetchall_arrayref();
  next if $vf->[0]->[0] ne "";

  ## Insert into variation_set_variation
  my @sets = split/\,/, $sets;
  foreach my $set (@sets){
    next unless (grep { /$set/ } @stable);
    print FH $id->[0]->[0] . "\t" . $set . "\n";
  }

}

close FH;

# Dump all variation_set_variation
my $vsv_ins_sth = $dbh->prepare(qq[ LOAD DATA LOCAL INFILE \"${TMP_DIR}/${TMP_FILE}.dump\" INTO TABLE ${tmp_vset_table} ]);
$vsv_ins_sth->execute() || die "Error dumping local rsid / set file\n";

# Re-enable keys
my $vsv_enable_sth = $dbh->prepare(qq[ ALTER TABLE $tmp_vset_table ENABLE KEYS ]);
$vsv_enable_sth->execute();

# Update variation_feature
my $temp_table = 'variation_feature_with_vs';

$dbh->do(qq{DROP TABLE IF EXISTS $temp_table})
    or die "Failed to drop pre-existing temp table";
    
$dbh->do(qq{CREATE TABLE $temp_table LIKE variation_feature})
    or die "Failed to create temp table";

## remove unneccessary non-null columns (for EGenomes)

$dbh->do(qq{ALTER TABLE $temp_table  drop seq_region_id,
                                     drop seq_region_start,    drop seq_region_end,
                                     drop seq_region_strand,   drop source_id,
                                     drop map_weight })
    or die "Failed to alter temp table";

$dbh->do(qq[
    INSERT INTO $temp_table (variation_id, variation_set_id)
    SELECT variation_id, GROUP_CONCAT(DISTINCT(variation_set_id))
    FROM $tmp_vset_table
    GROUP BY variation_id
  ]);


$dbh->do(qq{
    UPDATE  variation_feature vf, $temp_table tvf
    SET     vf.variation_set_id = tvf.variation_set_id
    WHERE   vf.variation_id = tvf.variation_id
}) or die "Failed to update vf table";

$dbh->do(qq{DROP TABLE $temp_table})
    or die "Failed to drop temp table";

# Remove temp files
system("rm $TMP_DIR/$TMP_FILE $TMP_DIR/$TMP_FILE.not_empty $TMP_DIR/$TMP_FILE.dump");

# Rename table
$dbh->do(qq{ DROP TABLE IF EXISTS variation_set_variation }) or die "Failed to drop variation_set_variation table";
$dbh->do(qq{ ALTER TABLE $tmp_vset_table RENAME TO variation_set_variation }) or die "Failed to rename table";

# Add Failed to variation set
$dbh->do(qq{ 
  INSERT IGNORE INTO variation_set_variation (variation_id, variation_set_id)
  SELECT DISTINCT variation_id, 1 
  FROM failed_variation; 
}) or die "Failed to add failed to variation_set_variation table";

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

sub usage {

  die "\n\tUsage: update_old_sets.pl -registry [registry file] -release [release number] -tmp [temp folder]
\tOptional: -old_registry [old database registry file]\n\n";
}
use strict;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

## Connect to Databases

my ($reg_file, $table_name, $columns, $TMP_DIR, $help);

GetOptions ("registry=s"      => \$reg_file,
            "table_name=s"    => \$table_name,
            "tmp=s"           => \$TMP_DIR,
            "columns=s"       => \$columns,
            "help|h"          => \$help,
);

usage() unless defined $reg_file && defined $table_name && defined $columns && defined $TMP_DIR;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($reg_file);
my $db_adaptor = $registry->get_DBAdaptor("homo_sapiens", "variation");
my $dbh = $db_adaptor->dbc;
my $dbname = $dbh->dbname;

# Get Variation_feature id
my $sql = qq{ SELECT MIN(${table_name}_id), MAX(${table_name}_id) FROM ${table_name} };
my $sth = $dbh->prepare( $sql );
$sth->execute();
my $vf = $sth->fetchall_arrayref();
my $min_id = $vf->[0]->[0];
my $max_id = $vf->[0]->[1];
my $TMP_FILE = "tmp.txt";
my $chunk = 1000000;

sub dumpPreparedSQL {
  my $dbVar = shift;
  my $tmp_num = shift;
  my $chunk = shift;
  my $id = shift;

  my $start = $chunk * $tmp_num;
  my $end = $chunk + $start;
  $end = $end < $id ? $end : $id;


  my $sql = qq{SELECT * FROM $table_name where ${table_name}_id > $start AND ${table_name}_id <= $end};

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

for my $tmp_num (map { $_ } $min_id/$chunk .. $max_id/$chunk) {
  dumpPreparedSQL($dbh, $tmp_num, $chunk, $max_id);
}

my @cols = split /,/, $columns;
my $cols_string = "";

foreach my $col (@cols){
    $cols_string .= "-k${col},${col}n "
}

system("sort --parallel=8 --buffer-size=30G ${cols_string}-o $TMP_DIR/sorted_$TMP_FILE $TMP_DIR/$TMP_FILE");
$dbh->do( "DROP TABLE IF EXISTS ${table_name}_tmp" ) or die "Failed to drop pre-existing temp table";
$dbh->do( "CREATE TABLE ${table_name}_tmp LIKE $table_name") or die "Failed to create temp table";
$dbh->do( "ALTER TABLE ${table_name}_tmp DISABLE KEYS") or die "Failed to disable keys in tmp table";
$dbh->do( "LOAD DATA LOCAL INFILE \"${TMP_DIR}/sorted_${TMP_FILE}\" INTO TABLE ${table_name}_tmp") or die "Failed to load infile in tmp table";
$dbh->do( "ALTER TABLE ${table_name}_tmp ENABLE KEYS") or die "Failed to re-enable keys in tmp table";

system("rm ${TMP_DIR}/*${TMP_FILE}*");

warn("TIP!!! REMEMBER TO EXECUTE SQL IN SERVER:\nDROP TABLE ${dbname}.${table_name};\nRENAME TABLE ${dbname}.${table_name}_tmp TO ${dbname}.${table_name};\n");

sub usage{

  die "\n\tUsage: sort_table_by_cols.pl
  \t\t-registry [registry file]
  \t\t-columns [comma separated numbers]
  \t\t-table_name [table name]
  \t\t-tmp [temp folder]\n\n";
}
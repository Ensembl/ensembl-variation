use strict;
use warnings;

package ImportUtils;

use Exporter;

our @ISA = ('Exporter');

our @EXPORT_OK = qw(dumpSQL debug create_and_load load);

our $TMP_DIR = "/ecs2/scratch5/ensembl/mcvicker/dbSNP";
our $TMP_FILE = 'tabledump.txt';

# successive dumping and loading of tables is typical for this process
# dump does effectively a select into outfile without server file system access
sub dumpSQL {
  my $db  = shift;
  my $sql = shift;

  local *FH;

  open FH, ">$TMP_DIR/$TMP_FILE";

  my $sth = $db->prepare( $sql, { mysql_use_result => 1 });
  $sth->execute();
  my $first;
  while ( my $aref = $sth->fetchrow_arrayref() ) {
    my @a = map {defined($_) ? $_ : '\N'} @$aref;
    print FH join("\t", @a), "\n";
  }

  close FH;

  $sth->finish();
}



# load imports a table, optionally not all columns
# if table doesnt exist, create a varchar(255) for each column
sub load {
  my $db = shift;
  my $tablename = shift;
  my @colnames = @_;


  my $cols = join( ",", @colnames );

  local *FH;
  open FH, "<$TMP_DIR/$TMP_FILE";
  my $sql;

  if ( @colnames ) {

    $sql = qq{
              LOAD DATA LOCAL INFILE '$TMP_DIR/$TMP_FILE'
              INTO TABLE $tablename( $cols )
             };
  } else {
    $sql = qq{
              LOAD DATA LOCAL INFILE '$TMP_DIR/$TMP_FILE'
              INTO TABLE $tablename
             };
  }

  $db->do( $sql );

  unlink( "$TMP_DIR/$TMP_FILE" );
}


#
# creates a table with specified columns and loads data that was dumped
# to a tmp file into the table.
#
# by default all columns are VARCHAR(255), but an 'i' may be added after the
# column name to make it an INT.  Additionally a '*' means add an index to
# the column.
#
# e.g.  create_and_load('mytable', 'col0', 'col1 *', 'col2 i', 'col3 i*');
#
sub create_and_load {
  my $db = shift;
  my $tablename = shift;
  my @cols = @_;

  my $sql = "CREATE TABLE $tablename ( ";

  my @col_defs;
  my @idx_defs;
  my @col_names;

  foreach my $col (@cols) {
    my ($name, $type) = split(/\s+/,$col);
    push @col_names, $name;

    if(defined($type) && $type =~ /i/) {
      push @col_defs, "$name INT";
    } else {
      push @col_defs, "$name VARCHAR(255)";
    }

    if(defined($type) && $type =~ /\*/) {
      push @idx_defs, "KEY ${name}_idx($name)";
    }
  }

  my $create_cols = join( ",\n", @col_defs, @idx_defs);


  $sql .= $create_cols.")";

  $db->do( $sql );

  load( $db, $tablename, @col_names );
}


sub debug {
  print STDERR @_, "\n";
}



1;

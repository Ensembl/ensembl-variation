use strict;
use warnings;

package ImportUtils;

use Exporter;

our @ISA = ('Exporter');

our @EXPORT_OK = qw(dumpSQL debug create_and_load load get_create_statement make_xml_compliant);

our $TMP_DIR = "/ecs4/scratch4/yuan/tmp";
our $TMP_FILE = 'tabledump.txt';

#ÊThis will strip non-xml-compliant characters from an infile, saving a backup in {infile name}.bak
# If no infile was specified, will use the tempfile but no backup will be kept. If no xml version was specified, will default to 1.1
# If the xml version was not recognized, will do nothing.
sub make_xml_compliant {
  my $infile = shift;
  my $version = shift;
  
  my $keep_backup = defined($infile);
  
  $infile ||= $TMP_DIR . "/" . $TMP_FILE;
  $version ||= "1.1";
  
  my @ARGV_bak = @ARGV;
  @ARGV = ($infile);
  $^I = '.bak';
  while (<>) {
    if ($version =~ m/1\.1/) {
      s/[^\x01-\x{D7FF}\x{E000}-\x{FFFD}\x{10000}-\x{10FFFF}]//go;
      s/[\x01-\x08\x0B-\x0C\x0E-\x1F\x7F-\x84\x86-\x9F]//go;
    }
    elsif ($version =~ m/1\.0/) {
      s/[^\x09\x0A\x0D\x20-\x{D7FF}\x{E000}-\x{FFFD}\x{10000}-\x{10FFFF}]//go;
    }
    print;
  }
  # Remove the backup if an input file wasn't supplied
  unlink($infile . ".bak") if (!$keep_backup);
  # Restore the @ARGV variable
  @ARGV = @ARGV_bak;
}

# successive dumping and loading of tables is typical for this process
# dump does effectively a select into outfile without server file system access
sub dumpSQL {
  my $db  = shift;
  my $sql = shift;
  my $sth = $db->prepare( $sql );
  dumpPreparedSQL($sth);
  $sth->finish();
}

sub dumpPreparedSQL {
  my $sth = shift;
  
  local *FH;
  my $counter = 0;
  open( FH, ">$TMP_DIR/$TMP_FILE" )
      or die( "Cannot open $TMP_DIR/$TMP_FILE: $!" );
      
  $sth->{mysql_use_result} = 1;
  $sth->execute();
  my $first;
  while ( my $aref = $sth->fetchrow_arrayref() ) {
    my @a = map {defined($_) ? $_ : '\N'} @$aref;
    print FH join("\t", @a), "\n";
  }

  close FH;
}

# load imports a table, optionally not all columns
# if table doesnt exist, create a varchar(255) for each column
sub load {
  my $db = shift;
  my $tablename = shift;
  my @colnames = @_;

  my $cols = join( ",", @colnames );

  my $table_file = "$TMP_DIR/$tablename\_$$\.txt";
  rename("$TMP_DIR/$TMP_FILE", $table_file);
   
#  my $host = $db->host();
#  my $user = $db->user();
#  my $pass = $db->pass();
#  my $port = $db->port();
#  my $dbname = $db->dbname();

#  my $call = "mysqlimport -c $cols -h $host -u $user " .
#    "-p$pass -P$port $dbname $table_file";

#  system($call);

#  unlink("$TMP_DIR/$tablename.txt");


##### Alternative way of doing same thing
  my $sql;

  #need to find out if possible use the LOCAL option

#  my $local_option = 'LOCAL'; #by default, use the LOCAL option
#  if( `hostname` =~ /^bc/ ){ # No LOCAL on bcs nodes
#    $local_option = '';
#  }
#  elsif( ! -e $table_file ){ # File is not on local filesystem
#    $local_option = '';
#  }

  my $local_option = 'LOCAL';

 # my $host = `hostname`;
 # chop $host;
 #  $host =~ /(ecs\d+)/; #get the machine, only use LOCAL in ecs machines (ecs2, ecs4)
 #  my $local_option = '';
 #  #the script is running in ecs machine, let's find out if the file is in the same machine, too
 #  if ($1){
 #      if ($table_file =~ /$1/){
 #	  $local_option = 'LOCAL';
 #      }
 #  }

  
   if ( @colnames ) {

     $sql = qq{
               LOAD DATA $local_option INFILE "$table_file"
               INTO TABLE $tablename( $cols )
              };
   } else {
     $sql = qq{
               LOAD DATA $local_option INFILE "$table_file"
               INTO TABLE $tablename
              };
   }
   $db->do( $sql );

   unlink( "$table_file" );
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
    }
    elsif (defined($type) && $type =~ /f/) {
      push @col_defs, "$name FLOAT";
    }
    elsif (defined($type) && $type =~ /l/) {
      push @col_defs, "$name TEXT";
    } 
    else {
      push @col_defs, "$name VARCHAR(255)";
    }

    if(defined($type) && $type =~ /\*/) {
      push @idx_defs, "KEY ${name}_idx($name)";
    }
  }

  my $create_cols = join( ",\n", @col_defs, @idx_defs);


  $sql .= $create_cols.")";

  $sql .= " MAX_ROWS = 100000000" if ($tablename =~ /^tmp.*gty$/); #need to make bigger this table for human
  $db->do( $sql );

  load( $db, $tablename, @col_names );
}

#
#ÊGets the create statement to create the desired table from the master_schema_variation database
#
sub get_create_statement {
  my $dbc = shift;
  my $table = shift;
  my $db_name = shift;

  if (defined($db_name)) {
    $table = $db_name . '.' . $table;
  }

  my $stmt = qq{
    SHOW CREATE TABLE
      $table
  };
  my $result = $dbc->db_handle->selectall_arrayref($stmt)->[0][1];
  
  return $result;
}

sub debug {
  print STDERR @_, "\n";
}



1;

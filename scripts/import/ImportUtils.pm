=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

use strict;


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use warnings;

package ImportUtils;

use Exporter;

our @ISA = ('Exporter');

our @EXPORT_OK = qw(dumpSQL debug create create_and_load load loadfile get_create_statement make_xml_compliant update_table);

our $TMP_DIR = "/tmp";
our $TMP_FILE = 'tabledump.txt';

# This will strip non-xml-compliant characters from an infile, saving a backup in {infile name}.bak
# If no infile was specified, will use the tempfile but no backup will be kept. If no xml version was specified, will default to 1.1
# A replacement character or string can be passed
# If the xml version was not recognized, will do nothing.
sub make_xml_compliant {
  my $infile = shift;
  my $version = shift;
  my $replacement = shift;
  
  my $keep_backup = defined($infile);
  
  $infile ||= $TMP_DIR . "/" . $TMP_FILE;
  $version ||= "1.1";
  $replacement ||= "";
  
  my @ARGV_bak = @ARGV;
  @ARGV = ($infile);
  $^I = ".bak";
  while (<>) {
    if ($version =~ m/1\.1/) {
      s/[^\x01-\x{D7FF}\x{E000}-\x{FFFD}\x{10000}-\x{10FFFF}]/$replacement/go;
      s/[\x01-\x08\x0B-\x0C\x0E-\x1F\x7F-\x84\x86-\x9F]/$replacement/go;
    }
    elsif ($version =~ m/1\.0/) {
      s/[^\x09\x0A\x0D\x20-\x{D7FF}\x{E000}-\x{FFFD}\x{10000}-\x{10FFFF}]/$replacement/go;
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
  my $dbe = shift;  ## handle postgreSQL differently

  if($dbe =~/pg|postgreSQL/i ){
    dumpSQL_PG($db, $sql);
  }
  else{
   my $sth = $db->prepare( $sql );
   dumpPreparedSQL($sth);
   $sth->finish();
  }
 
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

## use postgresql cursors to avoid memory issues on large exports
sub dumpSQL_PG {
 my $dbh = shift; 
 my $sql = shift;

  local *FH;
  my $counter = 0;
  open( FH, ">$TMP_DIR/$TMP_FILE" )
      or die( "Cannot open $TMP_DIR/$TMP_FILE: $!" );

   $dbh->db_handle->begin_work();
   $dbh->do("DECLARE csr CURSOR  FOR $sql");
   while (1) {
     my $sth = $dbh->prepare("fetch 100000 from csr");
     $sth->execute;
     last if 0 == $sth->rows;
   
     while ( my $aref = $sth->fetchrow_arrayref() ) {
       my @a = map {defined($_) ? $_ : '\N'} @$aref;
       print FH join("\t", @a), "\n";
     }
  }
   $dbh->do("CLOSE csr");
   $dbh->db_handle->rollback();
   close FH;
}



# load imports a table, optionally not all columns
# if table doesnt exist, create a varchar(255) for each column
sub load {
  loadfile("$TMP_DIR/$TMP_FILE",@_);
}

sub loadfile {
  my $loadfile = shift;
  my $db = shift;
  my $tablename = shift;
  my @colnames = @_;
  
  my $cols = join( ",", @colnames );

  my $table_file = "$TMP_DIR/$tablename\_$$\.txt";
  rename($loadfile, $table_file);
   
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
# creates a table with specified columns
#
# by default all columns are VARCHAR(255), but an 'i' may be added after the
# column name to make it an INT.  Additionally a '*' means add an index to
# the column.
#
# e.g.  create($db, 'mytable', 'col0', 'col1 *', 'col2 i', 'col3 i*');
#

sub create {
  my $db = shift;
  my $tablename = shift;
  my @cols = @_;

  my $sql = "CREATE TABLE $tablename ( ";

  my @col_defs;
  my @idx_defs;
  my @col_names;

  foreach my $col (@cols) {
    my ($name, $type, $nullable, $unsigned) = split(/\s+/,$col);
    push @col_names, $name;

    my $null ="";
    if (defined($nullable) && $nullable =~/not_null/){$null =" NOT NULL";}

    if(defined($type) && $type =~ /i/ && defined $unsigned ) {
      push @col_defs, "$name INT unsigned $null";
    }
    elsif(defined($type) && $type =~ /i/) {
      push @col_defs, "$name INT  $null";
    }
    elsif (defined($type) && $type =~ /f/) {
      push @col_defs, "$name FLOAT $null";
    }
    elsif (defined($type) && $type =~ /l/) {
      push @col_defs, "$name TEXT $null";
    } 
    elsif (defined($type) && $type =~ /d/) {
       push @col_defs, "$name DOUBLE $null";
    }
    elsif (defined($type) && $type =~ /v/) {
      $type =~ s/v//;
      my $len = 255;
       push @col_defs, "$name VARCHAR($len) $null";
    }
    else {
      push @col_defs, "$name VARCHAR(255) $null";
    }

    if(defined($type) && $type =~ /\*/) {
      push @idx_defs, "KEY ${name}_idx($name)";
    }
  }

  my $create_cols = join( ",\n", @col_defs, @idx_defs);


  $sql .= $create_cols.")";

  $sql .= " MAX_ROWS = 100000000" if ($tablename =~ /^tmp.*gty$/); #need to make bigger this table for human

  $sql .= " ENGINE = 'MyISAM' "; ##may not be default engine

  $db->do( $sql );

  return @col_names;
}

#
# creates a table with specified columns and loads data that was dumped
# to a tmp file into the table.
#
# by default all columns are VARCHAR(255), but an 'i' may be added after the
# column name to make it an INT.  Additionally a '*' means add an index to
# the column.
#
# e.g.  create_and_load($db, 'mytable', 'col0', 'col1 *', 'col2 i', 'col3 i*');
#

sub create_and_load {
  my $db = shift;
  my $tablename = shift;
  my @cols = @_;

  my @col_names = create($db,$tablename,@cols);

  load( $db, $tablename, @col_names );
}

#
#ÊGets the create statement to create the desired table from the master_schema_variation database
#
sub get_create_statement_from_db {
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

#ÊGet the create statement for a table from e.g. the table.sql schema definition file
sub get_create_statement {
  my $table = shift;
  my $sql_file = shift;
  
  # Parse the file
  open(FH,'<',$sql_file) or die("Could not open $sql_file for reading");
  
  # Parse the file into a string
  my $contents = "";
  while (<FH>) {
    chomp;
    $contents .= "$_ ";
  }
  close(FH);
  
  # Grab the correct create statement
  my ($stmt) = $contents =~ m/(create table $table [^\;]+)\;/i;
  
  if (!$stmt) {
    warn("Could not find CREATE TABLE statement for $table");
  }
  
  return $stmt;
}

sub update_table {
  my $db = shift;
  my $source_table = shift;
  my $target_table = shift;
  my $source_key = shift;
  my $target_key = shift;
  my $source_col = shift;
  my $target_col = shift;
  my $clean = shift;
  
  die("Incorrect arguments supplied to update_table\n") unless
    defined($db) &&
    defined($source_table) &&
    defined($target_table) &&
    defined($source_key) &&
    defined($target_key) &&
    defined($source_col) &&
    defined($target_col);
  
  my ($stmt, $sth);
  
  # get columns of source table
  $stmt = qq{
    DESCRIBE $source_table
  };
  $sth = $db->prepare($stmt);
  $sth->execute;
  
  my @source_cols = map {$_->[0]} @{$sth->fetchall_arrayref};
  $sth->finish;
  
  die("No columns found in table $source_table\n") unless @source_cols;
  die("Key column $source_key not found in source table $source_table\n") unless grep {$_ eq $source_key} @source_cols;
  die("Data column $source_col not found in source table $source_table\n") unless grep {$_ eq $source_col} @source_cols;
  
  # get columns of target table
  $stmt = qq{
    DESCRIBE $target_table
  };
  $sth = $db->prepare($stmt);
  $sth->execute;
  
  my @target_cols = map {$_->[0]} @{$sth->fetchall_arrayref};
  $sth->finish;
  
  die("No columns found in table $target_table\n") unless @target_cols;
  die("Key column $target_key not found in target table $target_table\n") unless grep {$_ eq $target_key} @target_cols;
  die("Data column $target_col not found in target table $target_table\n") unless grep {$_ eq $target_col} @target_cols;
  
  # construct columns to select
  my (@select_cols, $select_cols);
  
  foreach my $col(@target_cols) {
    if($col eq $target_col) {
      if($clean) {
        push @select_cols, $source_table.'.'.$source_col;
      }
      else {
        push @select_cols, qq{if($source_table.$source_col is null, $target_table.$col, $source_table.$source_col)};
      }
    }
    else {
      push @select_cols, $target_table.'.'.$col;
    }
  }
  
  $select_cols = join ",", @select_cols;
  
  # create tmp table
  my $tmp_table = $target_table.'_tmp_'.$$;
  
  $stmt = qq{
    CREATE TABLE $tmp_table LIKE $target_table
  };
  $db->do($stmt);
  
  # construct SQL
  $stmt = qq{
    INSERT INTO $tmp_table
    SELECT $select_cols
    FROM $target_table
    LEFT JOIN $source_table
    ON $target_table.$target_key = $source_table.$source_key
  };
  
  $db->do($stmt);
  
  # check things went OK
  $sth = $db->prepare("SELECT COUNT(*) FROM $target_table");
  $sth->execute();
  my $old_count = $sth->fetchall_arrayref()->[0]->[0];
  $sth->finish();
  
  $sth = $db->prepare("SELECT COUNT(*) FROM $tmp_table");
  $sth->execute();
  my $new_count = $sth->fetchall_arrayref()->[0]->[0];
  $sth->finish();
  
  die("New tmp table ($new_count) does not have the same number of rows as original target table ($old_count)\n") unless $old_count == $new_count;
  
  # rename and drop
  $db->do(qq{DROP TABLE $target_table});
  $db->do(qq{RENAME TABLE $tmp_table TO $target_table});
}

sub debug {
  print STDERR @_, "\n";
}



1;

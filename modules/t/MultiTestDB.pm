 

=pod

=head1 NAME - EnsTestDB

=head1 SYNOPSIS

=head1 DESCRIPTION


=head1 METHODS

=cut

package MultiTestDB;

use vars qw(%ENV);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use strict;

use DBI;
use Data::Dumper;


#homo sapiens is used if no species is specified
my $DEFAULT_SPECIES  = 'homo_sapiens';

#configuration file extension appended onto species name
my $FROZEN_CONF_EXT  = '.MultiTestDB.frozen.conf';

my $CONF_FILE    = 'MultiTestDB.conf';

my $DUMP_DIR = 'test-genome-DBs';

sub new {
  my( $pkg, $species ) = @_;

  my $self = bless {}, $pkg;

  # go and grab the current directory and store it away
  my $curr_dir = $ENV{'PWD'} . "/".__FILE__;
  $curr_dir =~ s/MultiTestDB.pm$//;

  $self->curr_dir($curr_dir);

  unless($species) {
    $species = $DEFAULT_SPECIES;
  }

  $self->species($species);


  if ( -e $self->curr_dir . $species . $FROZEN_CONF_EXT) {
    $self->load_config;
  }
  else {
    #load the databases and generate the conf hash
    $self->load_databases;

    #freeze configuration in a file
    $self->store_config;
  }


  #generate the db_adaptors from the $self->{'conf'} hash
  $self->create_adaptors;

  return $self;
} 


#
# load config into $self->{'conf'} hash
#
sub load_config {
  my $self = shift;

  my $conf = $self->curr_dir . $self->species . $FROZEN_CONF_EXT;
    
  eval {
    $self->{'conf'} = do $conf; #reads file into $self->{'conf'}
  };
    
  if($@) {
    die("Could not read frozen configuration file '$conf'\n");
  }
}



#
# Store $self->{'config'} hash into a file
#
sub store_config {
  my $self = shift;

  my $conf = $self->curr_dir . $self->species . $FROZEN_CONF_EXT;

  local *FILE;

  open(FILE, ">$conf") or die "Could not open config file ".$conf."\n";

  my $string = Dumper($self->{'conf'});

  #strip off leading '$VAR1 = '
  $string =~ s/^[\$]VAR1\s*=//;

  #store config in file
  print FILE $string;
  
  close FILE;
}



#create a set of adaptors based on the $self->{'conf'} hash
sub create_adaptors {
  my $self = shift;

  #establish a connection to each of the databases in the configuration
  foreach my $dbtype (keys %{$self->{'conf'}}) {

    my $db = $self->{'conf'}->{$dbtype};
    my $adaptor;
    my $module = $db->{'module'};

    #try to instantiate an adaptor for this database 
    eval {

      # require needs /'s rather than colons
      if ( $module =~ /::/ ) {
	$module =~ s/::/\//g;
      }
      require "${module}.pm";

      # but switch back for the new instantiation
      $module =~ s/\//::/g;

      $adaptor = $module->new(-dbname => $db->{'dbname'},
			      -user   => $db->{'user'},
			      -pass   => $db->{'pass'},
			      -port   => $db->{'port'},
			      -host   => $db->{'host'},
			      -driver => $db->{'driver'});
    };
	
    if ($@) {
      warning("WARNING: Could not instantiate $dbtype DBAdaptor:\n$@");
    } else {
      $self->{'db_adaptors'}->{$dbtype} = $adaptor;
    }
  }
}





sub load_databases {
  my ($self) = shift;

  print STDERR "\nTrying to load [$self->{'species'}] databases\n";

  #create database from conf and from zip files 
  my $db_conf = do $self->curr_dir . $CONF_FILE;

  my $port   = $db_conf->{'port'};
  my $driver = $db_conf->{'driver'};
  my $host   = $db_conf->{'host'};
  my $pass   = $db_conf->{'pass'};
  my $user   = $db_conf->{'user'};
  my $zip    = $db_conf->{'zip'};

  #create a config hash which will be frozen to a file
  $self->{'conf'} = {};

  #connect to the database
  my $locator = "DBI:".$driver.":host=".$host.";port=".$port;
  my $db = DBI->connect($locator, $user, $pass, {RaiseError => 1});

  unless($db) {
    warning("Can't connect to database $locator");
    return;
  }

  #only unzip if there are non-preloaded datbases
 UNZIP: foreach my $dbtype (keys %{$db_conf->{'databases'}}) {
    if (( ! exists $db_conf->{'preloaded'}->{$dbtype} ) ||
        ( ! _db_exists( $db, $db_conf->{'preloaded'}{$dbtype}) )) {
      #unzip database files
      $self->unzip_test_dbs($self->curr_dir . $zip);
      last UNZIP;
    }
  }


  #create a database for each database specified
  foreach my $dbtype (keys %{$db_conf->{'databases'}}) {
    #don't create a database if there is a preloaded one specified
    if (( $db_conf->{'preloaded'}->{$dbtype} ) &&
        ( _db_exists( $db,$db_conf->{'preloaded'}->{$dbtype} ))) {
      #copy the general config into a dbtype specific config 
      $self->{'conf'}->{$dbtype} = {};
      %{$self->{'conf'}->{$dbtype}} = %$db_conf;
      $self->{'conf'}->{$dbtype}->{'module'} = $db_conf->{'databases'}->{$dbtype};

      # it's not necessary to store the databases and zip bits of info
      delete $self->{'conf'}->{$dbtype}->{'databases'};
      delete $self->{'conf'}->{$dbtype}->{'zip'};


      #store the temporary database name in the dbtype specific config
      $self->{'conf'}->{$dbtype}->{'dbname'} = $db_conf->{'preloaded'}->{$dbtype};
      $self->{'conf'}->{$dbtype}->{'preloaded'} = 1;
    } else {

      $self->{'conf'}->{$dbtype} = {};
      %{$self->{'conf'}->{$dbtype}} = %$db_conf;
      $self->{'conf'}->{$dbtype}->{'module'} = $db_conf->{'databases'}->{$dbtype};

      # it's not necessary to store the databases and zip bits of info
      delete $self->{'conf'}->{$dbtype}->{'databases'};
      delete $self->{'conf'}->{$dbtype}->{'zip'};

      #create a unique random dbname    
      my $dbname = $db_conf->{'preloaded'}->{$dbtype};
      if ( ! defined $dbname ) {
        $dbname = $self->_create_db_name($dbtype);
        delete $self->{'conf'}->{$dbtype}->{'preloaded'};
      } else {
        $self->{'conf'}->{$dbtype}->{'preloaded'} = 1;
      }

      #store the temporary database name in the dbtype specific config
      $self->{'conf'}->{$dbtype}->{'dbname'} = $dbname;

      print STDERR "\nCreating db $dbname";
      
      unless($db->do("CREATE DATABASE $dbname")) {
        warning("Could not create database [$dbname]");
        return;
      }

      #copy the general config into a dbtype specific config 

      $db->do("use $dbname");
    
      #load the database with data
      my $dir = $self->curr_dir . "$DUMP_DIR/".$self->species."/$dbtype";
      local *DIR;

      unless(opendir(DIR, $dir)) {
        warning("could not open dump directory '$dir'");
        return;
      }

      my @files = readdir DIR;

      local *FILE;

      #read in table creat statements from *.sql files and process them with DBI

      foreach my $sql_file (grep /\.sql$/, @files) {
	
        $sql_file = "$dir/$sql_file";
	
        unless(-f $sql_file && -r $sql_file) {
          warning("could not read SQL file '$sql_file'\n");
          next;
        }

        open(FILE, $sql_file);
	
        my $sql_com ='';
	
        while (<FILE>) {
          next if ( /^#/ );  # ignore comments
          next unless ( /\S/ ); # ignore lines of white spaces

          $sql_com .= $_;
        }
        $sql_com =~ s/;$//; # chop off the last ;

        $db->do($sql_com);

        close FILE;

        #import data from the txt files of the same name
        $sql_file  =~ /.*\/(.*)\.sql/;
        my $tablename = $1;

        (my $txt_file = $sql_file) =~ s/\.sql$/\.txt/;

        unless(-f $txt_file && -r $txt_file) {
          warning("could not read data file '$txt_file'\n");
          next;
        }

        $db->do( "load data local infile '$txt_file' into table $tablename" );

      }
    }
    print STDERR "\n";
    closedir DIR;

  }

  $db->disconnect;

}


sub unzip_test_dbs {
  my ($self, $zipfile) = @_;

  if (-e $self->curr_dir . $DUMP_DIR) {
    warning("Test genome dbs already unpacked\n");
    return;
  }

  unless($zipfile) {
    throw("zipfile argument is required\n");
  }

  unless(-f $zipfile) {
    warning("zipfile could not be found\n");
    return;
  }

  # unzip the zip file quietly

  system ( "unzip -q $zipfile -d ". $self->curr_dir );
}




sub get_DBAdaptor {
  my ($self, $type) = @_;

  unless($type) {
    die('type arg must be specified\n');
  }

  unless($self->{'db_adaptors'}->{$type}) {
    warning("dbadaptor of type $type is not available\n");
    return undef;
  }

  return $self->{'db_adaptors'}->{$type};
}



=head2 hide

  Arg [1]    : string $dbtype
               The type of the database containing the hidden table
  Arg [2]    : string $table
               The name of the table to hide
  Example    : $multi_test_db->hide('core', 'gene', 'transcript', 'exon');
  Description: Hides the contents of specific table(s) in the specified db.
               The table(s) are first renamed and an empty table are created 
               in their place by reading the table schema file.
  Returntype : none
  Exceptions : thrown if the adaptor for dbtype is not available
               thrown if both arguments are not defined
               warning if a table is already hidden
               warning if a table cannot be hidden because its schema file 
               cannot be read
  Caller     : general

=cut

sub hide {
  my ($self, $dbtype, @tables) = @_;

  unless($dbtype && @tables) {
    die("dbtype and table args must be defined\n");
  }

  my $adaptor = $self->get_DBAdaptor($dbtype);

  unless($adaptor) {
    die "adaptor for $dbtype is not available\n";
  }

  foreach my $table (@tables) {

    if($self->{'conf'}->{$dbtype}->{'hidden'}->{$table}) {
      warning("table '$table' is already hidden and cannot be hidden again\n");
      next;
    }

    my $hidden_name = "_hidden_$table";

    #copy contents of table into a temp table

    my $sth =
      $adaptor->prepare("CREATE TABLE $hidden_name " .
                        "SELECT * FROM $table");

    $sth->execute();
    $sth->finish();

    #delete the contents of the original table
    $sth = $adaptor->prepare("DELETE FROM $table");
    $sth->execute();
    $sth->finish();

    #update the hidden table config
    $self->{'conf'}->{$dbtype}->{'hidden'}->{$table} = $hidden_name;
  }
}



=head2 restore

  Arg [1]    : (optional) $dbtype 
               The dbtype of the table(s) to be restored. If not specified all
               hidden tables in all the databases are restored.
  Arg [2]    : (optional) @tables
               The name(s) of the table to be restored.  If not specified all
               hidden tables in the database $dbtype are restored.
  Example    : $self->restore('core', 'gene', 'transcript', 'exon');
  Description: Restores a list of hidden tables. The current version of the
               table is discarded and the hidden table is renamed.
  Returntype : none
  Exceptions : thrown if the adaptor for a dbtype cannot be obtained
  Caller     : general

=cut

sub restore {
  my ($self, $dbtype, @tables) = @_;

  if(!$dbtype) {
    #restore all of the tables in every dbtype

    foreach my $dbtype (keys %{$self->{'conf'}}) {
      $self->restore($dbtype);
    }

    #lose the hidden table details
    delete $self->{'conf'}->{'hidden'};

    return;
  }

  my $adaptor = $self->get_DBAdaptor($dbtype);
  unless($adaptor) {
    die "Adaptor for $dbtype is not available";
  }

  if(!@tables) {
    #restore all of the tables for this db
    @tables = keys %{$self->{'conf'}->{$dbtype}->{'hidden'}};
  }

  foreach my $table (@tables) {
    my $hidden_name = $self->{'conf'}->{$dbtype}->{'hidden'}->{$table};
	
    #delete current contents of table
    my $sth = $adaptor->prepare("delete from $table");
    $sth->execute();
    $sth->finish();

    #copy contents of tmp table back into main table
    $sth = $adaptor->prepare("insert into $table " .
                             "select * from $hidden_name");
    $sth->execute();
    $sth->finish();

    #drop temp table
    $sth = $adaptor->prepare("drop table $hidden_name");
    $sth->execute();
    $sth->finish();

    #delete value from hidden table config
    delete $self->{'conf'}->{$dbtype}->{'hidden'}->{$table};
  }
}

=head2 save

  Arg [1]    : string $dbtype
               The type of the database containing the hidden/saved table
  Arg [2]    : string $table
               The name of the table to save
  Example    : $multi_test_db->save('core', 'gene', 'transcript', 'exon');
  Description: Saves the contents of specific table(s) in the specified db.
               The table(s) are first renamed and an empty table are created 
               in their place by reading the table schema file.  The contents
               of the renamed table(s) are then copied back into the newly
               created tables.  The method piggy-backs on the hide method
               and simply adds in the copying/insertion call.
  Returntype : none
  Exceptions : thrown if the adaptor for dbtype is not available
               warning if a table cannot be copied if the hidden table does not 
               exist
  Caller     : general

=cut

sub save {
  my ($self, $dbtype, @tables) = @_;

  # use the hide method to build the basic tables
  $self->hide($dbtype, @tables);

  my $adaptor = $self->get_DBAdaptor($dbtype);

  unless($adaptor) {
    die "adaptor for $dbtype is not available\n";
  }

  my $hidden_name = '';
  foreach my $table (@tables) {

    # only do if the hidden table exists
    if($self->{'conf'}->{$dbtype}->{'hidden'}->{$table}) {

      $hidden_name = "_hidden_$table";

      #copy the data from the hidden table into the new table
      my $sth = $adaptor->prepare("insert into $table " .
                                  "select * from $hidden_name"); 
      $sth->execute;
    }
    else {
      warning("hidden table '$hidden_name' does not exist so saving is not possible\n");
    }
  }
}


sub _db_exists {
  my ( $db, $db_name ) = @_;

  my $db_names = $db->selectall_arrayref( "show databases" );
  for my $db_name_ref ( @$db_names ) {
    if( $db_name_ref->[0] eq $db_name ) {
      return 1;
    }
  }

  return 0;
}



sub compare {
  my ($self, $dbtype, $table) = @_;

  warning("save method not yet implemented\n");

}


sub species {
  my ($self, $species) = @_;

  if($species) {
    $self->{'species'} = $species;
  }

  return $self->{'species'};
}



sub curr_dir {
  my ($self, $cdir) = @_;

  if($cdir) {
    $self->{'_curr_dir'} = $cdir;
  }

  return $self->{'_curr_dir'};
}



sub _create_db_name {
    my( $self, $dbtype ) = @_;

    my @t_info = localtime;

    my $date = join ( "_", $t_info[3],$t_info[4]+1);  
    my $time = join ( "", $t_info[2],$t_info[1],$t_info[0]);  

    my $species = $self->species;

    # create a unique name using host and date / time info
    my $db_name = "_test_db_${species}_${dbtype}_".$ENV{'USER'}."_".$date."_".$time;

    return $db_name;
}





sub cleanup {
  my $self = shift;

  #delete the unpacked schema and data files
  $self->_delete_files($self->curr_dir . $DUMP_DIR);

  #remove all of the handles on dbadaptors
  foreach my $dbtype (keys %{$self->{'db_adaptors'}}) {
    delete $self->{'db_adaptors'}->{$dbtype};
  }

  #delete each of the created temporary databases
  foreach my $dbtype (keys %{$self->{'conf'}}) {

    my $db_conf = $self->{'conf'}->{$dbtype};
    my $host   = $db_conf->{'host'};
    my $user   = $db_conf->{'user'};
    my $pass   = $db_conf->{'pass'};
    my $port   = $db_conf->{'port'};
    my $driver = $db_conf->{'driver'};
    my $dbname = $db_conf->{'dbname'};
    
    if( $db_conf->{'preloaded'} ) {
      next;
    }
    #connect to the database
    my $locator = "DBI:".$driver.":host=".$host.";port=".$port;

    my $db = DBI->connect($locator, $user, $pass, {RaiseError => 1});
      
    unless($db) {
      die "Can't connect to database $locator";
    }
    
    print STDERR "Dropping db $dbname \n";
    $db->do("DROP database $dbname");
    $db->disconnect;
  }

  my $conf_file = $self->curr_dir . $self->species . $FROZEN_CONF_EXT;
  
  #delete the frozen config file
  if(-e $conf_file && -f $conf_file) {
    unlink $conf_file;
  }
}


sub _delete_files {
  my ($self, $dir) = @_;

  local *DIR;
  opendir DIR, $dir;

  #ignore files starting with '.'

  my @files = grep !/^\./, readdir DIR;

  foreach my $file (@files) {

    $file = $dir ."/". $file;
    if(-d $file) {

      #call recursively on subdirectories
      $self->_delete_files($file);

    } else {
      unlink $file;
    }
  }
  closedir DIR;

  rmdir $dir;
}


sub DESTROY {
    my( $self ) = shift;

    if($ENV{'RUNTESTS_HARNESS'}) {
      #restore tables, do nothing else we want to use the database for 
      #the other tests as well
      $self->restore;
    } else {
      #we are runnning a stand-alone test, cleanup created databases
      print STDERR "\nCleaning up....\n";
      #restore database state since we may not actually delete it in the
      #cleanup - it may be defined as a preloaded db
      $self->restore;
      $self->cleanup;
    }
}




1;


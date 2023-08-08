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

debug($config, "Importing variation set starting now");

my $reg_file = $config->{registry};
$registry->load_all($reg_file);
my $db_adaptor = $registry->get_DBAdaptor("homo_sapiens", "variation");
my $dbh = $db_adaptor->dbc;
my $dbname;
$dbname = $dbh->dbname;
debug($config, "Connected to the new database $dbname");

my $TMP_DIR = $config->{tmp};

my $new_vf_file = "variation_feature.txt";
my $tmp_vset = "variation_name_set.txt";
my $tmp_merged = "new_variation_set_variation.txt";
my $tmp_vs_file = "vset_concat.txt";
my $old_dbh;
my $old_reg_file = $config->{old_registry} if (defined ($config->{old_registry}));
my $release =  $config->{release} if (defined ($config->{release}));

 my $old_dbname; 

if ($old_reg_file) {
  $registry->load_all($old_reg_file);
  $db_adaptor = $registry->get_DBAdaptor("homo_sapiens", "variation");
  # Connect to DBI
  $old_dbh = $db_adaptor->dbc;
  $old_dbname = $old_dbh->dbname;
  debug($config, "Connected to the old database $old_dbname");
} else {
  my $old_release = $release - 1;
  $old_dbname = $dbh->dbname =~ s/_${release}_/_${old_release}_/gr;
  my $old_host = $dbh->host;
  my $old_port = $dbh->port;
  #Connect to DBI
  $old_dbh = DBI->connect("DBI:mysql:database=${old_dbname};host=${old_host};port=${old_port}", "ensro", "");
  debug($config, "Connected to the old database $old_dbname");
}

debug($config, "Selecting minimum and maximum variation feature");
my $sql = qq{ SELECT MIN(variation_feature_id), MAX(variation_feature_id) FROM variation_feature};
my $sth = $dbh->prepare($sql);
$sth->execute();
my $vf = $sth->fetchall_arrayref();
my $min_id = $vf->[0]->[0];
my $chunk = 1000000;
my $max_id = $vf->[0]->[1];


debug($config, "Dumping the variation sets and the new variation feature into files");
for my $tmp_num (map { $_ } $min_id/$chunk .. $max_id/$chunk) {
  dump_old_sql_variation_sets($old_dbh, $tmp_num, $chunk, $max_id);
}

for my $tmp_num (map { $_ } $min_id/$chunk .. $max_id/$chunk) {
  dump_new_variation_feature($dbh, $tmp_num, $chunk, $max_id);
}

debug($config, "Sorting the files");
system(sort -u $new_vf_file);
system(sort -u $tmp_vset);

debug($config, "Merging the files based on the new variation_id and the old variation set id");
create_merged_file($tmp_vset, $new_vf_file, $tmp_merged);

debug($config, "Altering variation set variation table");
temp_table($dbh);

debug($config, "Loading new variation sets from merged file");
load_all_variation_sets($dbh, $tmp_merged);


debug($config, "Recalculating the variation sets"); # takes from the merged file and recalculates and creates the concatenate file that will be used to update variation feature
recalculate($tmp_merged, $tmp_vs_file);

debug($config, "Updating the variation feature table");
update_variation_feature_table($dbh, $tmp_vs_file);


debug($config, "Adding failed variation to variation set");
$dbh->do( qq{ 
  INSERT IGNORE INTO variation_set_variation (variation_id, variation_set_id)
  SELECT DISTINCT variation_id, 1 
  FROM failed_variation; 
}) or die "Failed to add failed to variation_set_variation table";

$dbh->do(qq{
  ALTER TABLE variation_set_variation ENABLE keys;
}) or die "Failed to alter variation_set_variation keys";

 

sub temp_table { 
  my $dbhvar = shift; 
   
  # creating a temp table and altering the table by disabling the keys 
  my $alter_sql = $dbhvar->prepare(qq{ALTER TABLE variation_set_variation DISABLE keys});
  #executing the sql 
  $alter_sql->execute();

  $alter_sql->finish();
}

sub create_merged_file {
  # this would merge the files variation_feature.txt and variation_name_set.txt using the column they share in common which is the var_name and creates a new file 
  # which is variation_name_set.txt containing the var_id and var_set_id which would be loaded to the new variation_set table
  my $tmp_vset = shift; 
  my $new_vf_file = shift; 
  my $tmp_merged = shift;

  open FILE1, "<", "$tmp_vset" or die "Cannot open $tmp_vset: $!";
  open FILE2, "<", "$new_vf_file" or die "Cannot open $new_vf_file: $!";
  open OUTPUT, ">", "$tmp_merged" or die "Cannot open $tmp_merged: $!";
  my %data;

  while (<FILE1>) {
    chomp;
    my $key = (split)[0]; # to access the variation name
    my @value = (split)[1]; # to access set id
    # checking if the key exists in the dictionary and if it exists it appends the array reference
    if (exists $data{$key}) {
       push @{$data{$key}}, \@value; # push the array reference to the hash
    } else { # if it does  not exist, it just creates a new key and an array
      $data{$key} = [\@value];
    }
  }
  while (<FILE2>) {
    chomp;
    my $var_name = (split)[1];
    my $var_id = (split)[0];
    if(exists $data{$var_name}) { # checking if the variation_name exists in the dictionary if it does 
      my $array_ref = $data{$var_name}; # creating an array reference if the var_name exists in the created hash 
      foreach my $set_id (@$array_ref) {
        print OUTPUT "$var_id\t@$set_id[0]\n"; # to access the values in the array reference and the variation_id to create a new file
      }
    }
  }
  close FILE1;
  close FILE2;
  close OUTPUT;

}

sub load_all_variation_sets {
  # this would load the created file from the function create_merged_file to a temp table that was created using temp_table function
  my $dbhvar = shift;
  my $load_file = shift;

  my $sql = qq{INSERT INTO variation_set_variation (variation_id, variation_set_id ) VALUES (?, ?)};
  my $sth = $dbhvar->prepare($sql);

  open my $load_vs, "<", "$load_file" or die "Can not open $load_file: $!";
  while (<$load_vs>) {
    chomp;
    my $var_id = (split)[0];
    my $var_set_id = (split)[1];
    
    $sth->execute($var_id, $var_set_id);
  }
  close $load_vs;
  $sth->finish();

}

sub update_variation_feature_table {
  # this function after populating the variation_feature_backup table created by inserting from the original table would then update the variation_set_id column uaing the file from the recalculate
  # using the parents and parents set id to update the variation feature table 
  my $dbhvar = shift; 
  my $load_file = shift;


  my $update_temp_vf = $dbhvar->prepare(q{ UPDATE variation_feature SET variation_set_id = ? 
                                          WHERE variation_id = ? });
  
  #my %var_data;
  open my $load_fh, "<", "$load_file" or die "Can not open $load_file: $!";
  while (<$load_fh>) {
    chomp;
    my @fields = split("\t");
    my $var_id = $fields[0];
    my $var_set_id = $fields[1];
    

    my @sets_array; 

    my @sets_array = unique(split(',', $var_set_id));

    my @sorted_array = sort { $a<=>$b } @sets_array;
    my $values = join(',', @sorted_array);
    $values =~ s/\s*,\s*/,/g; # to eliminate spaces and stuff 
    $values =~ s/^\s+//;  #to eliminate spaces and stuff 

    $update_temp_vf->execute($values, $var_id); # creating a hash which has the var_id has the key and the set var_set_id has the values 
  }
  
  close $load_fh;
  $update_temp_vf->finish();
}

sub dump_old_sql_variation_sets {
  # this would dump old variation_sets table from the old database and create a file called variation_name_set.txt with two columns var_id and variation_set_id
  my $dbhvar = shift;
  my $tmp_num = shift;
  my $chunk = shift;
  my $size = shift;

  my $start = $chunk * $tmp_num;
  my $end = $chunk + $start;
  $end = $end < $size ? $end : $size;


  my $sql = qq{SELECT v.name, vs.variation_set_id from variation_set_variation vs LEFT JOIN variation v
               ON v.variation_id = vs.variation_id where vs.variation_set_id != 1 AND v.source_id = 1 AND v.variation_id > $start AND v.variation_id <= $end  };
  my $sth = $dbhvar->prepare($sql);

  open (my $dump_fh, ">>$TMP_DIR/$tmp_vset" )
      or die( "Cannot open $TMP_DIR/$tmp_vset: $!" );

  $sth->execute();

  while ( my $aref = $sth->fetchrow_arrayref() ) {
    my @a = map {defined($_) ? $_ : '\N'} @$aref;
    print $dump_fh join("\t", @a), "\n";
  }
  
  $sth->finish();
  close $dump_fh;

}
 
sub dump_new_variation_feature {
  # this would dump the variation_id and name from the variation_feature table creating a file called variation_feature.txt with two columns contains var_id and name
  my $dbh = shift;
  my $tmp_num = shift;
  my $chunk = shift;
  my $size = shift;

  my $start = $chunk * $tmp_num;
  my $end = $chunk + $start;
  $end = $end < $size ? $end : $size;


  my $sql = qq{SELECT DISTINCT(variation_id), variation_name from variation_feature where variation_id > $start AND variation_id <= $end };
  my $sth = $dbh->prepare($sql);
  
  open (my $dump_var, ">>$TMP_DIR/$new_vf_file" )
      or die( "Cannot open $TMP_DIR/$new_vf_file: $!" );

  $sth->execute();

  while ( my $aref = $sth->fetchrow_arrayref() ) {
    my @a = map {defined($_) ? $_ : '\N'} @$aref;
    print $dump_var join("\t", @a), "\n";
  }
  
  $sth->finish();
  close $dump_var;

  
}

sub get_structure {
  my $dbhvar = shift;

  my $get_struc_sth = $dbhvar->prepare(qq[ select variation_set_sub, variation_set_super from variation_set_structure]);

  my %parent;
  $get_struc_sth->execute() ||die;
  my $dat = $get_struc_sth->fetchall_arrayref();
  foreach my $l(@{$dat}){
    $parent{$l->[0]} = $l->[1] ;
  }
  return \%parent;

}

sub recalculate {
  my $input_file = shift;
  my $output_file = shift;
  
  my $parent = get_structure($dbh);
  my %concat_sets;
  
  open my $re_input, "<", "$input_file" or die "Can not open $input_file: $!";
  

  while (<$re_input>) {
    chomp;
    my @fields = split("\t");
    my $var_id = $fields[0];
    my $var_set_id = $fields[1];
    my @sets;
    
    push @{$concat_sets{$var_id}}, $var_set_id;
    push @{$concat_sets{$var_id}}, $parent->{$var_set_id} if exists  $parent->{$var_set_id}; #pushing parents and var_set_id 
    push @{$concat_sets{$var_id}}, $parent->{$parent->{$var_set_id}} if exists $parent->{$parent->{$var_set_id}}; 
  
  }


  open(my $re_output, '>', $output_file) or die "Could not open file '$output_file': $!";
  foreach my $var_id (keys %concat_sets) {  
    my $values_str = join(", ", @{$concat_sets{$var_id}});
    print $re_output "$var_id\t$values_str\n"; # adding the values str to it 
  }
  
  # Close the file
  close $re_output;
  close $re_input;
}

sub unique {
 my @array = shift;
 
 my %u;
  map { $u{$_} = 1; } @_;
  return keys %u;

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
  
  # to define temporary directory if it is not defined 
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

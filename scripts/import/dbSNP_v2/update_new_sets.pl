use strict;
use DBI;
use Bio::EnsEMBL::Registry;
use Getopt::Long;
use Cwd qw(cwd);


my ($reg_file, $old_reg_file, $release, $tmp, $help);

GetOptions ("registry=s"        => \$reg_file,
            "old_registry:s"    => \$old_reg_file,
            "release=s"         =>  \$release,
            "tmp=s"             =>  \$tmp,
            "help|h"           =>   \$help,
);

# to define temporary directory if it is not defined 
if (!defined $tmp){
  $tmp = cwd;
}
usage() unless defined $reg_file && defined $release && defined $tmp;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($reg_file);

my $db_adaptor = $registry->get_DBAdaptor("homo_sapiens", "variation");
my $dbh = $db_adaptor->dbc;

my $TMP_DIR = $tmp;
my $TMP_FILE = "variation_feature.txt";
my $tmp_vset = "variation_name_set.txt";
my $tmp_merged = "new_variation_set.txt";
my $tmp_vs_file = "vset_concat.txt";
my $old_dbh;

# define variants and move through the list
#if ($old_reg_file) {
 # $registry->load_all($old_reg_file);
  #$db_adaptor = $registry->get_DBAdaptor("homo_sapiens", "variation");
  # Connect to DBI
  #$old_dbh = $db_adaptor->dbc;
#} else {
 # my $old_release = $release - 1;
  #my $old_dbname = $dbh->dbname =~ s/_${release}_/_${old_release}_/gr;
  #my $old_host = $dbh->host;
  #my $old_port = $dbh->port;
  # Connect to DBI
  #$old_dbh = DBI->connect("DBI:mysql:database=${old_dbname};host=${old_host};port=${old_port}", "ensro", "");
#}

#my $sql = qq{ SELECT MIN(variation_id), MAX(variation_id) FROM variation };
#my $sth = $old_dbh->prepare( $sql );
#$sth->execute();
#my $vf = $sth->fetchall_arrayref();
#my $min_id = $vf->[0]->[0];
#my $chunk = 1000000;
##my $max_id = $vf->[0]->[1];
#for my $tmp_num (map { $_ } $min_id/$chunk .. $max_id/$chunk) {
 # dump_old_sql_variation_sets($old_dbh, $tmp_num, $chunk, $max_id);
  #dump_new_variation_feature($dbh, $tmp_num, $chunk, $max_id);
#}

#system(sort -u $TMP_FILE );
#system(sort -u $tmp_vset);

#create_merged_file($tmp_vset, $TMP_FILE, $tmp_merged);



#temp_table($dbh);
#load_all_variation_sets($dbh, $tmp_merged);
update_variation_feature_table($dbh, $tmp_vs_file);




sub temp_table { 
  my $dbhvar = shift; 
   
  # creating a temp table and altering the table by disabling the keys 
  my $create_sql = $dbhvar->prepare(qq{CREATE TABLE if not exists temp_variation_set LIKE variation_set_variation});
  my $alter_sql = $dbhvar->prepare(qq{ALTER TABLE temp_variation_set DISABLE keys});
  
  #my $create_backup_vf = $dbhvar->prepare(q{CREATE table if not exists variation_feature_bk like variation_feature});

  #executing the sql 
  $create_sql->execute();
  $alter_sql->execute();

  $create_sql->finish();
  $alter_sql->finish();

  #$create_backup_vf->execute();
  #$create_backup_vf->finish();
}

sub create_merged_file {
  my $file = shift; 
  my $second_file = shift;
  my $third_file = shift;

  open FILE1, "<", "$file" or die "Cannot open $file: $!";
  open FILE2, "<", "$second_file" or die "Cannot open $second_file: $!";
  open OUTPUT, ">", "$third_file" or die "Cannot open $third_file: $!";
  my %data;

  while (<FILE1>) {
    chomp;
    my $key = (split)[0]; # to access the variation name
    my @value = (split)[1]; # to access set ids
    if (exists $data{$key}) {
       push @{$data{$key}}, \@value; # push the array reference to the hash
    } else {
      push $data{$key} = [\@value];
    }
  }
  while (<FILE2>) {
    chomp;
    my $var_name = (split)[1];
    my $var_id = (split)[0];
    if(exists $data{$var_name}) {
      my $array_ref = $data{$var_name}; # creating an array reference if the var_name exists in the created hash 
      foreach my $set_id (@$array_ref) {
        print OUTPUT "$var_id\t@$set_id[0]\n"; # to access the values in the array reference
      }
    }
  }
  close FILE1;
  close FILE2;
  close OUTPUT;

}

sub load_all_variation_sets {
  my $dbhvar = shift;
  my $load_file = shift;
  
  my $sql = qq{INSERT INTO temp_variation_set (variation_id, variation_set_id ) VALUES (?, ?)};
  my $sth = $dbhvar->prepare($sql);
  open FH, "<", "$load_file" or die "Can not open $load_file: $!";
  while (<FH>) {
    chomp;
    my $var_id = (split)[0];
    my $var_set_id = (split)[1];
    
    $sth->execute($var_id, $var_set_id);
  }
  close FH;
  $sth->finish();

  #dump variation_sets_id to populate the variation_sets table 

  my $dump_vs = $dbhvar->prepare(qq{  SELECT variation_id, GROUP_CONCAT(DISTINCT(variation_set_id)) FROM temp_variation_set  GROUP BY variation_id;})
}

sub update_variation_feature_table {
  my $dbhvar = shift; 
  my $load_file = shift;
  
  my $insert_temp_vf = $dbhvar->prepare(q{ INSERT into variation_feature_bk SELECT * from variation_feature });
  
  $insert_temp_vf->execute();
  $insert_temp_vf->finish();
  
  my $update_temp_vf = $dbhvar->prepare(q{ UPDATE variation_feature_bk SET variation_set_id = ? 
                                          WHERE variation_id = ?});
  
  my %var_data;
  open FH, "<", "$load_file" or die "Can not open $load_file: $!";
  while (<FH>) {
    chomp;
    my $var_id = (split)[0];
    my $var_set_id = (split)[1];
    
    $var_data{$var_id} = $var_set_id;
  }
  close FH;

  foreach my $key (keys %var_data){
    $update_temp_vf->execute($var_data{$key}, $key);
  }
  $update_temp_vf->finish();

}

sub dump_new_variation_sets {

  my $dump_vs = $dbh->prepare(qq{  SELECT variation_id, GROUP_CONCAT(DISTINCT(variation_set_id)) FROM variation_set_variation GROUP BY variation_id});
  open (FH, ">>$TMP_DIR/$tmp_vs_file" )
      or die( "Cannot open $TMP_DIR/$tmp_vs_file: $!" );
  $dump_vs->execute();

  while ( my $aref = $dump_vs->fetchrow_arrayref() ) {
    my @a = map {defined($_) ? $_ : '\N'} @$aref;
    print FH join("\t", @a), "\n";
  }
  
  $dump_vs->finish();
  close FH;
}
sub dump_old_sql_variation_sets {
  my $dbhvar = shift;
  my $tmp_num = shift;
  my $chunk = shift;
  my $size = shift;

  my $start = $chunk * $tmp_num;
  my $end = $chunk + $start;
  $end = $end < $size ? $end : $size;


  my $sql = qq{SELECT v.name, vs.variation_set_id from variation_set_variation vs LEFT JOIN variation v
               ON v.variation_id = vs.variation_id where variation_set_id != 1 AND v.variation_id > $start AND v.variation_id <= $end LIMIT 50 };
  my $sth = $dbhvar->prepare($sql);

  local *FH;
  open (FH, ">>$TMP_DIR/$tmp_vset" )
      or die( "Cannot open $TMP_DIR/$tmp_vset: $!" );

  $sth->execute();

  while ( my $aref = $sth->fetchrow_arrayref() ) {
    my @a = map {defined($_) ? $_ : '\N'} @$aref;
    print FH join("\t", @a), "\n";
  }
  
  $sth->finish();
  close FH;

}
 
sub dump_new_variation_feature {
  my $dbh = shift;
  my $tmp_num = shift;
  my $chunk = shift;
  my $size = shift;

  my $start = $chunk * $tmp_num;
  my $end = $chunk + $start;
  $end = $end < $size ? $end : $size;


  my $sql = qq{SELECT variation_id, variation_name from variation_feature where variation_id > $start AND variation_id <= $end LIMIT 50};
  my $sth = $dbh->prepare($sql);
  
  local *FH;
  open (FH, ">>$TMP_DIR/$TMP_FILE" )
      or die( "Cannot open $TMP_DIR/$TMP_FILE: $!" );

  $sth->execute();

  while ( my $aref = $sth->fetchrow_arrayref() ) {
    my @a = map {defined($_) ? $_ : '\N'} @$aref;
    print FH join("\t", @a), "\n";
  }
  
  $sth->finish();
  close FH;

  
}

sub usage {

  die "\n\tUsage: update_old_sets.pl -registry [registry file] -release [release number] \tOptional:  -tmp [temp folder] or gets set based on current directory 
  -old_registry [old database registry file]\n\n";
}
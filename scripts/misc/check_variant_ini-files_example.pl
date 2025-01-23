# Copyright [2018-2025] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Script to generate an HTML page containing the variation sources of each species


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut


use Bio::EnsEMBL::Registry;
use DBI;
use strict;
use POSIX;
use Getopt::Long;

###############
### Options ###
###############
my ($e_version,$hlist,$user,$public_plugins_git_dir,$help);

usage() if (!scalar(@ARGV));
 
GetOptions(
     'v=s'       => \$e_version,
     'help!'     => \$help,
     'hlist=s'   => \$hlist,
     'user=s'    => \$user,
     'git_dir=s' => \$public_plugins_git_dir
);

if (!$e_version) {
  print "> Error! Please give an Ensembl version, using the option '-v' \n";
  usage();
}
if (!$hlist) {
  print "> Error! Please give the list of host names where the new databases are stored using the option '-hlist'\n";
  usage();
}
if (!$user) {
  print "> Error! Please give user name using the option '-user'\n";
  usage();
}
if (!$public_plugins_git_dir) {
  print "> Error! Please give the path to your clone of the public_plugins Git repo, using the option '-git_dir'\n";
  usage();
}

usage() if ($help);

my @hostnames = split /,/, $hlist;

# Settings
my $database = "";
my $ini_files_dir = 'ensembl/conf/ini-files';
my $pswd = "";
my $db_type = 'variation';
my $error_char = '******';

# Git dir
chdir "$public_plugins_git_dir/$ini_files_dir";
`git checkout master`;
`git pull origin master`;
print "Git dir: $public_plugins_git_dir/$ini_files_dir\n";

my $sql   = qq{SHOW DATABASES LIKE '%$db_type\_$e_version%'};

my $sql_var_2  = qq{SELECT variation_id FROM variation where display=1 and variation_id=?};
my $sql_sv_2   = qq{SELECT structural_variation_id FROM structural_variation where structural_variation_id=?
                       AND structural_variation_id NOT IN (SELECT structural_variation_id FROM failed_structural_variation)
                   };
my $sql_phe_2  = qq{SELECT description FROM phenotype where phenotype_id=?};

my $sql_var_3  = qq{SELECT variation_id FROM variation where variation_id=?};
my $sql_sv_3   = qq{SELECT structural_variation_id FROM structural_variation where structural_variation_id=?};
my $sql_phe_3  = qq{SELECT count(*) FROM phenotype_feature where phenotype_id=?};

my $V_PARAM = 'VARIATION_PARAM';
my $V_TEXT  = 'VARIATION_TEXT';

my $SV_PARAM = 'STRUCTURAL_PARAM';
my $SV_TEXT  = 'STRUCTURAL_TEXT';

my $P_PARAM = 'PHENOTYPE_PARAM';
my $P_TEXT  = 'PHENOTYPE_TEXT';

foreach my $hostname (@hostnames) {
  
  my $sth = get_connection_and_query($database, $hostname, $sql);
  
  # Loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname =~ /^master_schema/ || $dbname =~ /private/ || $dbname =~ /_variation_\d+_\d+_\w+$/ );
    
    print STDERR "$dbname\n";

    $dbname =~ /^(.+)_variation/;
    my $s_name = $1;
    
    print "\n# $s_name:\n";
    
    
    #### Test Variant Entries ####
    my $var_report = '  - Variation: ';
    
    # VARIATION_PARAM
    my $var_param = `grep $V_PARAM ./$s_name.ini`;
    # VARIATION_TEXT
    my $var_text =`grep $V_TEXT ./$s_name.ini`;
    
    if ($var_param && $var_text) {
      # Parse data
      $var_param =~ /$V_PARAM\s+=\s+(\w+)\s*$/;
      my $var_id = $1;

      $var_text =~ /$V_TEXT\s+=\s+(\w+)\s*$/;
      my $var_desc = $1;
      
      if ($var_id eq $var_desc) {
        # Check variant ID exist and is displayable
        my $sth_var_2 = get_connection_and_query($dbname, $hostname, $sql_var_2, [$var_id]);
        if ($sth_var_2) {
          $var_report .= "ID found and visible [ID: $var_id]";
        }
        else {
          my $sth_var_3 = get_connection_and_query($dbname, $hostname, $sql_var_3, [$var_id]);
          if ($sth_var_3) {
            $var_report .= "The ID $var_id has been flagged as failed! $error_char";
          }
          else {
            $var_report .= "Can't find an entry in the DB for the ID $var_id! $error_char";
          }
          $sth_var_3->finish();
        }
        $sth_var_2->finish();
      }
      else {
         $var_report .= "$var_id ($V_PARAM) and $var_desc ($V_TEXT) are different! $error_char";
      }
    }
    elsif ($var_param || $var_text) {
        $var_report .= "missing attribute $V_PARAM in the 'ini' file! $error_char" if (!$var_param);
        $var_report .= "missing attribute $V_TEXT in the 'ini' file! $error_char"  if (!$var_text);
      }
    else {
      $var_report .= "no variant entry in the 'ini' file for this species! $error_char";
    }
    print "$var_report\n";
    
    
    
    #### Test Structural Variant Entries ####
    my $sv_report = '  - SV: ';
    
    # STRUCTURAL_PARAM
    my $sv_param = `grep $SV_PARAM ./$s_name.ini`;
    # STRUCTURAL_TEXT
    my $sv_text =`grep $SV_TEXT ./$s_name.ini`;
    
    if ($sv_param && $sv_text) {
      # Parse data
      $sv_param =~ /$SV_PARAM\s+=\s+(\w+)\s*$/;
      my $sv_id = $1;
      
      $sv_text =~ /$SV_TEXT\s+=\s+(\w+)\s*$/;
      my $sv_desc = $1;
      
      if ($sv_id eq $sv_desc) {
        # Check structural variant ID exist and is not flagged as failed
        my $sth_sv_2 = get_connection_and_query($dbname, $hostname, $sql_sv_2, [$sv_id]);
        if ($sth_sv_2) {
          $sv_report .= "ID found and visible [ID: $sv_id]";
        }
        else {
          my $sth_sv_3 = get_connection_and_query($dbname, $hostname, $sql_sv_3, [$sv_id]);
          if ($sth_sv_3) {
            $sv_report .= "The ID $sv_id has been flagged as failed! $error_char";
          }
          else {
            $sv_report .= "Can't find an entry in the DB for the ID $sv_id! $error_char";
          }
          $sth_sv_3->finish();
        }
        $sth_sv_2->finish();
      }
      else {
         $sv_report .= "$sv_id ($SV_PARAM) and $sv_desc ($SV_TEXT) are different! $error_char";
      }
    }
    elsif ($sv_param || $sv_text) {
        $sv_report .= "missing attribute $SV_PARAM in the 'ini' file! $error_char" if (!$sv_param);
        $sv_report .= "missing attribute $SV_TEXT in the 'ini' file! $error_char"  if (!$sv_text);
      }
    else {
      $sv_report .= "no structural variant entry in the 'ini' file for this species";
    }
    print "$sv_report\n";
    
    
    
    #### Test Phenotype Entries ####
    my $phe_report = '  - Phenotype: ';
    
    # PHENOTYPE_PARAM
    my $phe_param = `grep $P_PARAM ./$s_name.ini`;
    # PHENOTYPE_TEXT
    my $phe_text =`grep $P_TEXT ./$s_name.ini`;
    
    if ($phe_param && $phe_text) {
      # Parse data
      $phe_param =~ /$P_PARAM\s+=\s+(\w+)\s*$/;
      my $phe_id = $1;
      
      $phe_text =~ /$P_TEXT\s+=\s+(\w.+)\s*$/;
      my $phe_desc = $1;
      
      # Check phenotype ID exist and phenotype description matches
      my $sth_phe_2 = get_connection_and_query($dbname, $hostname, $sql_phe_2, [$phe_id]);
      if ($sth_phe_2) {
        my $db_phe_desc = $sth_phe_2->fetchrow_array;
        $sth_phe_2->finish();
        if ($phe_desc eq $db_phe_desc) {
            my $sth_phe_3 = get_connection_and_query($dbname, $hostname, $sql_phe_3, [$phe_id]);
            my $count = $sth_phe_3->fetchrow_array;
            $sth_phe_3->finish();
            $phe_report .= "ID and description match [ID: $phe_id | DESC: $phe_desc | COUNT: $count]";
        }
        else {
          $phe_report .= "description doesn't match! [ID:$phe_id | INI: $phe_desc | DB: $db_phe_desc] $error_char";
        }
      }
      else {
        $phe_report .= "Can't find an entry in the DB for the ID $phe_id! $error_char";
      }
    }
    elsif ($phe_param || $phe_text) {
      $phe_report .= "missing attribute $P_PARAM in the 'ini' file! $error_char" if (!$phe_param);
      $phe_report .= "missing attribute $P_TEXT in the 'ini' file! $error_char"  if (!$phe_text);
    }
    else {
      $phe_report .= "no phenotype entry in the 'ini' file for this species";
    }
    print "$phe_report\n";
  }
}


# Connects and execute a query
sub get_connection_and_query {
  my $dbname = shift;
  my $hname  = shift;
  my $sql    = shift;
  my $params = shift;
  
  my ($host, $port) = split /\:/, $hname;

  # DBI connection 
  my $dsn = "DBI:mysql:$dbname:$host:$port";
  my $dbh = DBI->connect($dsn, $user, $pswd) or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  if ($params) {
    $sth->execute(join(',',@$params));
  }
  else {
    $sth->execute();
  }
  return $sth;
}


sub usage {
  
  print qq{
  Usage: perl check_variant_ini-files_example.pl [OPTION]
  
  Check variant, structural variant and phenotype entries in the ini-files (public_plugins)
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 92 (Required)   
    -hlist          The list of host names (with port) where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1:1234, ensembldb.ensembl.org2:1234 (Required)
    -user           MySQL user name (Required)
    -git_dir        Path to your local copy of the  public_plugins Git repository, to check the ini-files,
                    e.g. ~/git/public_plugins (Required)

  } . "\n";
  exit(0);
}


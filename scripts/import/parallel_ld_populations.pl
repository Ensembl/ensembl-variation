# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use strict;


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use warnings;

use Getopt::Long;
use Fcntl ':flock';
use DBI;
use DBH;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);

use ImportUtils qw(debug load);

my $default_population = 'CSHL-HAPMAP:HapMap-CEU'; #default population used to show by default on ldview

my ($TMP_DIR, $TMP_FILE, $LIMIT,$status_file, $population, $num_processes); #ld_global is a flag to know wether we run the ld calculation globally or by population

{
  my ($vhost, $vport, $vdbname, $vuser, $vpass,
      $chost, $cport, $cdbname, $cuser, $cpass,
      $limit);

  GetOptions('chost=s'   => \$chost,
             'cuser=s'   => \$cuser,
             'cpass=s'   => \$cpass,
             'cport=i'   => \$cport,
             'cdbname=s' => \$cdbname,
             'vhost=s'   => \$vhost,
             'vuser=s'   => \$vuser,
             'vpass=s'   => \$vpass,
             'vport=i'   => \$vport,
             'vdbname=s' => \$vdbname,
             'tmpdir=s'  => \$ImportUtils::TMP_DIR,
             'tmpfile=s' => \$ImportUtils::TMP_FILE,
             'limit=s'   => \$limit,
	     'population=i' => \$population,
	     'status_file=s' => \$status_file,
	     'num_processes=i' => \$num_processes
	     );

  #added default options
  $chost    ||= 'ecs2';
  $cuser    ||= 'ensro';
  $cport    ||= 3365;

  $vport    ||= 3306;
  $vuser    ||= 'ensadmin';
  
  $LIMIT = ($limit) ? " LIMIT $limit " : '';

  usage('-vdbname argument is required') if(!$vdbname);
  usage('-cdbname argument is required') if(!$cdbname);

 # usage('-num_processes must at least be 1') if ($num_processes == 0);
#  usage('-status_file argument is required') if (!$status_file);

  my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host   => $chost,
     -user   => $cuser,
     -pass   => $cpass,
     -port   => $cport,
     -dbname => $cdbname);

  my $dbVar = DBH->connect
    ("DBI:mysql:host=$vhost;dbname=$vdbname;port=$vport",$vuser, $vpass );
  die("Could not connect to variation database: $!") if(!$dbVar);


  $TMP_DIR  = $ImportUtils::TMP_DIR;
  $TMP_FILE = $ImportUtils::TMP_FILE;

  #if is the last process, delete the variation table and upload with the new coordinate information from the different tables
  last_process($dbCore,$dbVar);

}

#will know if it is the last process running counting the lines of the status_file.If so, load all data
sub last_process{
    my $dbCore = shift;
    my $dbVar = shift;
    debug("Importing pairwise data");
   #group all the fragments in 1 file
    my $dbname = $dbVar->dbname(); #get the name of the database to create the file
    my $call = "cat $TMP_DIR/$dbname.pairwise_ld*out* > $TMP_DIR/$TMP_FILE";
    system($call);
    unlink(<$TMP_DIR/$dbname.pairwise_ld_*>);    

    #and import the data in the database
    load($dbVar, qw(pairwise_ld variation_feature_id_1 variation_feature_id_2 sample_id seq_region_id seq_region_start seq_region_end r2 d_prime sample_count));

    update_meta_coord($dbCore, $dbVar, 'pairwise_ld');

    update_meta($dbVar);
    #and delete the status file
    unlink("$TMP_DIR/$status_file");
}

#
# updates the meta coord table
#
sub update_meta_coord {
  my $dbCore = shift;
  my $dbVar  = shift;
  my $table_name = shift;
  my $csname = shift || 'chromosome';

  my $csa = $dbCore->get_CoordSystemAdaptor();

  my $cs = $csa->fetch_by_name($csname);

  my $sth = $dbVar->prepare
    ('INSERT INTO meta_coord set table_name = ?, coord_system_id = ?');

  $sth->execute($table_name, $cs->dbID());

  $sth->finish();

  return;
}

sub update_meta{
    my $dbVar = shift;

    my $population_id;
    my $sth_pop = $dbVar->prepare('SELECT s.sample_id FROM sample s, population p where p.sample_id = s.sample_id AND s.name = ?');
    
    $sth_pop->execute($default_population);
    $sth_pop->bind_columns(\$population_id);
    $sth_pop->fetch();
    $sth_pop->finish();

    my $sth = $dbVar->prepare(qq{INSERT INTO meta (meta_key,meta_value) VALUES (?,?)
				 });
    $sth->execute('pairwise_ld.default_population', $population_id);

    $sth->finish();
    return;
}

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl ld_populations.pl <options>

options:
    -chost <hostname>    hostname of core Ensembl MySQL database (default = ecs2)
    -cuser <user>        username of core Ensembl MySQL database (default = ensro)
    -cpass <pass>        password of core Ensembl MySQL database
    -cport <port>        TCP port of core Ensembl MySQL database (default = 3364)
    -cdbname <dbname>    dbname of core Ensembl MySQL database
    -vhost <hostname>    hostname of variation MySQL database to write to
    -vuser <user>        username of variation MySQL database to write to (default = ensadmin)
    -vpass <pass>        password of variation MySQL database to write to
    -vport <port>        TCP port of variation MySQL database to write to (default = 3306)
    -vdbname <dbname>    dbname of variation MySQL database to write to
    -limit <number>      limit the number of rows for testing
    -tmpdir <dir>        temp directory to use (with lots of space!)
    -tmpfile <filename>   name of temp file to use
    -num_processes <number> number of processes that are running (default = 1)
    -status_file <filename> name of a temp file where all the processes write when they finish
EOF

  die("\n$msg\n\n");
}

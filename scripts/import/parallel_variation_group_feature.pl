use strict;
use warnings;

use Getopt::Long;
use DBI;
use DBH;
use ImportUtils qw(debug);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);

my ($TMP_DIR, $TMP_FILE, $LIMIT);

{
  my ($vhost, $vport, $vdbname, $vuser, $vpass,
      $chost, $cport, $cdbname, $cuser, $cpass,
      $limit, $num_processes);

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
             'limit=s'   => \$limit);

  #added default options
  $chost    ||= 'ecs2';
  $cuser    ||= 'ensro';
#  $cport    ||= 3364;

  $vport    ||= 3306;
  $vuser    ||= 'ensadmin';
  
  $num_processes ||= 1;

  $LIMIT = ($limit) ? " LIMIT $limit " : '';

  usage('-vdbname argument is required') if(!$vdbname);
  usage('-cdbname argument is required') if(!$cdbname);


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

  variation_group_feature($dbCore, $dbVar);
}

# Loads the variation_group_feature table by using the locations of
# the variation_features and their groupings in the variation_group_feature
# table
#

sub variation_group_feature {
  my $dbCore = shift;
  my $dbVar  = shift;

  ### TBD: possibly fix this: does not check to see that all the
  ### variation_features are actually on the same seq_region or in close
  ### proximity

  debug("Loading variation_group_feature table");

  $dbVar->do(qq{INSERT INTO variation_group_feature 
                       (seq_region_id, seq_region_start, seq_region_end,
                        seq_region_strand, variation_group_id, variation_group_name)
                SELECT vf.seq_region_id, min(vf.seq_region_start),
                       max(vf.seq_region_end), vf.seq_region_strand,
                       vgv.variation_group_id, vg.name
                FROM   variation_feature vf, variation_group_variation vgv, variation_group vg
                WHERE  vgv.variation_id = vf.variation_id
		AND    vg.variation_group_id = vgv.variation_group_id
                GROUP BY vgv.variation_group_id, vf.seq_region_id});


  update_meta_coord($dbCore, $dbVar, 'variation_group_feature');

  return;
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

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl parallel_variation_group_feature.pl <options>

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
EOF

  die("\n$msg\n\n");
}

#!/usr/bin/env perl
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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


use strict;
use warnings;

use Getopt::Long;
use DBH;
use FindBin qw( $Bin );

use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use ImportUtils qw(debug load);

my ($TMP_DIR, $TMP_FILE); #global variables for the tmp files and folder

my ($vhost, $vport, $vdbname, $vuser, $vpass,
    $chost, $cport, $cdbname, $cuser, $cpass,
    $species, $read_dir, $max_level,
    $ind_file, $bsub_wait_status);

  GetOptions(#'chost=s'     => \$chost,
             #'cuser=s'     => \$cuser,
             #'cpass=s'     => \$cpass,
             #'cport=i'     => \$cport,
             #'cdbname=s'   => \$cdbname,
             #'vhost=s'     => \$vhost,
             'vuser=s'     => \$vuser,
             #'vpass=s'     => \$vpass,
             #'vport=i'     => \$vport,
             #'vdbname=s'   => \$vdbname,
	     'species=s'     => \$species,
             'tmpdir=s'    => \$ImportUtils::TMP_DIR,
             'tmpfile=s'   => \$ImportUtils::TMP_FILE,
             'maxlevel=i'  => \$max_level,
	     'readdir=s' => \$read_dir,
	     'indfile=s' => \$ind_file,
             'bsubwaitstatus=s' => \$bsub_wait_status);

warn("Make sure you have a updated ensembl.registry file!\n");

my $registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');


my $dbVar = $vdba->dbc->db_handle;
my $dbCore = $cdba;

#added default options
$chost = $cdba->dbc->host; 
$cuser    ||= 'ensro';
$cport = $cdba->dbc->port ;
$cdbname = $cdba->dbc->dbname;

$vhost = $vdba->dbc->host;
$vport = $vdba->dbc->port;
$vuser    ||= 'ensadmin';
$vdbname = $vdba->dbc->dbname;
$vpass = $vdba->dbc->password;

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

$bsub_wait_status ||= 'done';

my $call;
my $i = 0;
foreach my $read_file (glob("$read_dir/*.mapped")){
    $call = "bsub -J read_coverage_job_$i -o $TMP_DIR/output_reads_$i.txt -q research /sw/arch/bin/env perl read_coverage.pl -chost $chost -cuser $cuser -cport $cport -cdbname $cdbname -vhost $vhost -vuser $vuser -vpass $vpass -vport $vport -vdbname $vdbname -tmpdir $TMP_DIR -tmpfile read_coverage_$i.txt -maxlevel $max_level -readfile $read_file";

    if ($ind_file)
    {
        $call .= " -indfile $ind_file";
    }

    system($call);
    $i++;
}
debug("Waiting for read_coverage jobs to finish");
$call = "bsub -K -w '$bsub_wait_status(read_coverage_job_*)' -J waiting_process sleep 1"; #waits until all reads have succesfully finished
system($call);
 debug("Ready to import coverage data");
 $call = "cat $TMP_DIR/read_coverage_* > $TMP_DIR/$TMP_FILE"; #concat all files 
 system($call);

load($dbVar,"read_coverage","seq_region_id","seq_region_start","seq_region_end","level","sample_id"); #load table with information

unlink(<$TMP_DIR/read_coverage_*>); #and delete the files with the data
update_meta_coord($dbCore,$dbVar,'read_coverage');


#
# Updates the meta coord table
#
sub update_meta_coord {
  my $dbCore = shift;
  my $dbVar  = shift;
  my $table_name = shift;
  my $csname = shift || 'chromosome';

  my $csa = $dbCore->get_CoordSystemAdaptor();

  my $cs = $csa->fetch_by_name($csname);
  my $max_length_ref = $dbVar->selectall_arrayref(qq{SELECT max(seq_region_end - seq_region_start+1) from read_coverage});
  my $max_length = $max_length_ref->[0][0];
  my $sth = $dbVar->prepare
    ('INSERT INTO meta_coord set table_name = ?, coord_system_id = ?, max_length = ?');

  $sth->execute($table_name, $cs->dbID(), $max_length);

  $sth->finish();

  return;
}


sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl read_coverage.pl <options>

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
    -tmpdir <dir>        temp directory to use (with lots of space!)
    -tmpfile <filename>  name of temp file to use
    -readdir <dirname>   name of dir with read information
EOF

  die("\n$msg\n\n");
}

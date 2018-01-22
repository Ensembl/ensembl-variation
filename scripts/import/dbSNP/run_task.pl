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
use Progress;
use Fcntl qw( LOCK_SH LOCK_EX );
use dbSNP::ImportTask;
use dbSNP::DBManager;

#This script is used to launch a small subtask during the import process. Typically, a bigger task is broken into pieces and submitted to the farm using this script

my @option_defs = (
  'species=s',
  'dbSNP_shared=s',
  'registry_file=s',
  'task=s',
  'task_parameters=s',
  'task_management_file=s',
  'chunksize=i',
  'file_prefix=s',
  'logfile=s',
  'tempdir=s',
  'tempfile=s',
  'source_engine:s',
  'schema_name:s',
);

my %options;
GetOptions(\%options,@option_defs);

my $species = $options{'species'};
my $dbSNP_shared = $options{'dbSNP_shared'};
my $registryfile = $options{'registry_file'};
my $task = $options{'task'};
my $task_parameters = $options{'task_parameters'};
my $task_management_file = $options{'task_management_file'};
my $chunksize = $options{'chunksize'};
my $file_prefix = $options{'file_prefix'};
my $logfile = $options{'logfile'};
my $tempdir = $options{'tempdir'};
my $tempfile = $options{'tempfile'};
my $source_engine = $options{'source_engine'};
my $schema_name = $options{'schema_name'};

$ImportUtils::TMP_DIR = $tempdir;
$ImportUtils::TMP_FILE = $tempfile;

my $dbm = dbSNP::DBManager->new(
    $registryfile,
    $species,
    $schema_name 
);
$dbm->dbSNP_shared($dbSNP_shared);

my $logh = *STDOUT;
if (defined($logfile)) {
  open(LOG,'>>',$logfile);
  #Turn on autoflush for the logfile
  {
    my $ofh = select LOG;
    $| = 1;
    select $ofh;
  }
  $logh = *LOG;
}

#If we are a job array element submitted to the farm and the task_start and task_end parameters haven't been set, we refer to the task management file to get our task parameters
if (defined($ENV{'LSB_JOBINDEX'}) && -e $task_management_file) {
    my $jobid = $ENV{'LSB_JOBINDEX'};
    open(FH,'<',$task_management_file);
    flock(FH,LOCK_SH);
    while (<FH>) {
        chomp;
        my @row = split;
        if (shift(@row) == $jobid) {
            $task_parameters = join("\t",@row);
            last;
        }
    }
    close(FH);
}

my $task_obj = dbSNP::ImportTask->new(
    $dbm,
    $task,
    $file_prefix,
    undef,
    $logh
);
my @args = split(/\s/,$task_parameters);
push @args, $source_engine if  $task =~ /calculate_gtype/;  ## mysql or mssql syntax
$task_obj->$task(@args);

close(LOG) if (defined($logfile));

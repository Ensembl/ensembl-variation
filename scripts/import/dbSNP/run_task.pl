#!/usr/local/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Progress;
use Fcntl qw( LOCK_SH LOCK_EX );
use dbSNP::ImportTask;
use dbSNP::DBManager;

#�This script is used to launch a small subtask during the import process. Typically, a bigger task is broken into pieces and submitted to the farm using this script

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
  'tempfile=s'
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

$ImportUtils::TMP_DIR = $tempdir;
$ImportUtils::TMP_FILE = $tempfile;

my $dbm = dbSNP::DBManager->new(
    $registryfile,
    $species
);
$dbm->dbSNP_shared($dbSNP_shared);

my $logh = *STDOUT;
if (defined($logfile)) {
  open(LOG,'>>',$logfile);
  #�Turn on autoflush for the logfile
  {
    my $ofh = select LOG;
    $| = 1;
    select $ofh;
  }
  $logh = *LOG;
}

#�If we are a job array element submitted to the farm and the task_start and task_end parameters haven't been set, we refer to the task management file to get our task parameters
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
$task_obj->$task(@args);

close(LOG) if (defined($logfile));

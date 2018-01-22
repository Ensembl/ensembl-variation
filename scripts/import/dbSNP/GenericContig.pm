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




=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


use strict;
use warnings;

#generic object for the dbSNP data. Contains the general methods to dump the data into the new Variation database. Any change in the methods
# will need to overload the correspondent method in the subclass for the specie

package dbSNP::GenericContig;

use POSIX;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use ImportUtils qw(dumpSQL debug create_and_load load loadfile get_create_statement);
use dbSNP::ImportTask;
use Bio::EnsEMBL::Variation::Utils::dbSNP qw(get_alleles_from_pattern);

use Progress;
use DBI qw(:sql_types);
use Fcntl qw( LOCK_SH LOCK_EX );
use List::Util qw ( min max );

#our $FARM_BINARY = "bsub -R 'select[gpfs]'  ";  ##EBI
our $FARM_BINARY = "bsub";

our %FARM_PARAMS = (
  
  'default' => {
    'script' => 'run_task.pl',
    'max_concurrent_jobs' => 5,
    'memory' => 8000,
    'queue' => 1,
    'wait_queue' => 0
  },
  
  'allele_table' => {
    'script' => 'run_task.pl',
    'max_concurrent_jobs' => 30, ## reduced from 40 using same mysql server for both 
    'memory' => 2000,
    'queue' => 1,
    'wait_queue' => 0
  },
  
  'allele_table_load' => {
    'script' => 'run_task.pl',
    'max_concurrent_jobs' => 5,
    'memory' => 6000, #was 8000 for human
    'queue' => 2,
    'wait_queue' => 0
  },
  ## individual_genotypes human 137:
  ## ran most of these max_concurrent_jobs=50  & memory=2000
  ## 58 failed =>rerun with memory=8000
  ## 3  failed =>rerun with memory=10000
  ## mouse - 20/29 failed 4G limit - 5 needed 9.5G
  'individual_genotypes' => {
    'script' => 'run_task.pl',
    'max_concurrent_jobs' => 30,
    'memory' => 4000,           
    'queue' => 1,
    'wait_queue' => 0
  },
  
  'individual_genotypes_load' => {
    'script' => 'run_task.pl',
    'max_concurrent_jobs' => 5,
     'memory' =>4000, #was 9000 for human
    'queue' => 2,
    'wait_queue' => 0
  },
  
);

#The queues in ascending run time limit order
our @FARM_QUEUES = (
  'small',
  'normal',
  'long',
  'basement'
);

#The maximum amount of farm memory that we will request
our $MAX_FARM_MEMORY = 15900;
our $FARM_MEMORY_INCREMENT = 4000;

#creates the object and assign the attributes to it (connections, basically)
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbm, $tmp_dir, $tmp_file, $limit, $mapping_file_dir, $dbSNP_VERSION, $ASSEMBLY_VERSION, $GROUP_TERM,$GROUP_LABEL, $skip_routines, $scriptdir, $logh, $source_engine, $schema_name) =
        rearrange([qw(DBMANAGER TMPDIR TMPFILE LIMIT MAPPING_FILE_DIR DBSNP_VERSION ASSEMBLY_VERSION GROUP_TERM GROUP_LABEL SKIP_ROUTINES SCRIPTDIR LOG SOURCE_ENGINE SCHEMA_NAME )],@_);


  my $dbSNP = $dbm->dbSNP()->dbc();
  my $dbCore = $dbm->dbCore()->dbc();
  my $dbVar = $dbm->dbVar()->dbc();
  my $dbInt = $dbm->dbInt()->dbc();
  my $snp_dbname = $dbSNP->dbname();
  my $species = $dbm->species();
  my $shared_db = $dbm->dbSNP_shared();
  my $registry_file = $dbm->registryfile();
  $shared_db .="." if defined $source_engine && $source_engine =~/mssql/ ;

  debug(localtime() . "\tThe shared database is $shared_db");

  return bless {'dbSNP'  => $dbSNP,
		'dbCore' => $dbCore,
		'dbInt'  => $dbInt,
		'dbVar'  => $dbVar, ##this is a dbconnection
		'snp_dbname' => $snp_dbname,
		'species' => $species,
		'tmpdir' => $tmp_dir,
		'tmpfile' => $tmp_file,
		'limit' => $limit,
                'mapping_file_dir' => $mapping_file_dir,
		'dbSNP_version' => $dbSNP_VERSION,
		'dbSNP_share_db' => $shared_db,
		'assembly_version' => $ASSEMBLY_VERSION,
		'skip_routines' => $skip_routines,
		'log' => $logh,
		'registry_file' => $registry_file,
		'scriptdir' => $scriptdir,
		'dbm' => $dbm,
		'group_term'  => $GROUP_TERM,
		'group_label' => $GROUP_LABEL,
		'source_engine' => $source_engine,
		'schema_name'   => $schema_name,

		}, $class;
}

#main and only function in the object that dumps all dbSNP data 
sub dump_dbSNP{ 

  my $self = shift;

  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  #the following steps need to be run when initial starting the job. If job failed for some reason and some steps below are already finished, then can comment them out
  my @subroutines = (
   'create_coredb',
   'source_table', 
   'population_table',
   'individual_table',

   'variation_table',
   'subsnp_synonyms',
   'archive_rs_synonyms',
   'dbSNP_annotations',
   'pubmed_citations',

   'parallelized_individual_genotypes',
   'population_genotypes',
   'parallelized_allele_table',
#   'flanking_sequence_table',
   'allele_string',
   'variation_feature',   
   'cleanup'
  );

  ## add refSeq HGVS as synonym for human only
  push @subroutines, 'import_hgvs' if  $self->{'dbm'}->dbCore()->species eq 'human';
 
  #The GenericContig object has an array where routines that should be skipped can be specified. For now, add create_coredb and cleanup by default
  push(@{$self->{'skip_routines'}},('create_coredb',
				    'cleanup'
       ));
  
  # When resuming after a crash, put already finished modules into this array 
  push(@{$self->{'skip_routines'}},());
  
  my $resume;
  
  my $clock = Progress->new();
  
  # Loop over the subroutines and run each one
  foreach my $subroutine (@subroutines) {
    
    #Check if this subroutine should be skipped
    if (grep($_ eq $subroutine,@{$self->{'skip_routines'}})) {
      debug(localtime() . "\tSkipping $subroutine");
      print $logh localtime() . "\tSkipping $subroutine\n";
      next;
    }
    
    $clock->checkpoint($subroutine);
    print $logh $clock->to_string($subroutine);
    $self->$subroutine($resume);
    print $logh $clock->duration();
  }

}

sub run_on_farm {
  my $self = shift;
  my $jobname = shift;
  my $file_prefix = shift;
  my $task = shift;
  my $task_manager_file = shift;
  my $start = shift;
  my $end = shift;
  my @args = @_;
  
  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};

  my $param_key = $jobname;
  if (!exists($FARM_PARAMS{$jobname})) {
    warn("No farm resource parameters defined for job $jobname, will use default parameters");
    $param_key = 'default';
  }
  
  my $script = $FARM_PARAMS{$param_key}->{'script'};
  my $max_jobs = $FARM_PARAMS{$param_key}->{'max_concurrent_jobs'};
  my $memory = $FARM_PARAMS{$param_key}->{'memory'};
  my $queue = $FARM_QUEUES[$FARM_PARAMS{$param_key}->{'queue'}];
  my $wait_queue = $FARM_QUEUES[$FARM_PARAMS{$param_key}->{'wait_queue'}];
  my $memory_long = $memory . '000';
  my $logfile_prefix = $file_prefix . "_" . $jobname . "-";
  my $array_options = "";
  #If just the start was defined, it should be a list of subtasks that should be re-run
  if (defined($start) && defined($end)) {
    if ($start < $end) {
      $array_options = qq{[$start-$end]%$max_jobs};
    }
    else {
      $array_options = qq{[$start]%$max_jobs};
    }
  }
  elsif (defined($start)) {
    $array_options = "[" . join(",",@{$start}) . qq{]%$max_jobs}; 
  }
  
  #Wrap the command to be executed into a shell script
  my $tempfile = $task . "_" . $self->{'tmpfile'};
  my $task_command = qq{perl $self->{'scriptdir'}/$script -species $self->{'species'} -dbSNP_shared $self->{'dbSNP_share_db'} -registry_file $self->{'registry_file'} -task $task -file_prefix $file_prefix -task_management_file $task_manager_file -tempdir $self->{'tmpdir'} -tempfile $tempfile  -source_engine $self->{source_engine} -schema_name $self->{schema_name}} . join(" ",@args);
  warn"Sending task_command $task_command\n";
  my $script_wrapper = $file_prefix . '_command.sh';
  open(CMD,'>',$script_wrapper);
  flock(CMD,LOCK_EX);
  print CMD qq{#!/usr/local/bin/bash\n};
  print CMD qq{$task_command\n};
  close(CMD);
  print $logh "Running $script with task management file: $task_manager_file\n";
  print $logh Progress::location();
  ### changed $memory_long to $memory for farm3
#  my $bsub_cmd = qq{$FARM_BINARY -R'select[mem>$memory\] rusage[mem=$memory\]' -M$memory_long -q $queue -J'$jobname$array_options' -o $logfile_prefix\%J.%I.out -e $logfile_prefix\%J.%I.err bash $script_wrapper};
my $bsub_cmd = qq{$FARM_BINARY -R'select[mem>$memory\] rusage[mem=$memory\]' -M$memory  -J'$jobname$array_options' -o $logfile_prefix\%J.%I.out -e $logfile_prefix\%J.%I.err bash $script_wrapper}; ### EBI
  warn "sending to farm:  $bsub_cmd\n";
  #Submit the job array to the farm
  my $submission = `$bsub_cmd`;
  my ($jobid) = $submission =~ m/^Job \<([0-9]+)\>/i;
  print $logh Progress::location();
  warn "Need to wait for job id $jobid\nSenfin:$FARM_BINARY -J $jobid\_waiting  -w'ended($jobid)' -K -o $file_prefix\_waiting.out sleep 1\n\n";
  #Submit a job that depends on the job array so that the script will halt  (added meme req after problems on farm3)
  system(qq{$FARM_BINARY -R"select[mem>2000] rusage[mem=2000]" -M2000 -J $jobid\_waiting  -w'ended($jobid)' -K -o $file_prefix\_waiting.out sleep 1});  ###EBI
#  system(qq{$FARM_BINARY -R"select[mem>2000] rusage[mem=2000]" -M2000 -J $jobid\_waiting -q $wait_queue -w'ended($jobid)' -K -o $file_prefix\_waiting.out sleep 1});
  print $logh Progress::location();
  
  #Check the error and output logs for each subtask. If the error file is empty, delete it. If not, warn that the task generated errors. If the output file doesn't say that it completed successfully, report the job as unseccessful and report which tasks that failed
  my $all_successful = 1;
  sleep(60);  ## checking files before they are written on farm3
  my %job_details;
  for (my $index = $start; $index <= $end; $index++) {
    my $outfile = "$logfile_prefix$jobid\.$index\.out";
    my $errfile = "$logfile_prefix$jobid\.$index\.err";
    
    # Is error file empty?
    if (-z $errfile) {
      unlink($errfile);
    }
    else {
      $job_details{$index}->{'generated_error'} = 1;
    }
    
    #Does the outfile say that the process exited successfully? Faster than querying bhist... (?)
    $job_details{$index}->{'success'} = 0;
    open(FH,'<',$outfile)||die "Failed to open $outfile : $!\n";;
    flock(FH,LOCK_SH);
    my $content = "";
    while (<FH>) {
      $content .= "$_ ";
    }
    close(FH);
    
    if ($content =~ m/Successfully completed/i) {
      $job_details{$index}->{'success'} = 1;
    }
    # Else, did we hit the memory limit?
    elsif ($content =~ m/TERM_MEMLIMIT/) {
      $job_details{$index}->{'fail_reason'} = 'OUT_OF_MEMORY';
    }
    # Else, time limit for the queue?
    elsif ($content =~ m/TERM_RUNLIMIT/) {
      $job_details{$index}->{'fail_reason'} = 'OUT_OF_TIME';
    }
    # Else, we don't know why it failed
    else {
      $job_details{$index}->{'fail_reason'} = 'UNKNOWN';
    }
  }
  print $logh Progress::location();
  
  my $message = "";
  if ((my $count = grep(!$job_details{$_}->{'success'},keys(%job_details)))) {
    $all_successful = 0;
    $message = qq{$count subtasks failed. You should re-run them before proceeding!\n};
  }
  if ((my $count = grep($job_details{$_}->{'generated_error'},keys(%job_details)))) {
    $message .= qq{$count subtasks generated error messages, please check the logfiles!\n};
  }
  
  my $result = {
    'success' => $all_successful,
    'jobid' => $jobid,
    'subtask_details' => \%job_details,
    'message' => $message
  };
  
  return $result;
}

sub rerun_farm_job {
  my $self = shift;
  my $iteration = shift;
  my $jobname = shift;
  my $file_prefix = shift;
  my @args = @_;
  
  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  my $param_key = $jobname;
  if (!exists($FARM_PARAMS{$jobname})) {
    warn("No farm resource parameters defined for job $jobname, will use default parameters");
    $param_key = 'default';
  }
  
  my $max_time = 0;
  my $max_memory = 0;
  
  #This is a preparation step for re-running a farm job. What we do is to re-submit to a longer queue (if one is available) and request more memory
  my $current_queue = $FARM_QUEUES[$FARM_PARAMS{$param_key}->{'queue'}];
  if (scalar(@FARM_QUEUES) > ($current_queue + 1)) {
    $FARM_PARAMS{$param_key}->{'queue'}++;
  }
  else {
    $max_time = 1;
  }
  
  # Increase the memory requirements in steps of $FARM_MEMORY_INCREMENT unless we're at the limit already
  if ($FARM_PARAMS{$param_key}->{'memory'} < $MAX_FARM_MEMORY) {
    $FARM_PARAMS{$param_key}->{'memory'} = min($MAX_FARM_MEMORY,$FARM_MEMORY_INCREMENT + $FARM_PARAMS{$param_key}->{'memory'});
  }
  else {
    $max_memory = 1;
  }
  
  #If we have already maxed out the resources, it won't help running the job again
  return undef if ($max_time && $max_memory);
  
  return $self->run_on_farm($jobname,$file_prefix . "_submission",$iteration,@args);
}

sub create_coredb {

  my $self = shift;
  
  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  my $coredb_name = $self->{'dbCore'}->dbname();
  $self->{'dbVar'}->do(qq{CREATE DATABASE $coredb_name});
  print $logh Progress::location();
  debug(localtime() . "\tmake sure create $coredb_name.coord_system");
  my $csid_ref = $self->{'dbCore'}->selectall_arrayref(qq{SELECT coord_system_id from coord_system WHERE name = 'chromosome' and attrib = 'default_version'});
  my $csid;
  if ($csid_ref->[0][0]) {
    $csid = $csid_ref->[0][0];
  }
  $self->{'dbVar'}->do(qq{create table coord_system(coord_system_id int(10) unsigned,species_id int(10) unsigned default 1)});
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{insert into coord_system(coord_system_id) values($csid)});
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{RENAME TABLE coord_system TO $coredb_name.coord_system});
  print $logh Progress::location();
}

sub source_table {
    my $self = shift;
    my $source_name = shift;
   
    #get the version of the dbSNP release
    my $dbSNP_db_name = $self->{'snp_dbname'};
    $dbSNP_db_name =~ s/dbSNP\_//; ## remove extra pre-fix on mysql dbs
    my ($species,$tax_id,$version) = $dbSNP_db_name  =~ m/^(.+)?\_([0-9]+)\_([0-9]+)$/;
    my $url = 'http://www.ncbi.nlm.nih.gov/projects/SNP/';

    $version =  $self->{'dbSNP_version'} unless defined $version;
    $version =~ s/b//;
    if(defined $source_name &&  $source_name=~ /Archive/){ 

        $self->{'dbVar'}->do(qq{INSERT IGNORE INTO source (source_id,name,version,description,url,somatic_status, data_types) VALUES (2, "$source_name",$version,"Former variants names imported from dbSNP", "$url", "mixed","variation_synonym")});
    }    
    else{

        my $dbname = 'dbSNP';

        $self->{'dbVar'}->do(qq{INSERT INTO source (source_id,name,version,description,url,somatic_status, data_types) VALUES (1,"$dbname",$version,"Variants (including SNPs and indels) imported from dbSNP", "$url", "mixed","variation")});
    }
    my $source_id = $self->{'dbVar'}->db_handle->last_insert_id(undef, undef, qw(source source_id))|| die "no insert id for source\n";
    return $source_id;
}

sub table_exists_and_populated {

    # check if a table is present in dbSNP, and that it has some data in it
    
    my $self = shift;
    my $table = shift;

    my $sql;
    if($self->{source_engine} =~/mssql|sqlserver/ ){
	$sql = qq[SELECT OBJECT_ID('$table')];
    }
    elsif($self->{source_engine} =~/pg|postgreSQL/i ){
        $table = "\L$table";
        $sql = qq[SELECT * FROM pg_catalog.pg_tables where tablename like '$table'];
    }
    else{
	$sql = qq[show tables like '$table'];
    }
    my ($obj_id) = $self->{'dbSNP'}->db_handle->selectrow_array($sql);
 
    if (defined $obj_id) {

        # the table exists

        my ($count) = $self->{'dbSNP'}->db_handle->selectrow_array(qq{SELECT COUNT(*) FROM $table});

        if (defined $count && $count > 0) {
            
            # the table is populated

            debug(localtime() . "\t$table table exists and is populated");

            return 1;
        }
    }

    debug(localtime() . "\t$table table either doesn't exist or is empty");

    return 0;
}


sub minor_allele_freq {
    my $self = shift;

    return unless $self->table_exists_and_populated('SNPAlleleFreq_TGP');

    my $logh = $self->{'log'};

    # dump the data from dbSNP and store in the temporary maf table in
    # the variation database  [2hrs to dump 137]

    print $logh Progress::location();
    debug(localtime() . "\tDumping global minor allele freqs");

    my $shared = $self->{'dbSNP_share_db'};

    my $stmt = qq{
        SELECT  af.snp_id, a.allele, af.freq, af.count, af.is_minor_allele
        FROM    SNPAlleleFreq_TGP af, $shared\.Allele a
        WHERE   af.allele_id = a.allele_id
        AND     af.freq <= 0.5 
        AND     af.is_minor_allele = 1
    };

    dumpSQL($self->{'dbSNP'}, $stmt, $self->{source_engine});
    
    print $logh Progress::location();
    debug(localtime() . "\tLoading global minor allele freqs");
    
    create_and_load($self->{'dbVar'}, "maf", "snp_id i* not_null", "allele l", "freq f", "count i", "is_minor_allele i");

    debug(localtime() . "\tComplete MAF update");
    # we don't delete the maf temporary table because we need it for post-processing MAFs = 0.5
}

sub suspect_snps {
    my $self = shift;

    return unless $self->table_exists_and_populated('SNPSuspect');

    my $logh = $self->{'log'};
    
    # dump the data into a temporary suspect table

    print $logh Progress::location();
    debug(localtime() . "\tDumping suspect SNPs");
   
    # create a table to store this data in the variation database,
    # we use a mysql SET column which means the dbSNP values which
    # are stored essentially as a bitfield will automatically be 
    # assigned the correct reasons. 
    #
    # XXX: If dbSNP change anything though this CREATE statement 
    # will need to be changed accordingly.
    #

    my $var_table_sql = qq{
        CREATE TABLE suspect (
            snp_id      INTEGER(10) NOT NULL DEFAULT 0,
            reason_code SET('Paralog','byEST','oldAlign','Para_EST','1kg_failed','','','','','','other') DEFAULT NULL,
            PRIMARY KEY (snp_id)
        )
    };

    $self->{'dbVar'}->do($var_table_sql);

    my $stmt = qq{
        SELECT  ss.snp_id, ss.reason_code
        FROM    SNPSuspect ss
    };

    dumpSQL($self->{'dbSNP'}, $stmt, $self->{source_engine});
    
    print $logh Progress::location();
    debug(localtime() . "\tLoading suspect SNPs");
    
    load($self->{'dbVar'}, "suspect", "snp_id", "reason_code");

    # fail the variations tagged as suspect by adding them to the failed_variation table

    # for the moment we just fail these all for the same reason (using the same 
    # failed_description_id), we don't group them by dbSNP's reason_code, but we 
    # keep this information in the suspect table in case we want to do this at 
    # some point in the future

    print $logh Progress::location();
    debug(localtime() . "\tFailing suspect variations");

    my $stmt = qq{
        INSERT IGNORE INTO failed_variation (variation_id, failed_description_id)
        SELECT  v.variation_id, fd.failed_description_id
        FROM    suspect s, variation v, failed_description fd
        WHERE   fd.description = 'Flagged as suspect by dbSNP'
        AND     v.snp_id =  s.snp_id
    };
    my $fail_rs_ins_sth = $self->{'dbVar'}->prepare($stmt);
    $self->{'dbVar'}->do($stmt);

    # also fail any variants with synonyms, use INSERT IGNORE because dbSNP
    # redundantly fail both the refsnp and the synonym sometimes, and we have
    # a unique constraint on (variation_id, failed_description_id)
    debug(localtime() . "\tFailing suspect variations by synonym");
    
     $stmt = qq{
        INSERT IGNORE INTO failed_variation (variation_id, failed_description_id)
        SELECT  vs.variation_id, fd.failed_description_id
        FROM    suspect s, variation_synonym vs, failed_description fd
        WHERE   fd.description = 'Flagged as suspect by dbSNP'
        AND     vs.name = CONCAT('rs', s.snp_id)
    };

    my $fail_old_rs_ins_sth = $self->{'dbVar'}->prepare($stmt);
    $fail_rs_ins_sth->execute()||die;
    $fail_old_rs_ins_sth->execute()||die;

    $self->{'dbVar'}->do(qq{DROP TABLE suspect});
}

=head

Flag named elements (like (Z6867)) in zebrafish)
- not enough information provided for later QC or annotation

=cut
sub named_variants {

    my $self = shift;

    my $logh = $self->{'log'};
    
    print $logh Progress::location();
    debug(localtime() . "\tExtracting named variants");
  
    ## check if named variants present for the release
    my $named_ext_stmt = qq[ select SNP.snp_id
                             from SNP, $self->{'dbSNP_share_db'}.UniVariation uv, $self->{'dbSNP_share_db'}.SnpClassCode scc
                             where scc.abbrev ='Named'
                             and uv.univar_id = SNP.univar_id
                             and scc.code = uv.subsnp_class
                           ];

    my $named_ext_sth = $self->{'dbSNP'}->prepare($named_ext_stmt);
    $named_ext_sth->execute();

    my $named_rs_ids = $named_ext_sth->fetchall_arrayref();

    return unless defined $named_rs_ids->[0]->[0]; 

    ## get attrib id for sequence alteration
    my $attrib_ext_stmt = qq[ select attrib_id
                             from attrib
                             where value ='sequence_alteration'
                           ];

    my $attrib_ext_sth = $self->{'dbVar'}->prepare($attrib_ext_stmt);
    $attrib_ext_sth->execute() || die;
    my $attrib_id = $attrib_ext_sth->fetchall_arrayref();

    die "attribs to be loaded\n\n" unless defined $attrib_id->[0]->[0]; 


    ## update 
    my $var_upd_stmt = qq[ update variation
                           set class_attrib_id = ?
                           where snp_id = ?
                         ];

    my $var_upd_sth = $self->{'dbVar'}->prepare($var_upd_stmt);
    debug(localtime() . "\tUpdating named variants");
    foreach my $rs_id(@{$named_rs_ids}){

       $var_upd_sth->execute( $attrib_id->[0]->[0], $rs_id->[0]) || die;
    } 

}

## extract pubmed ids linked to refsnps
sub pubmed_citations{

    my $self = shift;
   
    return unless $self->table_exists_and_populated('SNPPubmed');
    my $logh = $self->{'log'};
       
    print $logh Progress::location();
    debug(localtime() . "\tExporting pubmed cited SNPs");
    
    ## create tmp table & populate with rs ids & pmids
    $self->{'dbVar'}->do(qq[ create table tmp_pubmed (
                             snp_id varchar(255) not null,
                             pubmed_id int(10) unsigned not null,
                             key snp_idx (snp_id )
                   )]);

    my $pubmed_ins_sth = $self->{'dbVar'}->prepare(qq[ insert into tmp_pubmed (snp_id, pubmed_id) values (?,?)]);

    my $sth = $self->{'dbSNP'}->prepare(qq[ SELECT snp_id, pubmed_id from SNPPubmed where pubmed_id is not null ]);
    $sth->execute();

    print $logh Progress::location();
    while (my $l = $sth->fetchrow_arrayref()){
       $pubmed_ins_sth->execute($l->[0], $l->[1]) ||die "Failed to enter pubmed_id for rs: $l->[0], PMID:$l->[1] \n";
    }
    ## move data to correct structure
    my $pubmed_ext_sth = $self->{'dbVar'}->prepare(qq [ select variation.variation_id, tmp_pubmed.pubmed_id
                                                        from variation,tmp_pubmed
                                                        where variation.snp_id = tmp_pubmed.snp_id]);

    ## linked to PMID via study
    my $publication_ext_sth    = $self->{'dbVar'}->prepare(qq[select publication_id from publication where pmid = ?]);

    my $publication_ins_sth    = $self->{'dbVar'}->prepare(qq[insert into publication (pmid ) values (?)]);

    my $citation_ins_sth = $self->{'dbVar'}->prepare(qq[insert into variation_citation ( variation_id, publication_id ) values (?,?)]);


    $pubmed_ext_sth->execute()||die "Failed to extract pubmed data from tmp table\n";;
    my $data2 = $pubmed_ext_sth->fetchall_arrayref();

    my %done;

    foreach my $l (@{$data2}){

	next if $done{$l->[0]}{$l->[1]};
	$done{$l->[0]}{$l->[1]} = 1;

	$publication_ext_sth->execute( $l->[1] )||die;
	my $publication_id = $publication_ext_sth->fetchall_arrayref();

	unless ( defined $publication_id->[0]->[0]){
            ## create one study per PMID
	    $publication_ins_sth->execute( $l->[1]);
	    $publication_ext_sth->execute( $l->[1] )||die;
	    $publication_id = $publication_ext_sth->fetchall_arrayref();
	}

	warn "problem adding new publication $l->[1]\n" unless defined $publication_id->[0]->[0];

	$citation_ins_sth->execute($l->[0], $publication_id->[0]->[0])||die;

    }
    ## remove tmp table
    $self->{'dbVar'}->do(qq [ drop table tmp_pubmed ]);
    debug(localtime() . "\tCompleted pubmed cited SNPs");


}


# filling of the variation table from SubSNP and SNP
# creating of a link table variation_id --> subsnp_id
sub variation_table {
    my $self = shift;
  
  #If this variable is set, variations with a subsnp_id below this integer will not be imported. This is useful for resuming when resuming an import that crashed at a particular subsnp_id. Also, any SQL statements preparing for the import will be skipped.
  my $resume_at_subsnp_id = -1;
  
  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  
    my $stmt;

  print $logh Progress::location();


   debug(localtime() . "\tDumping RefSNPs");
    

     $stmt = "SELECT ";
     if ($self->{'limit'}) {
       $stmt .= "TOP $self->{'limit'} ";
     }
    if( $self->{source_engine} =~/mssql/ ){
	$stmt .= qq{ 1, 
	             'rs'+LTRIM(STR(snp.snp_id)) AS sorting_id, 
	            snp.snp_id
		    FROM SNP snp 
		    WHERE exemplar_subsnp_id != 0
	       } ;
    }
    elsif( $self->{source_engine} =~/postgreSQL/i ){
	$stmt .= qq{ 1, 
	             'rs' || snp.snp_id AS sorting_id, 
	            snp.snp_id
		    FROM SNP snp 
		    WHERE exemplar_subsnp_id != 0
	       } ;
    }
    else{
	$stmt .= qq{ 1, 
	          CONCAT('rs', CAST(snp.snp_id AS CHAR)) AS sorting_id, 	          
	          snp.snp_id
		  FROM SNP snp 
		  WHERE exemplar_subsnp_id != 0
	       };
    }
    
    if ($self->{'limit'}) {
	$stmt .= qq{    
		   ORDER BY
		     sorting_id ASC  
	          };
    }
    dumpSQL($self->{'dbSNP'},$stmt, $self->{source_engine}) unless ($resume_at_subsnp_id > 0);

    
    debug(localtime() . "\tLoading RefSNPs into variation table");
    ### time test this
    unless ($resume_at_subsnp_id > 0){
      $self->{'dbVar'}->do( qq[ALTER TABLE variation add column snp_id int NOT NULL] ) ;
  
      $self->{'dbVar'}->do( qq[ALTER TABLE variation disable keys]);  

    }
    load( $self->{'dbVar'}, "variation", "source_id", "name", "snp_id" ) unless ($resume_at_subsnp_id > 0);
    $self->{'dbVar'}->do( "ALTER TABLE variation ADD INDEX snpidx( snp_id )" ) unless ($resume_at_subsnp_id > 0);
    debug(localtime() . "\tVariation table loaded");
    
    
    $self->{'dbVar'}->do( qq[ALTER TABLE variation enable keys] );
    
    
    debug(localtime() . "\tVariation table indexed");


    debug(localtime() . "\tStarting subSNPhandle");
    #create a subsnp_handle table
    $stmt = qq{
               SELECT 
                 s.subsnp_id, 
                 b.handle
	       FROM 
	         SubSNP s, 
	         Batch b
               WHERE 
                 s.batch_id = 
                 b.batch_id
              };
    dumpSQL($self->{'dbSNP'},$stmt, $self->{source_engine});
    load( $self->{'dbVar'}, "subsnp_handle", "subsnp_id", "handle");
    debug(localtime() . "\tFinished subSNPhandle");
   


    return unless  $self->{'dbm'}->dbCore()->species =~ /homo/i ;

    debug(localtime() . "\tLooking for Somatic variants");
    #Get somatic flag from dbSNP but only if the snp exclusively has somatic subsnps
    $stmt = qq{ SELECT sssl.snp_id
                FROM SubSNP ss JOIN
                SNPSubSNPLink sssl ON ( sssl.subsnp_id = ss.subsnp_id AND ss.SOMATIC_ind = 'Y')
                where not exists (select * from SubSNP ss2 JOIN SNPSubSNPLink sssl2 ON ( sssl2.subsnp_id = ss2.subsnp_id AND ss2.SOMATIC_ind = 'N')
                where  sssl2.snp_id = sssl.snp_id )
	};
	my $sth = $self->{'dbSNP'}->prepare($stmt);
  	print $logh Progress::location();
        $sth->execute();
  	print $logh Progress::location();
    
    my $snp_id;
    $sth->bind_columns(\$snp_id);

    
    # Loop over the somatic SNPs and set the flag in the variation table
    $stmt = qq{
    	UPDATE
    		variation
    	SET
    		somatic = 1
    	WHERE
    		snp_id = ?
    };
    my $up_sth = $self->{'dbVar'}->prepare($stmt);
    while ($sth->fetch()) {
    	$up_sth->execute($snp_id);
    }
    $sth->finish();
    $up_sth->finish();
    print $logh Progress::location();
    debug(localtime() . "\tFinished Somatic variants");


    return;
}

# import any clinical significance, global minor allele frequencies or suspect SNPs
# these subroutines all check if the table is present and populated before doing
# anything so should work fine on species without the necessary tables

sub dbSNP_annotations{

    my $self = shift;

    ## This is no longer run - ClinVar export has more information
    #debug(localtime() . "\tStarting clin_sig");
    #$self->clin_sig;
    if($self->{'dbm'}->dbCore()->species =~ /homo/i){
      debug(localtime() . "\tStarting MAF");
      $self->minor_allele_freq;
      debug(localtime() . "\tStarting suspect SNP");
      $self->suspect_snps;   
    }
    debug(localtime() . "\tStarting named variants");
    $self->named_variants;
    debug(localtime() . "\tFinished variation table with dbSNP annotations");
    return;
}

# extract ss id & strand relative to rs
# this mapping is used extensively in the import process 
# but discarded for human when import id complete due to Mart problems

sub subsnp_synonyms{

    my $self = shift;
    my $logh = $self->{'log'};

   # create a temp table of subSNP info
    
    debug(localtime() . "\tDumping SubSNPs");    

   my $stmt = "SELECT ";
   if ($self->{'limit'}) {
     $stmt .= "TOP $self->{'limit'} ";
   }
   $stmt .= qq{
                 subsnp.subsnp_id , 
                 subsnplink.snp_id, 
	         subsnplink.substrand_reversed_flag,
                 b.handle
	       FROM 
	         SubSNP subsnp, 
	         SNPSubSNPLink subsnplink, 
	         Batch b
	       WHERE 
	         subsnp.batch_id = b.batch_id AND   
	         subsnp.subsnp_id = subsnplink.subsnp_id
	      };
    if ($self->{'limit'}) {
      $stmt .= qq{    
	          ORDER BY
	            subsnp_id ASC  
	         };
    }
   dumpSQL($self->{'dbSNP'},$stmt,  $self->{source_engine} ) ; 
   create_and_load( $self->{'dbVar'}, "tmp_var_allele", "subsnp_id i*  not_null unsigned", "refsnp_id v* not_null", "substrand_reversed_flag i", "submitter_handle v* not_null");
  
  print $logh Progress::location(); 


  # load the synonym table with the subsnp identifiers
    
   debug(localtime() . "\tloading variation_synonym table with subsnps");


  print $logh Progress::location();
   
  # Subsnp_id interval to export at each go
  my $interval = 5e5;
    my $offset = 0;
  #Get the minimum and maximum subsnp_id. We will export rows based on subsnp_id rather than using limit and offset in order to avoid having the same subsnp_id split across two exports
  $stmt = qq{
    SELECT
      MIN(tv.subsnp_id) AS mn,
      MAX(tv.subsnp_id) AS mx
    FROM
      tmp_var_allele tv
  };

  


  my ($min_id,$max_id) = @{$self->{'dbVar'}->db_handle()->selectall_arrayref($stmt)->[0]};
  print $logh Progress::location();

  debug(localtime() . "\tDoing variation_synonym insert");
  ## remove indexes for quicker loading
  $self->{'dbVar'}->do("ALTER TABLE variation_synonym disable keys ");

  $self->{'dbVar'}->do(qq{ALTER TABLE variation_synonym add column substrand_reversed_flag tinyint});
  $self->{'dbVar'}->do(qq{ALTER TABLE variation_synonym add column submitter_handle varchar(25)});

  while ($offset < $max_id) {

    my $end = $offset + $interval;

    $self->{'dbVar'}->do( qq{ insert into variation_synonym (variation_id, subsnp_id, source_id, name, substrand_reversed_flag, submitter_handle)
               (SELECT  distinct v.variation_id, tv.subsnp_id, 1, CONCAT( 'ss', tv.subsnp_id), tv.substrand_reversed_flag, tv.submitter_handle
               FROM tmp_var_allele tv,
                    variation v USE INDEX (snpidx)
               WHERE tv.subsnp_id BETWEEN $offset  AND $end 
               AND v.snp_id = tv.refsnp_id
               ORDER BY   tv.subsnp_id ASC)
               });


 # Increase the offset
    $offset += ($interval + 1);
  }
  debug(localtime() . "\tIndexing variation_synonym table ");
  $self->{'dbVar'}->do("ALTER TABLE variation_synonym enable keys");

  debug(localtime() . "\tCreating subsnp_map table ");
  ## create subsnp_map table

    $self->{'dbVar'}->do(qq[ CREATE TABLE subsnp_map
                            ( variation_id  int(11) unsigned NOT NULL,
                              subsnp_id     int(11) unsigned DEFAULT NULL)
                            engine=MyISAM
                            ]);

    $self->{'dbVar'}->do( qq[ insert into subsnp_map (variation_id, subsnp_id)  select variation_id, subsnp_id from variation_synonym]);

    $self->{'dbVar'}->do(qq[ CREATE INDEX variation_idx on subsnp_map (variation_id) ]);
    print $logh Progress::location();
  
    return;
}

#
# loads the population table
#
# This subroutine produces identical results as the MySQL equivalence
#
sub population_table {
    my $self = shift;
  
  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  $self->{'dbVar'}->do("ALTER TABLE population ADD column pop_id int");   
  $self->{'dbVar'}->do("ALTER TABLE population ADD column pop_class_id int");   ## hold id of super pop in super pop


  print $logh Progress::location();

  # load PopClassCode data as populations 
  #    - these are super_populations  like 'MULTI-NATIONAL' and 'NORTH/EAST AFRICA & MIDDLE EAST'

  debug(localtime() . "\tDumping population class data" );

    my $stmt;
    if($self->{source_engine} =~/postgreSQL/i){
	$stmt = qq{
                SELECT 
                  btrim(pop_class), 
                  (pop_class_id), 
                  btrim(pop_class_text)
		FROM 
		  $self->{'dbSNP_share_db'}.PopClassCode
               };
    }
    else{
	$stmt = qq{
                SELECT 
                  RTRIM(pop_class), 
                  RTRIM(pop_class_id), 
                  RTRIM(pop_class_text)
		FROM 
		  $self->{'dbSNP_share_db'}.PopClassCode
               };
    }
  dumpSQL($self->{'dbSNP'}, $stmt, $self->{source_engine} );

  load($self->{'dbVar'}, 'population', 'name', 'pop_class_id', 'description');
  $self->{'dbVar'}->do(qq{ALTER TABLE population ADD INDEX pop_class_id (pop_class_id)}); 
  print $logh Progress::location();



  debug(localtime() . "\tDumping population data" );

  # load Population data as populations
    my $concat_syntax ;
    if($self->{source_engine} =~/mssql/){
	$concat_syntax = qq[ p.handle+':'+p.loc_pop_id ];
    }
    elsif($self->{source_engine} =~/postgreSQL/i){
	$concat_syntax = qq[ p.handle || ':' || loc_pop_id ];
    }
    else{
	$concat_syntax = qq[ CONCAT(p.handle,':',p.loc_pop_id) ];
    }
    
    $stmt = qq{
            SELECT DISTINCT 
              $concat_syntax,
	      p.pop_id, 
	      pc.pop_class_id, 
	      pl.line,
	      pl.line_num
	    FROM   
	      Population p LEFT JOIN 
	      $self->{'dbSNP_share_db'}.PopClass pc ON p.pop_id = pc.pop_id LEFT JOIN 
	      PopLine pl ON p.pop_id = pl.pop_id
	    ORDER BY
	      p.pop_id ASC,
	      pc.pop_class_id ASC,
	      pl.line_num ASC
            };	 #table size is small, so no need to change
    dumpSQL($self->{'dbSNP'}, $stmt,  $self->{source_engine} );

    debug(localtime() . "\tLoading population data");

    create_and_load( $self->{'dbVar'}, "tmp_pop", "name", "pop_id i*", "pop_class_id i*", "description l", "line_num i*" );
  print $logh Progress::location();
                   
    #populate the population table with the populations
    $self->{'dbVar'}->do("SET SESSION group_concat_max_len = 10000");
    $self->{'dbVar'}->do(qq{INSERT INTO population (name, pop_id,description)
                 SELECT tp.name, tp.pop_id, GROUP_CONCAT(description ORDER BY tp.pop_class_id ASC, tp.line_num ASC)
                 FROM   tmp_pop tp
                 GROUP BY tp.name, tp.pop_id
                 });	#table size is small, so no need to change
  print $logh Progress::location();

    $self->{'dbVar'}->do(qq{ALTER TABLE population ADD INDEX pop_id (pop_id)});
  print $logh Progress::location();

    
     debug(localtime() . "\tLoading population_synonym table");

     # build super/sub population relationships
     $self->{'dbVar'}->do(qq{INSERT INTO population_structure (super_population_id,sub_population_id)
 				    SELECT DISTINCT p1.population_id, p2.population_id
 				    FROM tmp_pop tp, population p1, population p2
 				    WHERE tp.pop_class_id = p1.pop_class_id
 				    AND   tp.pop_id = p2.pop_id});
  print $logh Progress::location();
    

     #load population_synonym table with dbSNP population id
     $self->{'dbVar'}->do(qq{INSERT INTO population_synonym (population_id,source_id,name)
 				      SELECT population_id, 1, pop_id
 				      FROM population
 				      WHERE pop_id is NOT NULL
 				  });

    

  print $logh Progress::location();
    
     $self->{'dbVar'}->do("DROP TABLE tmp_pop");
  print $logh Progress::location();
}



# loads the individual table
#
# pre-e!70 samples were merged on dbSNP ind_id and a submitted name chosen at random
# post-e!70 samples are only merged if they have the same name and ind_id
# post e!85 samples are not merged

# gender and ped info are held on the ind_id

sub individual_table {
    my $self = shift;

  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  debug(localtime() . "\tStarting on Individual data");


 #decide which individual_type should this species be make sure it's correct when adding new species
  my $individual_type_id;
  if ($self->{'dbm'}->dbCore()->species =~ /homo|pan|anoph/i) {
      $individual_type_id = 3;
  }
  elsif ($self->{'dbm'}->dbCore()->species =~ /mus/i) {
      $individual_type_id = 1;
  }
  else {
      $individual_type_id = 2;
  }


   ## create temp table to hold dbSNP submitted_ind_id to simplify genotype loading
   $self->{'dbVar'}->do(qq{
                          CREATE TABLE tmp_ind (
                            submitted_ind_id int(10) unsigned not null,
                            sample_id int(10) unsigned not null,
                            primary key(submitted_ind_id),
                            key individual_id_idx (sample_id))
                         });


    ## get pedigree info - held at dbSNP clustered Individual level 
    my $ped = $self->get_ped_data();

    ## get individual names, submitters and dbSNP clustered Individual info
    my ($individuals, $samples) = $self->get_ind_data();

    ## get population ids for pre-loaded populations
    my $pop_ids = $self->get_pop_ids();
 

    ## prepare insert statements   

    my $ind_ins_sth = $self->{'dbVar'}->prepare(qq[ INSERT INTO individual ( name, description, individual_type_id) values (?,?,?)]);

    my $tmp_ins_sth = $self->{'dbVar'}->prepare(qq[ INSERT INTO tmp_ind (sample_id, submitted_ind_id) values (?,?)]);

    my $pop_ins_sth = $self->{'dbVar'}->prepare(qq[ INSERT INTO sample_population (sample_id, population_id) values (?,?)  ]);

    my $syn_ins_sth = $self->{'dbVar'}->prepare(qq[ INSERT INTO individual_synonym (individual_id,source_id,name)  values (?,?,?)  ]);

    my $sam_ins_sth = $self->{'dbVar'}->prepare(qq[ INSERT INTO sample ( individual_id, name, description) values (?,?,?)]);

    my $ind_gen_upd_sth = $self->{'dbVar'}->prepare(qq[ update individual set gender =? where individual_id = ? ]);
    my $ind_par_upd_sth = $self->{'dbVar'}->prepare(qq[ update individual set father_individual_id =?, mother_individual_id =?  where individual_id = ? ]);



    ## insert sample & individual data
    my %individual_id;

    foreach my $sam_id (keys %$samples){

        ## insert individual if novel - using first sample name seen as name
        ## a few species (eg human) have multiple sample names per individual; check later
        unless (defined $individual_id{ $samples->{$sam_id}->{ind}} ){
            $ind_ins_sth->execute(  $samples->{$sam_id}->{name}  , $individuals->{ $samples->{$sam_id}->{ind} }->{des}, $individual_type_id );
            $individual_id{ $samples->{$sam_id}->{ind} } =  $self->{'dbVar'}->db_handle->last_insert_id(undef, undef, 'individual', 'individual_id');

            ## insert dbSNP synonym 
            $syn_ins_sth->execute( $individual_id{ $samples->{$sam_id}->{ind} }, 1, $samples->{$sam_id}->{ind} );
        }


        ## add sample data 
        $sam_ins_sth->execute( $individual_id{ $samples->{$sam_id}->{ind} }, $samples->{$sam_id}->{name}, $individuals->{ $samples->{$sam_id}->{ind} }->{des});

        ## extract sample_id
        my $ens_sample_id =  $self->{'dbVar'}->db_handle->last_insert_id(undef, undef, 'sample', 'sample_id');

        ## save temp look up table of ensembl id => submitted_ind_id for genotype import
        $tmp_ins_sth->execute($ens_sample_id, $sam_id);

        ## attach the sample to a population if there is one
        $pop_ins_sth->execute( $ens_sample_id , $pop_ids->{ $samples->{$sam_id}{pid} } )
            if defined $samples->{$sam_id}{pid};
   }

   
    ## insert individual gender & ped info

    foreach my $ind (keys %$individuals){

        $ind_gen_upd_sth->execute( $ped->{$ind}{gender}, $individual_id{$ind})
          if defined $ped->{$ind}{gender} && $ped->{$ind}{gender} =~/ale/;;

        my $mother_id ;
        my $father_id ;

        ## get ensmebl individual ids for parents
	if(defined $ped->{$ind}{father} ){
            ## should have and ensembl individual id 
	    if(defined $individual_id{ $ped->{$ind}{father}}){
		$father_id = $individual_id{ $ped->{$ind}{father} };
	    }else{
		warn "No ensembl id found for father : $ped->{$ind}{father} from child id $ind\n";
	    }
	}
	if(defined $ped->{$ind}{mother} ){
	    if(defined $individual_id{ $ped->{$ind}{mother} }){
		$mother_id = $individual_id{ $ped->{$ind}{mother} };
	    }
            else{
		warn "No ensembl id found for mother : $ped->{$ind}{mother} from child id $ind\n";
	    }
	}

	next unless (defined $father_id && $father_id =~ /\d+/ || defined $mother_id && $mother_id =~ /\d+/);
        $ind_par_upd_sth->execute( $father_id,
                                   $mother_id,
                                   $individual_id{$ind}
	                         );
    }
    

    ## update population.size once all individuals loaded
    $self->update_population_size();

    debug(localtime() . "\tIndividual data loaded");
    print $logh Progress::location();

    return;
}

## export family and gender info for individual
## held at level of curated Individual rather than SubmittedIndividual
sub get_ped_data{

    my $self = shift;

    my %ped;

    my $all_ped_sth = $self->{'dbSNP'}->prepare(qq[ SELECT ind_id, 
                                                           pa_ind_id, 
                                                           ma_ind_id, 
                                                           sex
                                                    FROM  PedigreeIndividual
                                                  ])||die "ERROR preparing ss_sth: $DBI::errstr\n";  

    $all_ped_sth->execute()||die "ERROR executing: $DBI::errstr\n";
    my $ped = $all_ped_sth->fetchall_arrayref();

    foreach my $l(@{$ped}){

	if(defined $l->[3]){
	    $l->[3] =~ s/M/Male/;	   
	    $l->[3] =~ s/F/Female/;
	}

	$ped{$l->[0]}{father} = $l->[1];
	$ped{$l->[0]}{mother} = $l->[2];
	$ped{$l->[0]}{gender} = $l->[3];
    }
    return \%ped;

}

## Extract submitted individual data
sub get_ind_data{

   my $self = shift;

   my %samples;
   my %individuals;

   my $all_ind_sth = $self->{'dbSNP'}->prepare(qq[ SELECT si.submitted_ind_id,
                                                          si.loc_ind_id_upp,
                                                          i.descrip, 
                                                          i.ind_id ,
                                                          si.pop_id   
                                                   FROM  SubmittedIndividual si LEFT OUTER JOIN  Individual i on (si.ind_id = i.ind_id)   
                                                   ])||die "ERROR preparing ss_sth: $DBI::errstr\n";


   $all_ind_sth->execute()||die "ERROR executing: $DBI::errstr\n";
   my $inds = $all_ind_sth->fetchall_arrayref();

   foreach my $l(@{$inds}){

       ## individual description is often not useful
       $individuals{$l->[3]}{des}  = $l->[2]
         if defined $l->[2] && $l->[2] !~ /^unknown|byPopDesc/i ;

       ## sample info
       $samples{$l->[0]}{name} = $l->[1];
       $samples{$l->[0]}{ind}  = $l->[3] if defined $l->[3];
       $samples{$l->[0]}{pid}  = $l->[4] if defined $l->[4];

   }
   
   return (\%individuals, \%samples);
}

## look up population_id for dbSNP pop_id to link individuals to pops
sub get_pop_ids{

    my $self = shift;

    my $pop_ext_sth =   $self->{'dbVar'}->prepare(qq[ select population_id, pop_id from population where pop_id is not null ]);

    my %pop_id;
    $pop_ext_sth->execute();
    my $pop_link = $pop_ext_sth->fetchall_arrayref();
    foreach my $l (@{$pop_link}){
	$pop_id{$l->[1]} = $l->[0];
    }

    return \%pop_id;
}


## count and store the number of individuals in a population
sub update_population_size{

    my $self = shift;

    my $size_ext_sth = $self->{'dbVar'}->prepare(qq[ select population.population_id, count(*) 
                                                     from population, sample_population
                                                     where population.population_id = sample_population.population_id
                                                     group by population.population_id
                                                   ]);

    my $size_upd_sth = $self->{'dbVar'}->prepare(qq[ update population
                                                     set size =?
                                                     where population_id = ?
                                                    ]);


    $size_ext_sth->execute()||die "Failed to extact individual counts for populations\n";
    my $sizes =  $size_ext_sth->fetchall_arrayref();
    foreach my $l (@{$sizes}){

	$size_upd_sth->execute( $l->[1], $l->[0])||die "Failed to update individual counts for populations\n";
    }

    return;

}
sub parallelized_allele_table {
  my $self = shift;
  my $load_only = shift;

  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  my $stmt;

  debug(localtime() . "\tStarting allele table");

 ## revert schema to use letter rather than number coded alleles
 update_allele_schema($self->{'dbVar'});  
  
  
  # The tempfile to be used for loading
  my $file_prefix = $self->{'tmpdir'} . '/allele_table';
  my $loadfile = $file_prefix . '_loadfile.txt';

  # Sample file is used for caching the samples
  my $samplefile = $file_prefix . '_population_samples.txt';
  # Allele file is used for caching the alleles
  my $allelefile = $file_prefix . '_alleles.txt';
  
  # Do the extraction of the alleles unless we're resuming and will only load the alleles
  my $task_manager_file = $file_prefix . '_task_management.txt';

  my $jobindex;
  my $jobid;
  my $jobname;
  my $result;
 print $logh Progress::location();

  unless ($load_only) {

##=head hashed out to re-run failed job
    #First, get the population_id -> sample_id mapping and write it to the file. The subroutine can also get it but it's faster to get all at once since we expect many to be used. 
    $stmt = qq{
      SELECT DISTINCT
	p.pop_id,
	p.population_id
      FROM
	population p
      WHERE
	p.pop_id IS NOT NULL AND
	p.pop_id > 0
    };
    my $sth = $self->{'dbVar'}->prepare($stmt);
    $sth->execute();
    my @samples;
    while (my @row = $sth->fetchrow_array()) {
      push(@samples,@row);
    }
    my %s = (@samples);
    dbSNP::ImportTask::write_samples($samplefile,\%s);
    print $logh Progress::location();
  
    #Process the alleles in chunks based on the SubSNP id. This number should be kept at a reasonable level, ideally so that we don't need to request more than 4GB memory on the farm. If so, we'll have access to the most machines.
    #The limitation is that the results from the dbSNP query is read into memory. If necessary, we can dump it to a file (if the results are sorted)


    ## improve binning for sparsely submitted species
    #hash out task file creation if re-running 1 failed job
    if($self->{source_engine} =~/postgreSQL/ ){
        $jobindex =write_allele_task_file_pg($self->{'dbSNP'}->db_handle(),$task_manager_file, $loadfile,  $allelefile, $samplefile,$self->{limit});
    }
    else{
        $jobindex = write_allele_task_file($self->{'dbSNP'}->db_handle(),$task_manager_file, $loadfile,  $allelefile, $samplefile,$self->{limit});
    }
    my $start = 1 ;
    debug(localtime() . "\tAt allele table - written export task file");
    # Run the job on the farm
    $jobname = 'allele_table';
    $result = $self->run_on_farm($jobname,$file_prefix,'allele_table',$task_manager_file,$start,$jobindex); 
    $jobid = $result->{'jobid'};
    debug(localtime() . "\tAt allele table - export jobs run on farm");
    # Check if any subtasks generated errors
    if ((my @error_subtasks = grep($result->{'subtask_details'}{$_}{'generated_error'},keys(%{$result->{'subtask_details'}})))) {
      warn($result->{'message'});
      print $logh Progress::location() . "\t " . $result->{'message'};
      print $logh Progress::location() . "\tThe following subtasks generated errors for job $jobid:\n";
      foreach my $index (@error_subtasks) {
	print $logh qq{\t\t$jobname\[$index\] ($jobid\[$index\])\n};
      }
    }
    # Check the result, if anything failed
    if (!$result->{'success'}) {
      warn($result->{'message'});
      print $logh Progress::location() . "\t " . $result->{'message'};
      print $logh Progress::location() . "\tThe following subtasks failed for unknown reasons for job $jobid:\n";
      foreach my $index (grep($result->{'subtask_details'}{$_}{'fail_reason'} =~ m/UNKNOWN/,keys(%{$result->{'subtask_details'}}))) {
	print $logh qq{\t\t$jobname\[$index\] ($jobid\[$index\])\n};
      }
    
      print $logh Progress::location() . "\tThe following subtasks failed because they ran out of resources for job $jobid:\n";
      foreach my $index (grep($result->{'subtask_details'}{$_}{'fail_reason'} =~ m/OUT_OF_[MEMORY|TIME]/,keys(%{$result->{'subtask_details'}}))) {
	print $logh qq{\t\t$jobname\[$index\] ($jobid\[$index\])\n};
      }
    }
    # If we still have subtasks that fail, this needs to be resolved before proceeding
    die("Some subtasks are failing (see log output). This needs to be resolved before proceeding with the loading of alleles!") unless ($result->{'success'});
  }

# $jobindex =  2606;  ##Put number of subfiles here if running on load_only

  debug(localtime() . " Creating single allele file to load ");

  ## merge all data files into one for loading
  die "Exiting: Number of allele files to merge not known - is load_only set?\n" unless defined $jobindex;
  my $new_load_file_name =  $self->{'tmpdir'} . "/allele_load_file_new.txt";
  open my $new_load_file, ">", $new_load_file_name || die "Failed to open allele load file to write:$!\n";
  foreach my $n(1..$jobindex){
      my $cat_file = $self->{'tmpdir'} ."/allele_table_loadfile.txt_$n";

      open my $subfile, $cat_file ||die  "Failed to open $cat_file to read:$!\n";
      while(<$subfile>){print $new_load_file $_;} 
  }
  close $new_load_file;

  ### start loading process - run as single job

  $self->{'dbVar'}->do( qq[ LOAD DATA LOCAL INFILE "$new_load_file_name" INTO TABLE allele( variation_id,subsnp_id,population_id,allele,frequency,count,frequency_submitter_handle )]) || die "Erro loading allele data: $self->{'dbVar'}::errstr \n";




  ## add indexes post load 
  debug(localtime() . "\tAt allele table - data load complete");
  $self->{'dbVar'}->do( qq[ alter table allele enable keys]);
  debug(localtime() . "\tAt allele table - indexing complete");
  debug(localtime() . "\tEnding allele table");

}

## this can run to extract data from dbSNP while alleles loaded into ensembl
sub allele_string{

  my $self = shift;
  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  my $stmt;


  debug(localtime() . "\tAt allele_string table - starting");

  #Create the allele_string table needed for variation_feature    ## change to store A/T rather than seperate rows per allele
  $stmt = qq{
    SELECT 
      snp.snp_id, 
      uv.var_str
    FROM   
      SNP snp JOIN 
      $self->{'dbSNP_share_db'}.UniVariation uv ON (
	uv.univar_id = snp.univar_id
      )
  };
  dumpSQL($self->{'dbSNP'},$stmt, $self->{source_engine} );
  print $logh Progress::location();
  create_and_load($self->{'dbVar'},"tmp_allele_string","snp_name * not_null","allele");
  print $logh Progress::location();
    
  $stmt = qq{
    CREATE TABLE
      allele_string
    SELECT
      v.variation_id AS variation_id,
      tas.allele AS allele_string
    FROM
      variation v,
      tmp_allele_string tas
    WHERE
      v.name = CONCAT(
	"rs",
	tas.snp_name
      )
  };
  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();

  $stmt = qq{
    ALTER TABLE
      allele_string
    ADD INDEX
      variation_idx (variation_id)
  };

  $self->{'dbVar'}->do($stmt);

  print $logh Progress::location();
  debug(localtime() . "\tEnding allele string table");

}


## revert schema to use letter rather than number coded alleles
sub  update_allele_schema{

  my $dbVar_dbh = shift;

  $dbVar_dbh->do( qq[ drop table allele ] ); 

  $dbVar_dbh->do( qq [ CREATE TABLE allele (
                          allele_id int(10) unsigned NOT NULL AUTO_INCREMENT,
                          variation_id int(10) unsigned NOT NULL,
                          subsnp_id int(15) unsigned NOT NULL,
                          allele varchar(25000) NOT NULL,
                          frequency float DEFAULT NULL,
                          population_id int(10) unsigned DEFAULT NULL,
                          count int(10) unsigned DEFAULT NULL,
                          frequency_submitter_handle int(10) DEFAULT NULL,
                          PRIMARY KEY (allele_id),
                          KEY subsnp_idx (subsnp_id),
                          KEY variation_idx (variation_id))
                          engine=MyISAM
                        ]);  

  ## disable keys for quick loading
  $dbVar_dbh->do( qq[ alter table allele disable keys]);

  return;
}

sub write_allele_task_file{

    my ($dbh, $task_manager_file, $loadfile,  $allelefile, $samplefile, $limit) = @_;
    #warn "Starting allele_task file \n";
    ### previously binning at 500,000, switched to 400,000
    my ($first, $previous,  $jobindex, $ssid);

   open(MGMT,'>',$task_manager_file) || die "Failed to open allele table task management file ($task_manager_file): $!\n";;
    my $stmt = "SELECT ";
    if ($limit) {
	$stmt .= "TOP $limit ";
    }
    $stmt .= qq[ ss.subsnp_id
                 FROM  SubSNP ss
                 order by ss.subsnp_id
		 ];

    my $counter = 0;
    my $ss_extract_sth = $dbh->prepare($stmt);
    $ss_extract_sth->execute() ||die "Error extracting ss ids for allele_table binning\n";

    $ss_extract_sth->bind_columns(\$ssid);

    while( $ss_extract_sth->fetchrow_arrayref()){
	
	$counter++;
	unless(defined $first){
	    ### set up first bin
            $first = $ssid;
#           warn "setting first to $ssid\n";
            $jobindex   = 1;
            next;
        }
        if($counter >=100000){
            ## end bin
             print MGMT qq{$jobindex $loadfile\_$jobindex $first $previous $allelefile $samplefile\n};

            ## start new bin
            $first       = $ssid;

            $counter     = 1;
            $jobindex++;
        }
        else{
            $previous = $ssid;
        }
        
    }


    ### write out last bin
    print MGMT qq{$jobindex $loadfile\_$jobindex $first $previous $allelefile $samplefile\n};
    close MGMT;
    debug(localtime() . "\tFinished allele task file");

    return ($jobindex);
}

sub write_allele_task_file_pg{

    my ($dbh, $task_manager_file, $loadfile,  $allelefile, $samplefile, $limit) = @_;
    #warn "Starting allele_task file \n";
    ### previously binning at 500,000, switched to 400,000
    my ($first, $previous,  $jobindex, $ssid);
    debug(localtime() . "\tAt write_allele_task_file - starting");
   open(MGMT,'>',$task_manager_file) || die "Failed to open allele table task management file ($task_manager_file): $!\n";;
    my $stmt = "SELECT ";
    if ($limit) {
        $stmt .= "TOP $limit ";
    }
    $stmt .= qq[ ss.subsnp_id
                 FROM  SubSNP ss
                 order by ss.subsnp_id
                 ];

    my $counter = 0;
    
    $dbh->begin_work();
    $dbh->do("DECLARE csr CURSOR  FOR $stmt");
    $stmt = "fetch 100000 from csr";
    while (1) {
        my $sth = $dbh->prepare($stmt);
            $sth->execute;
        last if 0 == $sth->rows;
        
        debug(localtime() . "\tAt write_allele_task_file - executed");
        while( $ssid = $sth->fetchrow_arrayref()){
            
            $counter++;
            unless(defined $first){
                ### set up first bin
                $first = $ssid->[0];
#           warn "setting first to $ssid\n";
                $jobindex   = 1;
                next;
            }
            if($counter >=100000){
                ## end bin
                print MGMT qq{$jobindex $loadfile\_$jobindex $first $previous $allelefile $samplefile\n};
                
                ## start new bin
                $first       = $ssid->[0];
                
                $counter     = 1;
                $jobindex++;
            }
            else{
                $previous = $ssid->[0];
            }       
        }
    }
    
    $dbh->do("CLOSE csr");
    $dbh->rollback();

    ### write out last bin
    print MGMT qq{$jobindex $loadfile\_$jobindex $first $previous $allelefile $samplefile\n};
    close MGMT;
    debug(localtime() . "\tFinished allele task file");

    return ($jobindex);
}


#
# loads the flanking sequence table
#
sub flanking_sequence_table {
  my $self = shift;

  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  my $stmt;

  $self->{'dbVar'}->do(qq{CREATE TABLE tmp_seq (variation_id int NOT NULL,
						      subsnp_id int NOT NULL,
						      seq_type int NOT NULL,
						      line_num int,
						      type enum ('5','3'),
						      line varchar(255),
						      revcom tinyint)
				MAX_ROWS = 100000000});

  ## Create flanking sequence table as no longer part of production schema
  $self->{'dbVar'}->do(qq{ create table if not exists flanking_sequence (
                           variation_id int(10) unsigned not null,
                           up_seq text,
 	  	           down_seq text,
 	  	           up_seq_region_start int,
 	  	           up_seq_region_end   int,
 	  	           down_seq_region_start int,
 	  	           down_seq_region_end int,
 	  	           seq_region_id int(10) unsigned,
 	  	           seq_region_strand tinyint,
 	  	           primary key( variation_id )
 	  	           ) MAX_ROWS = 100000000});
  print $logh Progress::location();

# import both the 5prime and 3prime flanking sequence tables
## in human the flanking sequence tables have been partitioned
  if($self->{'dbm'}->dbCore()->species =~ /human|homo/i) {
	foreach my $type ('3','5') {
  
	  foreach my $partition('p1_human','p2_human','p3_human','ins') {
	      warn "Dumping $type flank from $partition\n";
		debug("Dumping $type\_$partition flanking sequence");
  ## TEST THIS
		$stmt = "SELECT ";
		if ($self->{'limit'}) {
		  $stmt .= "TOP $self->{'limit'} ";
		}
		$stmt .= qq{
		    seq.subsnp_id AS sorting_id,
		    seq.type,
		    seq.line_num,
		    seq.line,
        subsnplink.substrand_reversed_flag
		  FROM
		    SubSNPSeq$type\_$partition seq,
		    SNP snp,
        SNPSubSNPLink subsnplink
		  WHERE
		    snp.exemplar_subsnp_id = seq.subsnp_id AND
        seq.subsnp_id = subsnplink.subsnp_id
		};
		if ($self->{'limit'}) {
		  $stmt .= qq{    
			      ORDER BY
				sorting_id ASC  
			     };
		}
		
		dumpSQL($self->{'dbSNP'}, $stmt,  $self->{source_engine} );
  
  
		$self->{'dbVar'}->do(qq{CREATE TABLE tmp_seq_$type\_$partition (
									  subsnp_id int NOT NULL,
									  seq_type int NOT NULL,
									  line_num int,
									  line varchar(255),
										revcom tinyint
									 )
										MAX_ROWS = 100000000 });
                print $logh Progress::location();
  
               load($self->{'dbVar'}, "tmp_seq_$type\_$partition", "subsnp_id", "seq_type", "line_num", "line");
               print $logh Progress::location();
               $self->{'dbVar'}->do("CREATE INDEX subsnp_id_idx on tmp_seq_$type\_$partition (subsnp_id)");  


		# merge the tables into a single tmp table
		$self->{'dbVar'}->do(qq{INSERT INTO tmp_seq (variation_id, subsnp_id, seq_type, line_num, type, line, revcom)
					  SELECT vs.variation_id, ts.subsnp_id, ts.seq_type, ts.line_num, '$type',
					  ts.line, ts.revcom
					  FROM   tmp_seq_$type\_$partition ts, variation_synonym vs
					  WHERE  vs.subsnp_id = ts.subsnp_id});
  print $logh Progress::location();
		#drop tmp table to free space
		$self->{'dbVar'}->do(qq{DROP TABLE tmp_seq_$type\_$partition});
  print $logh Progress::location();
	  }
	}
  }
  
  ## other species no partitions
  else {
  
    foreach my $type ('3','5') {
      debug(localtime() . "\tDumping $type' flanking sequence");
      
      $stmt = "SELECT ";
      if ($self->{'limit'}) {
	$stmt .= "TOP $self->{'limit'} ";
      }
      $stmt .= qq{
		    seq.subsnp_id AS sorting_id,
		    seq.type,
		    seq.line_num, 
		    seq.line,
        subsnplink.substrand_reversed_flag
		  FROM 
		    SubSNPSeq$type seq, 
		    SNP snp,
        SNPSubSNPLink subsnplink
		  WHERE 
		    snp.exemplar_subsnp_id = seq.subsnp_id AND
        seq.subsnp_id = subsnplink.subsnp_id
		 };
       if ($self->{'limit'}) {
	 $stmt .= qq{    
		     ORDER BY
		       sorting_id ASC  
		    };
       }
      dumpSQL($self->{'dbSNP'},$stmt,  $self->{source_engine} );
      
  
      $self->{'dbVar'}->do(qq{CREATE TABLE tmp_seq_$type (
								subsnp_id int NOT NULL,
								seq_type int NOT NULL,
								line_num int,
								line varchar(255),
								revcom tinyint,
								KEY subsnp_id_idx(subsnp_id)
							) 
							MAX_ROWS = 100000000 });
      print $logh Progress::location();
      
      load($self->{'dbVar'}, "tmp_seq_$type", "subsnp_id", "seq_type", "line_num", "line", "revcom");
      print $logh Progress::location();
  
      # merge the tables into a single tmp table
      $self->{'dbVar'}->do(qq{INSERT INTO tmp_seq (variation_id, subsnp_id, seq_type,
							 line_num, type, line, revcom)
				    SELECT vs.variation_id, ts.subsnp_id, ts.seq_type, ts.line_num, '$type',
				    ts.line, ts.revcom
				    FROM   tmp_seq_$type ts, variation_synonym vs
				    WHERE  vs.subsnp_id = ts.subsnp_id});
  print $logh Progress::location();
      #drop tmp table to free space
      $self->{'dbVar'}->do(qq{DROP TABLE tmp_seq_$type});
  print $logh Progress::location();
    }
  }
      
  $self->{'dbVar'}->do("ALTER TABLE tmp_seq ADD INDEX idx (subsnp_id, type, seq_type, line_num)");

  print $logh Progress::location();

  my $sth = $self->{'dbVar'}->prepare(qq{SELECT ts.variation_id, ts.subsnp_id, ts.type,
					       ts.line, ts.revcom
					       FROM   tmp_seq ts FORCE INDEX (idx)
					       ORDER BY ts.subsnp_id, ts.type, ts.seq_type, ts.line_num},{mysql_use_result => 1});
  
  $sth->execute();

  my ($vid, $ssid, $type, $line, $revcom);

  $sth->bind_columns(\$vid, \$ssid, \$type, \$line, \$revcom);

  open(FH, ">" . $self->{'tmpdir'} . "/" . $self->{'tmpfile'});
  my $upstream = '';
  my $dnstream = '';
  my $cur_vid;
  my $cur_revcom;

#=head Not flipping & merging flanks as part of import process

  debug(localtime() . "\tRearranging flanking sequence data");


  # dump sequences to file that can be imported all at once
  while($sth->fetch()) {
    if(defined($cur_vid) && $cur_vid != $vid) {
      # if subsnp in reverse orientation to refsnp,
      # reverse compliment flanking sequence
      if($cur_revcom) {
	($upstream, $dnstream) = ($dnstream, $upstream);
	reverse_comp(\$upstream);
	reverse_comp(\$dnstream);
      }

      print FH join("\t", $cur_vid, $upstream, $dnstream), "\n";
      
      $upstream = '';
      $dnstream = '';
    }

    $cur_vid   = $vid;
    $cur_revcom = $revcom;
    
    if($type == 5) {
      $upstream .= $line;
    } else {
      $dnstream .= $line;
    }
    
  }

  # do not forget last row...
  if($cur_revcom) {
    ($upstream, $dnstream) = ($dnstream, $upstream);
    reverse_comp(\$upstream);
    reverse_comp(\$dnstream);
  }
  print FH join("\t", $cur_vid, $upstream, $dnstream), "\n";
  
  $sth->finish();
  
  close FH;
  #$self->{'dbVar'}->do("DROP TABLE tmp_seq");
  print $logh Progress::location();

  debug(localtime() . "\tLoading flanking sequence data");

  # import the generated data
  load($self->{'dbVar'},"flanking_sequence","variation_id","up_seq","down_seq");
  print $logh Progress::location();

  unlink($self->{'tmpdir'} . "/" . $self->{'tmpfile'});
#=cut
  return;
}



sub variation_feature {
  my $self = shift;

  ## new mysql version errors with empty not null columns
  ## switch to null allowable here then back to non null in QC process 
  $self->{'dbVar'}->do("alter table variation_feature modify column map_weight int default null");

  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  debug(localtime() . "\tDumping seq_region data");

  dumpSQL($self->{'dbCore'}->db_handle, qq{SELECT sr.seq_region_id, sr.name, sr.coord_system_id
 				      FROM   seq_region sr, coord_system cs
                                      WHERE  sr.coord_system_id = cs.coord_system_id
                                      AND    cs.attrib like "%default_version%"});
    
  debug(localtime() . "\tLoading seq_region data");
  load($self->{'dbVar'}, "seq_region", "seq_region_id", "name");
  print $logh Progress::location();

  debug(localtime() . "\tDumping SNPLoc data");

  my ($tablename1,$tablename2,$row);

  my $version = substr($self->{'dbSNP_version'} ,1);
  if(  $version < 137){ ## table rename for dbSNP - keeping this temporarily for backwards comparibility

## 201206 - dbSNP no longer using assembly version in table names
## - leaving this temporarily for backwards compatibility

 # my ($assembly_version) =  $self->{'assembly_version'} =~ /^[a-zA-Z]+(\d+)\.*.*$/;
# override for platypus
 # $assembly_version = 1 if $self->{'dbm'}->dbCore()->species =~ /ornith/i;
  #$assembly_version = 3 if $self->{'dbm'}->dbCore()->species =~ /rerio/i;

        my $stmt = qq{
                SELECT 
                  name 
                FROM 
                  $self->{'snp_dbname'}..sysobjects 
                WHERE 
                  name LIKE '$self->{'dbSNP_version'}\_SNPContigLoc\__'
        };
        my $sth = $self->{'dbSNP'}->prepare($stmt);
        $sth->execute();
        
        while($row = $sth->fetchrow_arrayref()) {
            next if $row->[0] =~/Locus/;
            $tablename1 = $row->[0];
        }
        

        $stmt = qq{
             SELECT 
               name 
             FROM 
               $self->{'snp_dbname'}..sysobjects 
             WHERE 
               name LIKE '$self->{'dbSNP_version'}\_ContigInfo%'
          };
	my $sth1 = $self->{'dbSNP'}->prepare($stmt);
	$sth1->execute();
	
	while($row = $sth1->fetchrow_arrayref()) {
	    $tablename2 = $row->[0];
	}
    }
  else{
      $tablename1 = $self->{'dbSNP_version'} . "_SNPContigLoc" ;
      $tablename2 = $self->{'dbSNP_version'} . "_ContigInfo";      
  }
	##                     SNPContigLoc               ContigInfo
	debug(localtime() . "\ttable_name1 is $tablename1 table_name2 is $tablename2");
  #note the contig based cordinate is 0 based, ie. start at 0, lc_ngbr+2, t1.rc_ngbr
	
  my $stmt = "SELECT ";
	if ($self->{'limit'}) {
    $stmt .= "TOP $self->{'limit'} ";
  }
  $stmt .= qq{
                t1.snp_id AS sorting_id, 
                t2.contig_acc,
                t1.lc_ngbr+2,t1.rc_ngbr,
                CASE WHEN
                  t1.orientation = 1
                THEN
                  -1
                ELSE
                  1
                END,
                t1.aln_quality
              FROM   
                $tablename1 t1, 
                $tablename2 t2
              WHERE 
                t1.ctg_id=t2.ctg_id
	     };
     if ($self->{'limit'}) {
       $stmt .= qq{    
		   ORDER BY
		     sorting_id ASC  
	          };
     }
      dumpSQL($self->{'dbSNP'},$stmt, $self->{source_engine});
    
    
  debug(localtime() . "\tLoading SNPLoc data");

  create_and_load($self->{'dbVar'}, "tmp_contig_loc", "snp_id i* not_null", "contig * not_null", "start i", 
		  "end i", "strand i", "aln_quality d");
  print $logh Progress::location();
    
  debug(localtime() . "\tCreating genotyped variations");
    
  #creating the temporary table with the genotyped variations
  $self->{'dbVar'}->do(qq{CREATE TABLE tmp_genotyped_var SELECT DISTINCT variation_id FROM tmp_sample_genotype_single_bp});
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX variation_idx ON tmp_genotyped_var (variation_id)});
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{INSERT IGNORE INTO  tmp_genotyped_var SELECT DISTINCT variation_id FROM sample_genotype_multiple_bp});
  print $logh Progress::location();

    
  debug(localtime() . "\tCreating tmp_variation_feature data");
    
  dumpSQL($self->{'dbVar'},qq{SELECT v.variation_id, ts.seq_region_id, tcl.start, tcl.end, tcl.strand, v.name, v.source_id, tcl.aln_quality, v.somatic, v.class_attrib_id
					  FROM variation v, tmp_contig_loc tcl, seq_region ts
					  WHERE v.snp_id = tcl.snp_id
					  AND ts.name = tcl.contig});

  create_and_load($self->{'dbVar'},'tmp_variation_feature',"variation_id i* not_null","seq_region_id i", "seq_region_start i", "seq_region_end i", "seq_region_strand i", "variation_name", "source_id i not_null", "aln_quality d", "somatic i", "class_attrib_id i");
  print $logh Progress::location();

  debug(localtime() . "\tDumping data into variation_feature table");
  $self->{'dbVar'}->do(qq{INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,variation_name, flags, source_id, alignment_quality, somatic, class_attrib_id)
				SELECT tvf.variation_id, tvf.seq_region_id, tvf.seq_region_start, tvf.seq_region_end, tvf.seq_region_strand,tvf.variation_name,IF(tgv.variation_id,'genotyped',NULL), tvf.source_id, tvf.aln_quality, tvf.somatic, tvf.class_attrib_id
				FROM tmp_variation_feature tvf LEFT JOIN tmp_genotyped_var tgv ON tvf.variation_id = tgv.variation_id
				 });
  print $logh Progress::location();
    
    $self->{'dbVar'}->do("DROP TABLE tmp_contig_loc");
    $self->{'dbVar'}->do("DROP TABLE tmp_genotyped_var");
    $self->{'dbVar'}->do("DROP TABLE tmp_variation_feature");
}


#The task of getting the genotype will be chunked up and distributed to the farm in small pieces. Each result will be written to a loadfile that in the end
# will be used to populate the tmp_individual.. tables after everything has finished. We will chunk up the results by a) chromosome and b) individual. Because
#the number of genotypes per individual varies _drastically_ with 1000 genomes data having lots of genotypes and most other data only having a few, we need
#to determine the best way to chunk up the individuals in order to get an even distribution of the workload
sub parallelized_individual_genotypes {
  my $self = shift;
  my $load_only = shift;

  debug(localtime() . "\tStarting parallelized_individual_genotypes");
  my $genotype_table = 'tmp_sample_genotype_single_bp';
  my $multi_bp_gty_table = 'sample_genotype_multiple_bp';
  my $jobindex;


 #Get the create statement for tmp_sample_genotype_single_bp from master schema. We will need this to create the individual chromosome tables
  my $ind_gty_stmt = get_create_statement($genotype_table,$self->{'schema_file'});
  ## set default db engine
  $ind_gty_stmt =~ s/\;//;
  $ind_gty_stmt .= " ENGINE = MyISAM ;" ;

 
  my $failure_recovery ;   ## set to allow re-running of individual export jobs after memory failure

  #Put the log filehandle in a local variable
  my $logh = $self->{'log'}; 


  my $task_manager_file = 'individual_genotypes_task_management.txt';

 # Get the SubInd tables that may be split by chromosomes
  my $sql;
  if($self->{source_engine} =~/mssql|sqlserver/ ){
      $sql = qq[ SELECT name 
                 FROM $self->{'snp_dbname'}..sysobjects 
                 WHERE name LIKE 'SubInd%'SELECT OBJECT_ID('SubInd%')];
    }
  elsif($self->{source_engine} =~/postgreSQL/ ){
      $sql = qq[ SELECT tablename 
                 FROM pg_catalog.pg_tables
                 WHERE schemaname = '$self->{schema_name}' 
                 AND tablename LIKE 'subind%' ];
      warn "doing $sql\n";
    }
    else{
	$sql = qq[show tables like 'SubInd%'];
    }
  my @subind_tables = map {$_->[0]} @{$self->{'dbSNP'}->db_handle()->selectall_arrayref($sql)};
  print $logh "Found subind tables for genotypes: " . join ",", @subind_tables . "\n";

 # Prepared statement to find out if a table exists
 my $table2_ext_stmt = qq{  SHOW TABLES LIKE ?  };
 my $table_sth = $self->{'dbVar'}->prepare($table2_ext_stmt);
 
 # Use one common file for alleles and one for samples.
 my $file_prefix = $self->{'tmpdir'} . '/sample_genotypes';
 #Multi-bp genotypes will be written to a separate loadfile
 my $multi_bp_gty_file = $file_prefix . '_multi_bp_gty';
 # Allele file is used for caching the alleles
 my $allele_file = $file_prefix . '_alleles.txt';
 # Sample file is used for caching the samples
 my $sample_file = $file_prefix . '_individual_samples.txt';
  
    # For each SubInd_ch.. table, determine the best chunking of the individuals. We aim to get $target_rows number of rows to work on for each subtask but this will just be an approximate number.
  my $target_rows = 1e6; ## reduced from 10e6 

  my $genotype_counts;
  my %gty_tables;
  #my @skip_tables;
  my $sth;
  print $logh Progress::location() . "\tDividing the import task into suitable chunks\n";
  foreach my $subind_table (@subind_tables) {
    
    print $logh Progress::location() . "\t\tProcessing $subind_table\n";
    #The subtable to store the data for this subind table and the loadfile to use. The mapping file is used to temporarily store the subsnp_id -> variation_id mapping
    my $dst_table = "tmp_sample_genotype_single_bp\_$subind_table";
    my $loadfile = $file_prefix . '_' . $dst_table . '.txt';
    my $mapping_file = $file_prefix . '_subsnp_mapping_' . $subind_table . '.txt';
    #warn "setting up gty_tables: $subind_table=> $dst_table,$loadfile,$mapping_file\n";
    $gty_tables{$subind_table} = [$dst_table,$loadfile,$mapping_file];
    
    # If the subtable already exists, warn about this but skip the iteration (perhaps we are resuming after a crash)
    $table_sth->execute($dst_table);
    print $logh Progress::location();
    if (defined($table_sth->fetchrow_arrayref())) {
      warn("Table $dst_table already exists in destination database. Will skip importing from $subind_table");
     # push(@skip_tables,$subind_table);
      delete $gty_tables{$subind_table};
      next;
    }
  } 
   ### write task file
   unless($load_only){   

    unless(defined $failure_recovery  && $failure_recovery == 1){
     
     ## set up task file unless rerunning
     ($jobindex,$task_manager_file) = $self->create_parallelized_individual_genotypes_task_file(\%gty_tables, $target_rows);   

     return unless defined  $jobindex;  ## there may be nothing to do
    }

    $task_manager_file = "individual_genotypes_task_management.txt";
    # Run the job on the farm
    print $logh Progress::location() . "\tSubmitting the importing of the genotypes to the farm\n";
    my $jobname = 'individual_genotypes';

    my $result = $self->run_on_farm($jobname,$file_prefix,'calculate_gtype',$task_manager_file,1,$jobindex);

    my $jobid = $result->{'jobid'};
      
    # Check if any subtasks generated errors
    if ((my @error_subtasks = grep($result->{'subtask_details'}{$_}{'generated_error'},keys(%{$result->{'subtask_details'}})))) {
      warn($result->{'message'});
      print $logh Progress::location() . "\t " . $result->{'message'};
      print $logh Progress::location() . "\tThe following subtasks generated errors for job $jobid:\n";
      foreach my $index (@error_subtasks) {
	print $logh qq{\t\t$jobname\[$index\] ($jobid\[$index\])\n};
      }
    }

    # Check the result, if anything failed
    if (!$result->{'success'}) {
      warn($result->{'message'});
      print $logh Progress::location() . "\t " . $result->{'message'};
      print $logh Progress::location() . "\tThe following subtasks failed for unknown reasons for job $jobid:\n";
      foreach my $index (grep($result->{'subtask_details'}{$_}{'fail_reason'} =~ m/UNKNOWN/,keys(%{$result->{'subtask_details'}}))) {
	print $logh qq{\t\t$jobname\[$index\] ($jobid\[$index\])\n};
      }
    
      print $logh Progress::location() . "\tThe following subtasks failed because they ran out of resources for job $jobid:\n";
      foreach my $index (grep($result->{'subtask_details'}{$_}{'fail_reason'} =~ m/OUT_OF_[MEMORY|TIME]/,keys(%{$result->{'subtask_details'}}))) {
	print $logh qq{\t\t$jobname\[$index\] ($jobid\[$index\])\n};
      }
    }
    
    # If we still have subtasks that fail, this needs to be resolved before proceeding
    die("Some subtasks are failing (see log output). This needs to be resolved before proceeding with the loading of genotypes!") unless ($result->{'success'});
   }

  ## merge seperate export files into table specific load files
  ## db storage of dump file => load file link ( expect loading of larger files to be more efficient) 
  my $file_data_ext_sth = $self->{'dbVar'}->prepare(qq[SELECT sub_file, destination_file from tmp_indgeno_file]);
  $file_data_ext_sth->execute()||die "Problem extracting list of files to merge\n";
  my $file_list =  $file_data_ext_sth->fetchall_arrayref();
#  print "Cat'ing files\n"; 
  foreach my $pair (@{$file_list}){
      unless (-e $pair->[0] ){
        warn "No file of name $pair->[0] created\n" unless $pair->[0] =~/multi/;  ## many bins have multi files missing - report only missing singles 
        next;
      }
      open my $load_file, ">>", $pair->[1] || die "Failed to open $pair->[1] to write: $!\n";
#      print "Cat'ing $pair->[0] to $pair->[1]\n";
      open my $subfile, "<", $pair->[0] ||die  "Failed to open $pair->[0] to read: $!\n";
      while(<$subfile>){print $load_file $_ ;} 
      close $load_file;  
      close $subfile;
  }
  print "Load files written\n";

  # Loop over the subtables and load each of them
  my @subtables;
  my $jobindex = 0;  ## reset for loading jobs
  my $task_manager_file = $file_prefix . '_task_management_load.txt';
  open(MGMT,'>',$task_manager_file);
  foreach my $subind_table (keys(%gty_tables)) {
    
    my $dst_table = $gty_tables{$subind_table}->[0];
    my $loadfile = $gty_tables{$subind_table}->[1];
    
    # Skip if necessary
   # next if (grep(/^$subind_table$/,@skip_tables));
    
    # Create the sub table
    my $stmt = $ind_gty_stmt;
    $stmt =~ s/$genotype_table/$dst_table/;
    $self->{'dbVar'}->do($stmt);
    print $logh Progress::location();
    
    # If the loadfile doesn't exist, we can't load anything
    next unless (-e $loadfile);
    
    $jobindex++;
    print "Writing to load file $jobindex $loadfile $dst_table variation_id subsnp_id sample_id allele_1 allele_2\n";
    print MGMT qq{$jobindex $loadfile $dst_table variation_id subsnp_id sample_id allele_1 allele_2\n};
  }
  # Include the multiple bp genotype file here as well
  if (-e $multi_bp_gty_file) {
    $jobindex++;
    print MGMT qq{$jobindex $multi_bp_gty_file $multi_bp_gty_table variation_id subsnp_id sample_id allele_1 allele_2\n};
  }
  close(MGMT);
  
  # Run the job on the farm
  print $logh Progress::location() . "\tSubmitting loading of the genotypes to the farm\n";
  my $jobname = 'individual_genotypes_load';
  my $result = $self->run_on_farm($jobname,$file_prefix,'load_data_infile',$task_manager_file,1,$jobindex);
  my $jobid = $result->{'jobid'};
  
  # Check if any subtasks generated errors
  if ((my @error_subtasks = grep($result->{'subtask_details'}{$_}{'generated_error'},keys(%{$result->{'subtask_details'}})))) {
    warn($result->{'message'});
    print $logh Progress::location() . "\t " . $result->{'message'};
    print $logh Progress::location() . "\tThe following subtasks generated errors for job $jobid:\n";
    foreach my $index (@error_subtasks) {
      print $logh qq{\t\t$jobname\[$index\] ($jobid\[$index\])\n};
    }
  }
  # Check the result, if anything failed
  if (!$result->{'success'}) {
    warn($result->{'message'});
    print $logh Progress::location() . "\t " . $result->{'message'};
    print $logh Progress::location() . "\tThe following subtasks failed for unknown reasons for job $jobid:\n";
    foreach my $index (grep($result->{'subtask_details'}{$_}{'fail_reason'} =~ m/UNKNOWN/,keys(%{$result->{'subtask_details'}}))) {
      print $logh qq{\t\t$jobname\[$index\] ($jobid\[$index\])\n};
    }
  
    print $logh Progress::location() . "\tThe following subtasks failed because they ran out of resources for job $jobid:\n";
    foreach my $index (grep($result->{'subtask_details'}{$_}{'fail_reason'} =~ m/OUT_OF_[MEMORY|TIME]/,keys(%{$result->{'subtask_details'}}))) {
      print $logh qq{\t\t$jobname\[$index\] ($jobid\[$index\])\n};
    }
  }
  # If we still have subtasks that fail, this needs to be resolved before proceeding
  die("Some subtasks are failing (see log output). This needs to be resolved before proceeding with the loading of genotypes!") unless ($result->{'success'});
  
  #Drop the tmp_individual.. table if it exists
 my $drop_tmp_stmt = qq{
    DROP TABLE IF EXISTS
      $genotype_table
  };
  $self->{'dbVar'}->do($drop_tmp_stmt);
  print $logh Progress::location();
  
  # Merge all the subtables into the final table, or if there is just one subtable, rename it to the final table name
  print $logh Progress::location() . "\tMerging the genotype subtables into a big $genotype_table table\n";
  my $merge_subtables = join(",",map {$gty_tables{$_}->[0]} keys(%gty_tables));
  my  $stmt;
 if (scalar(keys(%gty_tables)) > 1) {
    
    #Add an empty table where any subsequent inserts will end up
    my $extra_table = $genotype_table . '_extra';
    $stmt = $ind_gty_stmt;
    $stmt =~ s/$genotype_table/$extra_table/;
    debug(localtime() . "\tDoing  $stmt in individual_genotypes");
    $self->{'dbVar'}->do($stmt);
    print $logh Progress::location();
    $merge_subtables .= ",$extra_table";
    
    $stmt = $ind_gty_stmt;
    $stmt =~ s/ENGINE = MyISAM ;//;  
    $stmt .= " ENGINE=MERGE INSERT_METHOD=LAST UNION=($merge_subtables)";
  }
  else {
    $stmt = qq{
      RENAME TABLE
	$merge_subtables
      TO
	$genotype_table
    };
  }

  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();
  
=head no longer needed as all go to the same compressed table
  #Move multiple bp genotypes into the multiple bp table
  $stmt = qq{
      SELECT DISTINCT
          variation_id
      FROM
          $multi_bp_gty_table
  };
  dumpSQL($self->{'dbVar'},$stmt);
  print $logh Progress::location();
  
  create_and_load($self->{'dbVar'},"tmp_multiple_bp_gty_variations","variation_id i* not_null"); 
  print $logh Progress::location();

  $stmt = qq{
      INSERT INTO
          $multi_bp_gty_table (
              variation_id,
              subsnp_id,
              individual_id,
              allele_1,
              allele_2
          )
      SELECT
          s.variation_id,
          s.subsnp_id,
          s.individual_id,
          s.allele_1,
          s.allele_2
      FROM
          tmp_multiple_bp_gty_variations t JOIN
          $genotype_table s ON (
              s.variation_id = t.variation_id
          )
  };
  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();

  $stmt = qq{
      DELETE FROM
          s
      USING
          tmp_multiple_bp_gty_variations t JOIN
          $genotype_table s ON (
              s.variation_id = t.variation_id
          )
  };
  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();
  
  $stmt = qq{
      DROP TABLE
          tmp_multiple_bp_gty_variations
  };
  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();
=cut  

}

sub create_parallelized_individual_genotypes_task_file{

  my $self            = shift;
  my $gty_tables      = shift;
  my $target_rows     = shift;

  my $logh = $self->{'log'};

  print $logh Progress::location();
  my $task_manager_file  = 'individual_genotypes_task_management.txt';

  #Multi-bp genotypes will be written to a separate loadfile
  my $multi_bp_gty_file =  $self->{'tmpdir'} . '/sample_genotypes_multi_bp_gty';

  # Store the data for the farm submission in a hash with the data as key. This way, the jobs will be "randomized" w.r.t. chromosomes and chunks when submitted so all processes shouldn't work on the same tables/files at the same time.
  my %job_data;
  debug(localtime() . "\tCreating parallelized_individual_genotypes task file");
  ## temp table to hold files to create & load
  $self->{'dbVar'}->do(qq[ DROP TABLE IF EXISTS tmp_indgeno_file ]);
  $self->{'dbVar'}->do(qq[ CREATE TABLE tmp_indgeno_file (sub_file varchar(255), destination_file varchar(255) ) ] );
  my $file_ins_sth = $self->{'dbVar'}->prepare(qq[insert into tmp_indgeno_file (sub_file ,destination_file) values (?,?)]);

  foreach my $subind_table (keys(%$gty_tables)) {
    
  
  my $stmt = qq{
    SELECT
 	submitted_ind_id,
 	COUNT(*)
    FROM
 	$subind_table
    GROUP BY
 	submitted_ind_id
    ORDER BY
 	submitted_ind_id ASC
  };
  my $count_lines_sth = $self->{'dbSNP'}->prepare($stmt);
  $count_lines_sth->execute();
  print $logh Progress::location();
  my ($submitted_ind_id,$genotype_count);
  $count_lines_sth->bind_columns(\$submitted_ind_id,\$genotype_count);
  
  #Loop over the counts and when a large enough number of genotypes have been counted, split that off as a chunk to be submitted to the farm
  my $total_count = 0;
  my $start_id = -1;
  while ($count_lines_sth->fetch()) {
    $total_count += $genotype_count;
    # Set the start id if it's not specified
    $start_id = $submitted_ind_id if ($start_id < 0);
    #Break off the chunk if it's large enough
    if ($total_count >= $target_rows) {
 	$job_data{qq{$subind_table $start_id $submitted_ind_id}}++;
 	$total_count = 0;
 	$start_id = -1;
    }
  }
  #Add the rest of the individual_ids if a new chunk has been started
  $job_data{qq{$subind_table $start_id $submitted_ind_id}}++ if ($start_id >= 0);
  print $logh Progress::location();
  }
  
  # If there are no genotypes at all to be imported, just return here
  return unless (scalar(keys(%job_data)));
  
  # Print the job parameters to a file

  my $jobindex = 0;
  open(MGMT,'>',$task_manager_file);
  foreach my $params (keys(%job_data)) {
    $jobindex++;
    my @arr = split(/\s+/,$params);

    print $logh "writing genotype task management file : $params\n";  
    $file_ins_sth->execute("$$gty_tables{$arr[0]}->[1]\_$jobindex", $$gty_tables{$arr[0]}->[1])|| die "Problem entering geno file info\n";
    $file_ins_sth->execute("$multi_bp_gty_file\_$jobindex", $multi_bp_gty_file)|| die "Problem entering geno file info\n";
    print MGMT qq{$jobindex $arr[0] $$gty_tables{$arr[0]}->[1]\_$jobindex $multi_bp_gty_file\_$jobindex $arr[1] $arr[2] $$gty_tables{$arr[0]}->[2]\n};
  }
  close(MGMT);

  return ($jobindex,$task_manager_file);
    
}

#
# loads population genotypes into the population_genotype table
#
sub population_genotypes {
    my $self = shift;

  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};

  ## revert schema to letter rather than number coded genotype
  $self->{'dbVar'}->do( qq[drop table population_genotype]);

  $self->{'dbVar'}->do( qq[CREATE TABLE  population_genotype (
                          population_genotype_id int(10) unsigned NOT NULL AUTO_INCREMENT,
                          variation_id  int(10) unsigned NOT NULL,
                          subsnp_id int(15) unsigned DEFAULT NULL,
                          allele_1 varchar(25000) DEFAULT NULL,
                          allele_2 varchar(25000) DEFAULT NULL,
                          frequency float DEFAULT NULL,
                          population_id int(10) unsigned DEFAULT NULL,
                          count int(10) unsigned DEFAULT NULL,
                          PRIMARY KEY (population_genotype_id),
                          KEY variation_idx (variation_id),
                          KEY subsnp_idx (subsnp_id),
                          KEY population_idx (population_id))
                          engine=MyISAM
                         ]);


    my $stmt;
    my $allele_table_ref = $self->{'dbVar'}->db_handle->selectall_arrayref(qq{show tables like "tmp_rev_allele"});
    my $allele_table = $allele_table_ref->[0][0];
    if (! $allele_table) {
      debug(localtime() . "\tDumping allele data");
      $stmt = qq{
                 SELECT 
                   a1.allele_id, 
                   a1.allele, 
                   a2.allele
		 FROM 
		   $self->{'dbSNP_share_db'}.Allele a1, 
		   $self->{'dbSNP_share_db'}.Allele a2
                 WHERE 
                   a1.rev_allele_id = a2.allele_id
                };
      dumpSQL($self->{'dbSNP'},$stmt, $self->{source_engine});
      create_and_load($self->{'dbVar'}, "tmp_rev_allele", "allele_id i* not_null","allele * not_null", "rev_allele not_null");
      print $logh Progress::location();
      
    }
    debug(localtime() . "\tDumping GtyFreqBySsPop and UniGty data");
 
    $stmt = "SELECT ";
    if ($self->{'limit'}) {
      $stmt .= "TOP $self->{'limit'} ";
    }
#ROUND(gtfsp.cnt,0),
    $stmt .= qq{
                  gtfsp.subsnp_id AS sorting_id, 
                  gtfsp.pop_id, 
                  gtfsp.freq,
                  gtfsp.cnt,
                  a1.allele, 
                  a2.allele
	       	FROM   
	       	  GtyFreqBySsPop gtfsp, 
	       	  $self->{'dbSNP_share_db'}.UniGty ug, 
	       	  $self->{'dbSNP_share_db'}.Allele a1, 
	       	  $self->{'dbSNP_share_db'}.Allele a2
	       	WHERE  
	       	  gtfsp.unigty_id = ug.unigty_id AND    
	       	  ug.allele_id_1 = a1.allele_id AND    
	       	  ug.allele_id_2 = a2.allele_id
	       };
     if ($self->{'limit'}) {
       $stmt .= qq{    
		   ORDER BY
		     sorting_id ASC  
	          };
     }
    dumpSQL($self->{'dbSNP'},$stmt, $self->{source_engine});

     debug(localtime() . "\tloading population_genotype data");

     create_and_load($self->{'dbVar'}, "tmp_pop_gty", 'subsnp_id i* not_null unsigned', 'pop_id i* not_null', 'freq','count','allele_1', 'allele_2');
  print $logh Progress::location();

   $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX pop_genotype_idx ON population_genotype(variation_id,subsnp_id,frequency,population_id,allele_1(5),allele_2(5))});

  ## human is too big to do in one go - breaking up crudely 

  my $get_max_sth = $self->{'dbVar'}->prepare(qq[select max(variation_synonym_id) from variation_synonym]);
  $get_max_sth->execute()||die;
  my $max = $get_max_sth->fetchall_arrayref();

    if($max > 2000000){
	pop_geno_batched($self, $max);
    }
    else{

    ## This is done to remove duplicates
    $self->{'dbVar'}->do(qq{INSERT IGNORE INTO population_genotype (variation_id,subsnp_id,allele_1, allele_2, frequency, population_id, count)
			                  SELECT vs.variation_id,vs.subsnp_id,tra1.rev_allele as allele_1,tra2.rev_allele as allele_2,tg.freq,p.population_id, tg.count
					  FROM   variation_synonym vs, tmp_pop_gty tg,tmp_rev_allele tra1,tmp_rev_allele tra2, population p
					  WHERE  vs.subsnp_id = tg.subsnp_id
                                          AND    tg.allele_1 = tra1.allele
                                          AND    tg.allele_2 = tra2.allele
                                          AND    vs.substrand_reversed_flag = 1
					  AND    p.pop_id = tg.pop_id});



     print $logh Progress::location(); 
     $self->{'dbVar'}->do(qq{INSERT IGNORE INTO population_genotype (variation_id,subsnp_id,allele_1, allele_2, frequency, population_id, count)
					  SELECT vs.variation_id,vs.subsnp_id,tg.allele_1,tg.allele_2,tg.freq,p.population_id, tg.count
					  FROM   variation_synonym vs, tmp_pop_gty tg,population p
					  WHERE  vs.subsnp_id = tg.subsnp_id
                                          AND    vs.substrand_reversed_flag = 0
					  AND    p.pop_id = tg.pop_id
                                          });
    }
   print $logh Progress::location(); 

    $self->{'dbVar'}->do(qq{DROP INDEX pop_genotype_idx ON population_genotype});
    $self->{'dbVar'}->do("DROP TABLE tmp_pop_gty");
  debug(localtime() . "\tFinished population_genotype");
}

## crude break up of population_genotype import for larger data sets
sub pop_geno_batched{

  my $self = shift;
  my $max = shift;

  my $logh = $self->{'log'};
  print $logh Progress::location(); 

  my $batch = 1000000;
  my $start = 1;
  debug(localtime() . "\tStarting binned pop geno");
  while( $start  < $max->[0]->[0] ){	
     
     my $end = $start + $batch ;   
     print $logh Progress::location(); 

    ## This is done to remove duplicates
    $self->{'dbVar'}->do(qq{INSERT IGNORE INTO population_genotype (variation_id,subsnp_id,allele_1, allele_2, frequency, population_id, count)
			                  SELECT vs.variation_id,vs.subsnp_id,tra1.rev_allele as allele_1,tra2.rev_allele as allele_2,tg.freq,p.population_id, tg.count
					  FROM   variation_synonym vs, tmp_pop_gty tg,tmp_rev_allele tra1,tmp_rev_allele tra2, population p
					  WHERE  vs.subsnp_id = tg.subsnp_id
                                          AND    tg.allele_1 = tra1.allele
                                          AND    tg.allele_2 = tra2.allele
                                          AND    vs.substrand_reversed_flag = 1
					  AND    p.pop_id = tg.pop_id
                                          AND variation_synonym_id between $start and $end
});

    
     $self->{'dbVar'}->do(qq{INSERT IGNORE INTO population_genotype (variation_id,subsnp_id,allele_1, allele_2, frequency, population_id, count)
					  SELECT vs.variation_id,vs.subsnp_id,tg.allele_1,tg.allele_2,tg.freq,p.population_id, tg.count
					  FROM   variation_synonym vs, tmp_pop_gty tg,population p
					  WHERE  vs.subsnp_id = tg.subsnp_id
                                          AND    vs.substrand_reversed_flag = 0
					  AND    p.pop_id = tg.pop_id
                                          AND variation_synonym_id between $start and $end
                                          });

     $start = $end + 1;
    }
  print $logh "Completed batched population_genotype import\n";
  print $logh Progress::location(); 
  debug(localtime() . "\tFinished binned pop geno");
}

## Extract old rs ids to add to synonym table

sub archive_rs_synonyms {

    my $self = shift;

    return unless $self->table_exists_and_populated('RsMergeArch');

    debug(localtime() . "\tLooking for old rs");  

    ## Add source description only if old refSNP names found
    my $source_id = $self->source_table("Archive dbSNP");

    my $logh = $self->{'log'};
    
    print $logh Progress::location();

    my $concat_sql;
    if($self->{source_engine} =~/mssql/ ){
	$concat_sql = qq[ SELECT 'rs' + CAST(rsHigh as VARCHAR(20)),
                          'rs' + CAST(rsCurrent as VARCHAR(20)) , 
                          orien2Current,
                          'rs' + CAST(rsLow as VARCHAR(20)) 
                          FROM RsMergeArch ] ;
    }
    elsif($self->{source_engine} =~/postgreSQL/ ){
	$concat_sql = qq[ SELECT 'rs' || rsHigh,
                          'rs' || rsCurrent, 
                          orien2Current,
                          'rs' || rsLow 
                          FROM RsMergeArch ] ;
    }
    else{
	$concat_sql = qq[ SELECT CONCAT('rs', CAST(rsHigh as CHAR) ),
                          CONCAT('rs', CAST(rsCurrent as CHAR)) , 
                          orien2Current,
                          CONCAT('rs', CAST(rsLow as CHAR) )
                          FROM RsMergeArch ];
    }
    # export old rs id from dbSNP
    debug(localtime() . "\tExporting old rs with $concat_sql");  
    dumpSQL($self->{'dbSNP'}, $concat_sql , $self->{source_engine} ) ;

   #loading it to variation database in temp rshist table
   create_and_load( $self->{'dbVar'}, "rsHist", "rsHigh * not_null", "rsCurrent * not_null","orien2Current not_null", "rsLow") ;
   print $logh Progress::location();

   debug(localtime() . "\tAdding old rs as synonyms"); 

   my $get_max_sth = $self->{'dbVar'}->prepare(qq[select min(snp_id), max(snp_id) from variation ]);
   $get_max_sth->execute()||die;
   my $range = $get_max_sth->fetchall_arrayref();

   my $batch_size = 100000;
   my $start      = $range->[0]->[0];
   my $max        = $range->[0]->[1];
 

   # append to synonym table with database 'Archive dbSNP' - slow query binned

   while ( $start < $max ){

      my $end =  $start + $batch_size;

      $self->{'dbVar'}->do(qq{INSERT INTO variation_synonym (variation_id, source_id, name)
                             (SELECT v.variation_id, $source_id, r.rsHigh
                              FROM variation v, rsHist r
                              WHERE v.name = r.rsCurrent
                              AND v.snp_id between $start and $end)
                             }); 

      ## take from rsLow if rsCurrent has no id
      if($self->{source_engine} =~/mssql/ ){
	  $self->{'dbVar'}->do(qq{INSERT INTO variation_synonym (variation_id, source_id, name)
                                 (SELECT v.variation_id, $source_id, r.rsHigh
                                 FROM variation v, rsHist r
                                 WHERE v.name = r.rsLow
                                 AND  r.rsCurrent = 'rs'
                                 AND v.snp_id between $start and $end)
                             }); 
      }
      else{
	  ## for missing data in mysql mirror
	  $self->{'dbVar'}->do(qq{INSERT INTO variation_synonym (variation_id, source_id, name)
                                 (SELECT v.variation_id, $source_id, r.rsHigh
                                  FROM variation v, rsHist r
                                  WHERE v.name = r.rsLow                              
                                  AND  r.rsCurrent not like 'rs%'
                                  AND v.snp_id between $start and $end)
                             }); 
      }

      $start =  $end +1;
   }
    ## clean up temp table
    $self->{'dbVar'}->do(qq[drop table rsHist]);

    debug(localtime() . "\tArchive rs synonyms done");

}



# cleans up some of the necessary temporary data structures after the
# import is complete
sub cleanup {
    my $self = shift;

  #Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  
    debug(localtime() . "\tIn cleanup...");
    #remove populations that are not present in the Individual or Allele table for the specie
    $self->{'dbVar'}->do('CREATE TABLE tmp_pop (population_id int PRIMARY KEY)'); #create a temporary table with unique populations
  print $logh Progress::location();
    $self->{'dbVar'}->do('INSERT IGNORE INTO tmp_pop SELECT distinct(population_id) FROM allele'); #add the populations from the alleles
  print $logh Progress::location();
    $self->{'dbVar'}->do('INSERT IGNORE INTO tmp_pop SELECT distinct(population_id) FROM population_genotype'); #add the populations from the population_genotype
  print $logh Progress::location();
    $self->{'dbVar'}->do('INSERT IGNORE INTO tmp_pop SELECT population_id FROM individual_population'); #add the populations from the individuals
  print $logh Progress::location();
    $self->{'dbVar'}->do(qq{INSERT IGNORE INTO tmp_pop SELECT super_population_id 
 				      FROM population_structure ps, tmp_pop tp 
 				      WHERE tp.population_id = ps.sub_population_id}); #add the populations from the super-populations
  print $logh Progress::location();
    

    ### delete population entries without alleles, population_genotypes or individuals
    my $sql = qq{DELETE FROM p USING population p
		  LEFT JOIN tmp_pop tp ON p.population_id = tp.population_id
		  WHERE tp.population_id is null };

    ### delete sample_synonym entries without alleles, population_genotypes or individuals
    my $sql_2 = qq{DELETE FROM ps USING population_synonym ps
		    LEFT JOIN tmp_pop tp ON ps.population_id = tp.population_id
		    WHERE tp.population_id is null };


    $self->{'dbVar'}->do($sql); #delete from population
  print $logh Progress::location();
    # populations not present
    $self->{'dbVar'}->do($sql_2); #delete from population_synonym

  print $logh Progress::location();
  


    $self->{'dbVar'}->do('DROP TABLE tmp_pop'); #and finally remove the temporary table
  print $logh Progress::location();

    $self->{'dbVar'}->do('ALTER TABLE variation_synonym DROP COLUMN substrand_reversed_flag');
  print $logh Progress::location();
    $self->{'dbVar'}->do('ALTER TABLE population DROP COLUMN pop_class_id, DROP COLUMN pop_id');
    print $logh Progress::location();

    ## link between ensembl sample_id and dbSNP submitted_ind_id
    $self->{'dbVar'}->do('DROP TABLE tmp_ind');
    
    ## list of genotype files to create, check and load
    $self->{'dbVar'}->do('DROP TABLE tmp_indgeno_file');

    ## subsnp strand rs info from synonym creation
    $self->{'dbVar'}->do('DROP TABLE tmp_var_allele');
}

## add dbSNP HGVS synonyms
sub import_hgvs{

  my $self = shift;
  my $hgvs_extr_stmt = "select snp_id, hgvs_name from snp_hgvs where hgvs_name like 'NP%' or hgvs_name like 'NM%' ";

  dumpSQL($self->{'dbSNP'}, $hgvs_extr_stmt , $self->{source_engine} ) ;

  #load it to variation database in temp dbsnp_hgvs table
  create_and_load( $self->{'dbVar'}, "dbsnp_hgvs", "snp_id * not_null", "hgvs_name * not_null") ;

  my $get_source_sth = $self->{'dbVar'}->prepare(qq[select source_id from source where name ='dbSNP HGVS' ]);
  $get_source_sth->execute()||die;
  my $synon_source = $get_source_sth->fetchall_arrayref();

  $self->{'dbVar'}->do(qq[INSERT IGNORE INTO variation_synonym (variation_id,name,source_id)  (select variation.variation_id, dbsnp_hgvs.hgvs_name, $synon_source->[0]->[0]
                          from variation, dbsnp_hgvs
                          where variation.snp_id = dbsnp_hgvs.snp_id)
                         ]);

}

sub sort_num{

return $a<=>$b;
}


1;

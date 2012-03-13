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

use Progress;
use DBI qw(:sql_types);
use Fcntl qw( LOCK_SH LOCK_EX );
use List::Util qw ( min max );

our $FARM_BINARY = 'bsub';

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
    'max_concurrent_jobs' => 10,
    'memory' => 8000,
    'queue' => 1,
    'wait_queue' => 0
  },
  
  'allele_table_load' => {
    'script' => 'run_task.pl',
    'max_concurrent_jobs' => 5,
    'memory' => 8000,
    'queue' => 2,
    'wait_queue' => 0
  },
  
  'individual_genotypes' => {
    'script' => 'run_task.pl',
    'max_concurrent_jobs' => 5,
    'memory' => 8000,
    'queue' => 1,
    'wait_queue' => 0
  },
  
  'individual_genotypes_load' => {
    'script' => 'run_task.pl',
    'max_concurrent_jobs' => 5,
    'memory' => 8000,
    'queue' => 2,
    'wait_queue' => 0
  },
  
);

#ÊThe queues in ascending run time limit order
our @FARM_QUEUES = (
  'small',
  'normal',
  'long',
  'basement'
);

#ÊThe maximum amount of farm memory that we will request
our $MAX_FARM_MEMORY = 15900;
our $FARM_MEMORY_INCREMENT = 4000;

#creates the object and assign the attributes to it (connections, basically)
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbm, $tmp_dir, $tmp_file, $limit, $mapping_file_dir, $dbSNP_BUILD_VERSION, $ASSEMBLY_VERSION, $GROUP_TERM,$GROUP_LABEL, $skip_routines, $scriptdir, $logh) =
        rearrange([qw(DBMANAGER TMPDIR TMPFILE LIMIT MAPPING_FILE_DIR DBSNP_VERSION ASSEMBLY_VERSION GROUP_TERM GROUP_LABEL SKIP_ROUTINES SCRIPTDIR LOG)],@_);

#  my $dbSNP_share_db = "dbSNP_$dbSNP_BUILD_VERSION\_shared";
#  $dbSNP_share_db =~ s/dbSNP_b/dbSNP_/;
#  my $dbSNP_share_db = 'new_shared';
  my $dbSNP = $dbm->dbSNP()->dbc();
  my $dbCore = $dbm->dbCore()->dbc();
  my $dbVar = $dbm->dbVar()->dbc();
  my $snp_dbname = $dbSNP->dbname();
  my $species = $dbm->species();
  my $shared_db = $dbm->dbSNP_shared();
  my $registry_file = $dbm->registryfile();
  
  debug(localtime() . "\tThe shared database is $shared_db");

  return bless {'dbSNP' => $dbSNP,
		'dbCore' => $dbCore,
		'dbVar' => $dbVar, ##this is a dbconnection
		'snp_dbname' => $snp_dbname,
		'species' => $species,
		'tmpdir' => $tmp_dir,
		'tmpfile' => $tmp_file,
		'limit' => $limit,
    'mapping_file_dir' => $mapping_file_dir,
		'dbSNP_version' => $dbSNP_BUILD_VERSION,
		'dbSNP_share_db' => $shared_db,
		'assembly_version' => $ASSEMBLY_VERSION,
		'skip_routines' => $skip_routines,
		'log' => $logh,
		'registry_file' => $registry_file,
		'scriptdir' => $scriptdir,
		'dbm' => $dbm,
		'group_term'  => $GROUP_TERM,
		'group_label' => $GROUP_LABEL
		}, $class;
}

#main and only function in the object that dumps all dbSNP data 
sub dump_dbSNP{ 

  my $self = shift;

  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  #the following steps need to be run when initial starting the job. If job failed for some reason and some steps below are already finished, then can comment them out
  my @subroutines = (
    'create_coredb',
    'source_table',
    'population_table',
    'individual_table',
    'variation_table',
    'parallelized_individual_genotypes',
    'population_genotypes',
    'parallelized_allele_table',
    'flanking_sequence_table',
    'variation_feature',
    'cleanup'
  );
  
  #ÊThe GenericContig object has an array where routines that should be skipped can be specified. For now, add create_coredb and cleanup by default
  push(@{$self->{'skip_routines'}},('create_coredb','cleanup'));
  
  # When resuming after a crash, put already finished modules into this array 
  push(@{$self->{'skip_routines'}},());
  
  my $resume;
  
  my $clock = Progress->new();
  
  # Loop over the subroutines and run each one
  foreach my $subroutine (@subroutines) {
    
    #ÊCheck if this subroutine should be skipped
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
  
  #ÊPut the log filehandle in a local variable
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
  #ÊIf just the start was defined, it should be a list of subtasks that should be re-run
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
  
  #ÊWrap the command to be executed into a shell script
  my $task_command = qq{perl $self->{'scriptdir'}/$script -species $self->{'species'} -dbSNP_shared $self->{'dbSNP_share_db'} -registry_file $self->{'registry_file'} -task $task -file_prefix $file_prefix -task_management_file $task_manager_file -tempdir $self->{'tmpdir'} -tempfile $self->{'tmpfile'} } . join(" ",@args);
  my $script_wrapper = $file_prefix . '_command.sh';
  open(CMD,'>',$script_wrapper);
  flock(CMD,LOCK_EX);
  print CMD qq{#!/usr/local/bin/bash\n};
  print CMD qq{$task_command\n};
  close(CMD);
  print $logh Progress::location();
  
  my $bsub_cmd = qq{$FARM_BINARY -R'select[mem>$memory\] rusage[mem=$memory\]' -M$memory_long -q $queue -J'$jobname$array_options' -o $logfile_prefix\%J.%I.out -e $logfile_prefix\%J.%I.err bash $script_wrapper};
  
  #ÊSubmit the job array to the farm
  my $submission = `$bsub_cmd`;
  my ($jobid) = $submission =~ m/^Job \<([0-9]+)\>/i;
  print $logh Progress::location();
  
  #ÊSubmit a job that depends on the job array so that the script will halt
  system(qq{$FARM_BINARY -J $jobid\_waiting -q $wait_queue -w'ended($jobid)' -K -o $file_prefix\_waiting.out sleep 1});
  print $logh Progress::location();
  
  #ÊCheck the error and output logs for each subtask. If the error file is empty, delete it. If not, warn that the task generated errors. If the output file doesn't say that it completed successfully, report the job as unseccessful and report which tasks that failed
  my $all_successful = 1;
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
    
    #ÊDoes the outfile say that the process exited successfully? Faster than querying bhist... (?)
    $job_details{$index}->{'success'} = 0;
    open(FH,'<',$outfile);
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
  
  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  my $param_key = $jobname;
  if (!exists($FARM_PARAMS{$jobname})) {
    warn("No farm resource parameters defined for job $jobname, will use default parameters");
    $param_key = 'default';
  }
  
  my $max_time = 0;
  my $max_memory = 0;
  
  #ÊThis is a preparation step for re-running a farm job. What we do is to re-submit to a longer queue (if one is available) and request more memory
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
  
  #ÊIf we have already maxed out the resources, it won't help running the job again
  return undef if ($max_time && $max_memory);
  
  return $self->run_on_farm($jobname,$file_prefix . "_submission",$iteration,@args);
}

sub create_coredb {

  my $self = shift;
  
  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  my $coredb_name = $self->{'dbCore'}->dbc->dbname();
  $self->{'dbVar'}->do(qq{CREATE DATABASE $coredb_name});
  print $logh Progress::location();
  debug(localtime() . "\tmake sure create $coredb_name.coord_system");
  my $csid_ref = $self->{'dbCore'}->dbc->selectall_arrayref(qq{SELECT coord_system_id from coord_system WHERE name = 'chromosome' and attrib = 'default_version'});
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
    
  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
    #my ($dbname,$version) = split /\_/,$self->{'dbSNP'}->dbname(); #get the version of the dbSNP release
#    my ($dbname,$version) = split /\_/,$self->{'snp_dbname'};
    my ($species,$tax_id,$version) = $self->{'snp_dbname'} =~ m/^(.+)?\_([0-9]+)\_([0-9]+)$/;
    my $dbname = 'dbSNP';
    my $url = 'http://www.ncbi.nlm.nih.gov/projects/SNP/';
    $self->{'dbVar'}->do(qq{INSERT INTO source (source_id,name,version,description,url,somatic_status) VALUES (1,"$dbname",$version,"Variants (including SNPs and indels) imported from dbSNP", "$url", "mixed")});
    print $logh Progress::location();

}

sub table_exists_and_populated {

    # check if a table is present in dbSNP, and that it has some data in it
    
    my $self = shift;
    my $table = shift;

    my ($obj_id) = $self->{'dbSNP'}->db_handle->selectrow_array(qq{SELECT OBJECT_ID('$table')});
 
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

sub clin_sig {
    my $self = shift;

    return unless $self->table_exists_and_populated('SNPClinSig');

    my $logh = $self->{'log'};
    
    # dump the data from dbSNP and store it in a temporary clin_sig table
    # in the variation database

    print $logh Progress::location();
    debug(localtime() . "\tDumping clinical significance");

    my $stmt = qq{
        SELECT  cs.snp_id, csc.descrip 
        FROM    SNPClinSig cs, ClinSigCode csc
        WHERE   cs.clin_sig_id = csc.code
    };
    
    dumpSQL($self->{'dbSNP'}, $stmt);
    
    print $logh Progress::location();
    debug(localtime() . "\tLoading clinical significance");
    
    create_and_load( $self->{'dbVar'}, "clin_sig", "snp_id i*", "descrip");

    # check that we know about all the possible values
    
    $stmt = qq{
        SELECT  count(descrip) 
        FROM    clin_sig 
        WHERE   descrip NOT IN (
            SELECT  a.value 
            FROM    attrib a, attrib_type att 
            WHERE   a.attrib_type_id = att.attrib_type_id 
            AND     att.code = 'dbsnp_clin_sig'
        )
    };

    my ($count) = $self->{'dbVar'}->db_handle->selectrow_array($stmt);

    if ($count > 0) {   
        die "There are unexpected (probably new) clinical significance types, add them to the attrib table first";
    }

    # update the variation table with the correct attrib_ids

    print $logh Progress::location();
    debug(localtime() . "\tUpdating variation table with clinical significance");

    $stmt = qq{
        UPDATE  variation v, clin_sig c, attrib a, attrib_type att
        SET     v.clinical_significance_attrib_id = a.attrib_id
        WHERE   att.code = 'dbsnp_clin_sig'
        AND     a.attrib_type_id = att.attrib_type_id
        AND     c.descrip = a.value
        AND     v.name = CONCAT('rs', c.snp_id)
    };

    $self->{'dbVar'}->do($stmt);

    # in case dbsnp also assign clin sigs to synonyms we also join to the variation_synonym table
    # (for build 135 they (redundantly) had both the original and the synonym rsID associated with 
    # a clinical significance, but this may  change in the future as it has for MAF and suspect 
    # SNPs and it doesn't do any harm to run this statement)

    $stmt = qq{ 
        UPDATE  variation v, variation_synonym vs, clin_sig c, attrib a, attrib_type att
        SET     v.clinical_significance_attrib_id = a.attrib_id
        WHERE   att.code = 'dbsnp_clin_sig'
        AND     a.attrib_type_id = att.attrib_type_id
        AND     c.descrip = a.value
        AND     vs.name = CONCAT('rs', c.snp_id)
        AND     vs.variation_id = v.variation_id
    };

    $self->{'dbVar'}->do($stmt);

    #$self->{'dbVar'}->do(qq{DROP TABLE clin_sig});
}

sub minor_allele_freq {
    my $self = shift;

    return unless $self->table_exists_and_populated('SNPAlleleFreq_TGP');

    my $logh = $self->{'log'};
    
    # dump the data from dbSNP and store in the temporary maf table in
    # the variation database

    print $logh Progress::location();
    debug(localtime() . "\tDumping global minor allele freqs");

    my $shared = $self->{'dbSNP_share_db'};

    my $stmt = qq{
        SELECT  af.snp_id, a.allele, af.freq, af.count, af.is_minor_allele
        FROM    SNPAlleleFreq_TGP af, ${shared}..Allele a
        WHERE   af.allele_id = a.allele_id
    };

    dumpSQL($self->{'dbSNP'}, $stmt);
    
    print $logh Progress::location();
    debug(localtime() . "\tLoading global minor allele freqs");
    
    create_and_load($self->{'dbVar'}, "maf", "snp_id i*", "allele l", "freq f", "count i", "is_minor_allele i");
    
    print $logh Progress::location();
    debug(localtime() . "\tUpdating variations with global minor allele frequencies");
   
    my $variation_sql = qq{
        UPDATE  variation v, maf m
        SET     v.minor_allele_freq = m.freq, v.minor_allele_count = m.count, v.minor_allele = m.allele
        WHERE   v.name = CONCAT('rs', m.snp_id)
    };

    my $synonym_sql = qq{
        UPDATE  variation v, maf m, variation_synonym vs
        SET     v.minor_allele_freq = m.freq, v.minor_allele_count = m.count, v.minor_allele = m.allele
        WHERE   vs.name = CONCAT('rs', m.snp_id)
        AND     v.variation_id = vs.variation_id
    };

    # update the variation table with the frequencies for only the minor alleles

    $self->{'dbVar'}->do($variation_sql . " AND m.is_minor_allele");

    # it seems dbSNP also store frequencies on synonyms but don't copy this across to the
    # merged refsnp for some reason, we want to do this though so we also join to the 
    # variation_synonym table and copy across the data
    
    $self->{'dbVar'}->do($synonym_sql . " AND m.is_minor_allele");

    # we also copy across data where the MAF = 0.5 which do not have the is_minor_allele 
    # flag set, we will ensure this is set on the non-reference allele later in the post-processing
    # when we have the allele string from the variation_feature table to check which is the 
    # reference allele

    $self->{'dbVar'}->do($variation_sql . " AND m.freq = 0.5");

    $self->{'dbVar'}->do($synonym_sql . " AND m.freq = 0.5");

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

    dumpSQL($self->{'dbSNP'}, $stmt);
    
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

    $stmt = qq{
        INSERT INTO failed_variation (variation_id, failed_description_id)
        SELECT  v.variation_id, fd.failed_description_id
        FROM    suspect s, variation v, failed_description fd
        WHERE   fd.description = 'Flagged as suspect by dbSNP'
        AND     v.name = CONCAT('rs', s.snp_id)
    };

    $self->{'dbVar'}->do($stmt);

    # also fail any variants with synonyms, use INSERT IGNORE because dbSNP
    # redundantly fail both the refsnp and the synonym sometimes, and we have
    # a unique constraint on (variation_id, failed_description_id)

    $stmt = qq{
        INSERT IGNORE INTO failed_variation (variation_id, failed_description_id)
        SELECT  v.variation_id, fd.failed_description_id
        FROM    suspect s, variation v, variation_synonym vs, failed_description fd
        WHERE   fd.description = 'Flagged as suspect by dbSNP'
        AND     vs.name = CONCAT('rs', s.snp_id)
        AND     v.variation_id = vs.variation_id
    };

    $self->{'dbVar'}->do($stmt);

    #$self->{'dbVar'}->do(qq{DROP TABLE suspect});
}

# filling of the variation table from SubSNP and SNP
# creating of a link table variation_id --> subsnp_id
sub variation_table {
    my $self = shift;
  
  #ÊIf this variable is set, variations with a subsnp_id below this integer will not be imported. This is useful for resuming when resuming an import that crashed at a particular subsnp_id. Also, any SQL statements preparing for the import will be skipped.
  my $resume_at_subsnp_id = -1;
  
  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
    my $stmt;
  $self->{'dbVar'}->do( "ALTER TABLE variation add column snp_id int" ) unless ($resume_at_subsnp_id > 0);
  print $logh Progress::location();

   debug(localtime() . "\tDumping RefSNPs");
    
   $stmt = qq{
              SELECT 
                COUNT(*) 
              FROM 
                SNPAncestralAllele
             };
  my $count = $self->{'dbSNP'}->db_handle()->selectall_arrayref($stmt);
   ###if we have data in table SNPAncestralAllele, like human, will use it, otherwise goto else
   ## Note that in setting the validation_status below, the elements from the validation_status are
   ## chosen according to their bitmapped decimal values. Therefore, when updating the set, new values
   ## cannot just be appended, they have to be added in the same order as they are in the dbSNP
   ## SnpValidationCode table 
   if ($count->[0][0]) {
     $stmt = "SELECT ";
     if ($self->{'limit'}) {
       $stmt .= "TOP $self->{'limit'} ";
     }
     $stmt .= qq{
	          1, 
	          'rs'+LTRIM(STR(snp.snp_id)) AS sorting_id, 
	          CASE WHEN
	            snp.validation_status = 0
	          THEN
	            NULL
	          ELSE
	            snp.validation_status
	          END, 
	          a.allele, 
	          snp.snp_id
		FROM 
		  SNP snp LEFT JOIN 
		  SNPAncestralAllele snpa ON snp.snp_id = snpa.snp_id LEFT JOIN 
		  $self->{'dbSNP_share_db'}..Allele a on snpa.ancestral_allele_id = a.allele_id
		WHERE 
		  exemplar_subsnp_id != 0
	       };
     if ($self->{'limit'}) {
       $stmt .= qq{    
		   ORDER BY
		     sorting_id ASC  
	          };
     }
    dumpSQL($self->{'dbSNP'},$stmt) unless ($resume_at_subsnp_id > 0);
   }
   else {
     $stmt = "SELECT ";
     if ($self->{'limit'}) {
       $stmt .= "TOP $self->{'limit'} ";
     }
     $stmt .= qq{
	          1, 
	          'rs'+LTRIM(STR(snp.snp_id)) AS sorting_id, 
	          CASE WHEN
	            snp.validation_status = 0
	          THEN
	            NULL
	          ELSE
	            snp.validation_status
	          END,
	          NULL,
	          snp.snp_id
		FROM 
		  SNP snp 
		WHERE 
		  exemplar_subsnp_id != 0
	       };
     if ($self->{'limit'}) {
       $stmt .= qq{    
		   ORDER BY
		     sorting_id ASC  
	          };
     }
     dumpSQL($self->{'dbSNP'},$stmt) unless ($resume_at_subsnp_id > 0);
   }
    
    debug(localtime() . "\tLoading RefSNPs into variation table");
  load( $self->{'dbVar'}, "variation", "source_id", "name", "validation_status", "ancestral_allele", "snp_id" ) unless ($resume_at_subsnp_id > 0);
  print $logh Progress::location();
  
  $self->{'dbVar'}->do( "ALTER TABLE variation ADD INDEX snpidx( snp_id )" ) unless ($resume_at_subsnp_id > 0);
  print $logh Progress::location();

	#ÊGet somatic flag from dbSNP but only if the snp exclusively has somatic subsnps
	$stmt = qq{
		SELECT DISTINCT
			sssl.snp_id
		FROM
			SubSNP ss JOIN
			SNPSubSNPLink sssl ON (
				sssl.subsnp_id = ss.subsnp_id AND
				ss.SOMATIC_ind = 'Y'
			)
		WHERE
		  NOT EXISTS (
		      SELECT 
		          * 
		      FROM 
		          SNPSubSNPLink sssl2 JOIN 
                  SubSNP ss2 ON (
                    sssl2.snp_id = sssl.snp_id AND 
                    sssl2.subsnp_id != sssl.subsnp_id AND 
                    ss2.subsnp_id = sssl.subsnp_id
                  ) 
              WHERE 
                ss2.SOMATIC_ind != ss.SOMATIC_ind
         )
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
    
# create table rsHist to store rs history
# Table vwRsMergeArch only exists for human. Filter this as a temporary solution and talk to Yuan about it
    if ($self->{'dbm'}->dbCore()->species =~ m/human|homo/i) {
      dumpSQL($self->{'dbSNP'},(qq{SELECT * FROM RsMergeArch})) unless ($resume_at_subsnp_id > 0);
#loading it to variation database
    create_and_load( $self->{'dbVar'}, "rsHist", "rsHigh *", "rsCurrent *","orien2Current") unless ($resume_at_subsnp_id > 0);
  print $logh Progress::location();
    }
    else {
      $stmt = qq{
	CREATE TABLE
	  rsHist (
	    rsHigh VARCHAR(255) DEFAULT NULL,
	    rsCurrent VARCHAR(255) DEFAULT NULL,
	    orien2Current TINYINT(4) NOT NULL,
	    UNIQUE KEY rsHigh_2 (rsHigh),
	    KEY rsCurrent (rsCurrent),
	    KEY rsHigh (rsHigh)
	  )
      };
      $self->{'dbVar'}->do($stmt) unless ($resume_at_subsnp_id > 0);
  print $logh Progress::location();
    }
    
#change rs_id to rs_name, i.e add rs to number
  $self->{'dbVar'}->do(qq{UPDATE rsHist SET rsHigh = CONCAT('rs',rsHigh), rsCurrent = CONCAT('rs',rsCurrent)}) unless ($resume_at_subsnp_id > 0);
  print $logh Progress::location();
    
    # create a temp table of subSNP info
    
    debug(localtime() . "\tDumping SubSNPs");
    

   $stmt = "SELECT ";
   if ($self->{'limit'}) {
     $stmt .= "TOP $self->{'limit'} ";
   }
   $stmt .= qq{
                 subsnp.subsnp_id AS sorting_id, 
                 subsnplink.snp_id, 
                 b.pop_id, 
                 a.allele,
	         subsnplink.substrand_reversed_flag, 
	         b.moltype,
		 uv.allele_id
	       FROM 
	         SubSNP subsnp, 
	         SNPSubSNPLink subsnplink, 
	         $self->{'dbSNP_share_db'}..ObsVariation ov, 
	         Batch b, 
	         $self->{'dbSNP_share_db'}..UniVariAllele uv, 
	         $self->{'dbSNP_share_db'}..Allele a
	       WHERE 
	         subsnp.batch_id = b.batch_id AND   
	         subsnp.subsnp_id = subsnplink.subsnp_id AND   
	         ov.var_id = subsnp.variation_id AND   
	         ov.univar_id = uv.univar_id AND   
	         uv.allele_id=a.allele_id
	      };
    if ($self->{'limit'}) {
      $stmt .= qq{    
	          ORDER BY
	            sorting_id ASC  
	         };
    }
   dumpSQL($self->{'dbSNP'},$stmt) unless ($resume_at_subsnp_id > 0);
   create_and_load( $self->{'dbVar'}, "tmp_var_allele", "subsnp_id i*", "refsnp_id i*","pop_id i", "allele", "substrand_reversed_flag i", "moltype", "allele_id i") unless ($resume_at_subsnp_id > 0);
  print $logh Progress::location();
    
    # load the synonym table with the subsnp identifiers
    
   debug(localtime() . "\tloading variation_synonym table with subsnps");

   $self->{'dbVar'}->do(qq{ALTER TABLE variation_synonym add column substrand_reversed_flag tinyint}) unless ($resume_at_subsnp_id > 0);
  print $logh Progress::location();
    
  # The statement below cause the server to run out of space in the /tmp directory. As a workaround (and speed-up), we dump out the contents and
  # do the group operation in the shell. We'll also do a reasonable amount of rows each time to avoid overloading the system
=head The old, problematic code:
   $self->{'dbVar'}->do( qq{INSERT INTO variation_synonym (variation_id, source_id, name,
						  moltype, subsnp_id, substrand_reversed_flag )
		       SELECT v.variation_id, 1,
		       CONCAT('ss',tv.subsnp_id), tv.moltype, tv.subsnp_id,
		       tv.substrand_reversed_flag
		       FROM tmp_var_allele tv, variation v
		       WHERE tv.refsnp_id = v.snp_id
		       GROUP BY tv.subsnp_id 			
		   }); #can use distinct instead
  print $logh Progress::location();
    
    $self->{'dbVar'}->do("ALTER TABLE variation_synonym ADD INDEX subsnp_id(subsnp_id)");
  print $logh Progress::location();
=cut

  # -- BEGIN WORKAROUND CODE --
  
  # Subsnp_id interval to export at each go
  my $interval = 1e6;
  my $offset = max(0,$resume_at_subsnp_id);
  #ÊGet the minimum and maximum subsnp_id. We will export rows based on subsnp_id rather than using limit and offset in order to avoid having the same subsnp_id split across two exports
  $stmt = qq{
    SELECT
      MIN(tv.subsnp_id) AS mn,
      MAX(tv.subsnp_id) AS mx
    FROM
      tmp_var_allele tv
  };
  my ($min_id,$max_id) = @{$self->{'dbVar'}->db_handle()->selectall_arrayref($stmt)->[0]};
  print $logh Progress::location();
  
  # Write the columns to be imported to a temporary file
  my $tmpfile = $self->{'tmpdir'} . '/' . $self->{'tmpfile'} . '_var_syn_import';
  # The write mode should be '>' unless we are resuming
  my $write_mode = ($resume_at_subsnp_id > 0 ? '>>' : '>');
  open(IMP,$write_mode,$tmpfile) or die ('Could not open temporary variation_synonym import file for writing');
  
  #ÊLoop until we've dumped all rows
  while ($offset < $max_id) {
    
    debug(localtime() . "\tDumping from tmp_var_allele where subsnp_id between $offset and " . ($offset + $interval));
    
    # Dump the statement to a tempfile
    $stmt = qq{
      SELECT
	v.variation_id,
        tv.subsnp_id,
	tv.moltype,
	tv.substrand_reversed_flag
      FROM
        tmp_var_allele tv USE INDEX (subsnp_id_idx, refsnp_id_idx),
	variation v USE INDEX (snpidx)
      WHERE
        tv.subsnp_id BETWEEN $offset AND ($offset + $interval) AND
	v.snp_id = tv.refsnp_id
      ORDER BY
        tv.subsnp_id ASC,
	tv.refsnp_id ASC,
	tv.pop_id ASC
    };
    dumpSQL($self->{'dbVar'},$stmt);
    print $logh Progress::location() . "\tDumped from tmp_var_allele for subsnp_ids between $offset and " . ($offset + $interval) . "\n";
    
    # Parse the dumpfile and write to the import file but skip rows with duplicated subsnp_id (emulate the group statement)
    open(FH,'<',$self->{'tmpdir'} . "/" . $self->{'tmpfile'}) or die ('Could not open tmp_var_allele dumpfile for reading');
    my $last = -1;
    while (<FH>) {
      chomp;
      my @arr = split(/\t/);
      next if ($arr[1] == $last);
      $last = $arr[1];
      print IMP $arr[0] . "\t1\tss" . $arr[1] . "\t" . $arr[2] . "\t" . $arr[1] . "\t" . $arr[3] . "\n";
    }
    close(FH);
    print $logh Progress::location();
    
    # Increase the offset
    $offset += ($interval + 1);
  }
  close(IMP);
  
  #ÊWe'll move the temp import file to the tempfile used by the ImportUtils package
  (system(("mv",$tmpfile,$ImportUtils::TMP_DIR . "/" . $ImportUtils::TMP_FILE)) == 0) or die ("Could not rename $tmpfile to " . $ImportUtils::TMP_DIR . "/" . $ImportUtils::TMP_FILE);
  print $logh Progress::location();
  
  debug(localtime() . "\tLoading into variation_synonym table");
  # Now, load the import file into the variation_synonym table
  load(
    $self->{'dbVar'},
    "variation_synonym",
    "variation_id",
    "source_id",
    "name",
    "moltype",
    "subsnp_id",
    "substrand_reversed_flag"
  );
  print $logh Progress::location();
  
  # -- END WORKAROUND CODE -- 
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
    dumpSQL($self->{'dbSNP'},$stmt);
    load( $self->{'dbVar'}, "subsnp_handle", "subsnp_id", "handle");
  print $logh Progress::location();
  
    # import any clinical significance, global minor allele frequencies or suspect SNPs
    # these subroutines all check if the table is present and populated before doing
    # anything so should work fine on species without the necessary tables

    $self->clin_sig;
    $self->minor_allele_freq;
    $self->suspect_snps;

    return;
}


#
# dumps subSNPs and associated allele information
#
sub dump_subSNPs {
    my $self = shift;

  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
    my $stmt = "SELECT ";
    if ($self->{'limit'}) {
      $stmt .= "TOP $self->{'limit'} ";
    }
    $stmt .= qq{
	          subsnp.subsnp_id AS sorting_id, 
	          subsnplink.snp_id, 
	          b.pop_id, 
	          ov.pattern,
	          subsnplink.substrand_reversed_flag, 
	          b.moltype
                FROM 
                  SubSNP subsnp, 
                  SNPSubSNPLink subsnplink, 
                  $self->{'dbSNP_share_db'}..ObsVariation ov, 
                  Batch b
                WHERE 
                  subsnp.batch_id = b.batch_id AND
                  subsnp.subsnp_id = subsnplink.subsnp_id AND   
                  ov.var_id = subsnp.variation_id
	       };
    if ($self->{'limit'}) {
      $stmt .= qq{    
	          ORDER BY
	            sorting_id ASC  
	         };
    }

    my $sth = $self->{'dbSNP'}->prepare($stmt,{mysql_use_result => 1});
    
    $sth->execute();

    open ( FH, ">" . $self->{'tmpdir'} . "/" . $self->{'tmpfile'} );
    
  my ($row);
  while($row = $sth->fetchrow_arrayref()) {
    my $prefix;
    my @alleles = split('/', $row->[3]);
    if ($row->[3] =~ /^(\(.*\))\d+\/\d+/) {
      $prefix = $1;
    }
    my @row = map {(defined($_)) ? $_ : '\N'} @$row;

    # split alleles into multiple rows
    foreach my $a (@alleles) {
      if ($prefix and $a !~ /\(.*\)/) {#incase (CA)12/13/14 CHANGE TO (CA)12/(CA)13/(CA)14
	$a = $prefix.$a;
	#$prefix = "";
      }
      $row[3] = $a;
      print FH join("\t", @row), "\n";
    }
  }

  $sth->finish();
  close FH;
}


#
# loads the population table
#
# This subroutine produces identical results as the MySQL equivalence
#
sub population_table {
    my $self = shift;
  
  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  $self->{'dbVar'}->do("ALTER TABLE sample ADD column pop_id int");   
  print $logh Progress::location();
  $self->{'dbVar'}->do("ALTER TABLE sample ADD column pop_class_id int"); 
  print $logh Progress::location();

  # load PopClassCode data as populations

  debug(localtime() . "\tDumping population class data");

  my $stmt = qq{
                SELECT 
                  RTRIM(pop_class), 
                  RTRIM(pop_class_id), 
                  RTRIM(pop_class_text)
		FROM 
		  $self->{'dbSNP_share_db'}..PopClassCode
               };
  dumpSQL($self->{'dbSNP'},$stmt);

  load($self->{'dbVar'}, 'sample', 'name', 'pop_class_id', 'description');
  print $logh Progress::location();

  $self->{'dbVar'}->do(qq{ALTER TABLE sample ADD INDEX pop_class_id (pop_class_id)});
  print $logh Progress::location();

  debug(localtime() . "\tDumping population data");

  # load Population data as populations
  
#  For compatibility with MS SQL, moved the group operations to the MySQL database
#    $self->{'dbSNP'}->do("SET SESSION group_concat_max_len = 10000");

  $stmt = qq{
            SELECT DISTINCT 
              p.handle+':'+p.loc_pop_id,
	      p.pop_id, 
	      pc.pop_class_id, 
	      pl.line,
	      pl.line_num
	    FROM   
	      Population p LEFT JOIN 
	      $self->{'dbSNP_share_db'}..PopClass pc ON p.pop_id = pc.pop_id LEFT JOIN 
	      PopLine pl ON p.pop_id = pl.pop_id
	    ORDER BY
	      p.pop_id ASC,
	      pc.pop_class_id ASC,
	      pl.line_num ASC
            };	 #table size is small, so no need to change
    dumpSQL($self->{'dbSNP'},$stmt);

    debug(localtime() . "\tLoading sample data");

    create_and_load( $self->{'dbVar'}, "tmp_pop", "name", "pop_id i*", "pop_class_id i*", "description l", "line_num i*" );
  print $logh Progress::location();
                   
    #populate the Sample table with the populations
    $self->{'dbVar'}->do("SET SESSION group_concat_max_len = 10000");
    $self->{'dbVar'}->do(qq{INSERT INTO sample (name, pop_id,description)
                 SELECT tp.name, tp.pop_id, GROUP_CONCAT(description ORDER BY tp.pop_class_id ASC, tp.line_num ASC)
                 FROM   tmp_pop tp
                 GROUP BY tp.pop_id
                 });	#table size is small, so no need to change
  print $logh Progress::location();

    $self->{'dbVar'}->do(qq{ALTER TABLE sample ADD INDEX pop_id (pop_id)});
  print $logh Progress::location();

    #and copy the data from the sample to the Population table
    debug(localtime() . "\tLoading population table with data from Sample");

    $self->{'dbVar'}->do(qq{INSERT INTO population (sample_id)
				  SELECT sample_id
			          FROM sample});
  print $logh Progress::location();

     debug(localtime() . "\tLoading population_synonym table");

     # build super/sub population relationships
     $self->{'dbVar'}->do(qq{INSERT INTO population_structure (super_population_sample_id,sub_population_sample_id)
 				    SELECT DISTINCT p1.sample_id, p2.sample_id
 				    FROM tmp_pop tp, sample p1, sample p2
 				    WHERE tp.pop_class_id = p1.pop_class_id
 				    AND   tp.pop_id = p2.pop_id});
  print $logh Progress::location();
    

     #load population_synonym table with dbSNP population id
     $self->{'dbVar'}->do(qq{INSERT INTO sample_synonym (sample_id,source_id,name)
 				      SELECT sample_id, 1, pop_id
 				      FROM sample
 				      WHERE pop_id is NOT NULL
 				  });
  print $logh Progress::location();
    
     $self->{'dbVar'}->do("DROP TABLE tmp_pop");
  print $logh Progress::location();
}



# loads the individual table
#
# This subroutine produces identical results as the MySQL equivalence
#
sub individual_table {
    my $self = shift;

  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  # load individuals into the population table
  my $stmt;
  debug(localtime() . "\tDumping Individual data");

  # a few submitted  individuals have the same individual or no individual
  # we ignore this problem with a group by
  #there were less individuals in the individual than in the individual_genotypes table, the reason is that some individuals do not have
  #assigned a specie for some reason in the individual table, but they do have in the SubmittedIndividual table
  #to solve the problem, get the specie information from the SubmittedIndividual table
  $stmt = qq{ 
    SELECT
      CASE WHEN
	si.loc_ind_alias IS NULL OR
	LEN(si.loc_ind_alias) = 0
      THEN
	LTRIM(STR(si.pop_id))+'_'+LTRIM(si.loc_ind_id)
      ELSE
	si.loc_ind_alias
      END, 
      i.descrip, 
      i.ind_id,
      si.submitted_ind_id
    FROM   
      SubmittedIndividual si, 
      Individual i
    WHERE  
      si.ind_id = i.ind_id
  };
  dumpSQL($self->{'dbSNP'}, $stmt);

  create_and_load($self->{'dbVar'}, 'tmp_ind_ungrouped', 'loc_ind_id', 'description', 'ind_id i*', 'submitted_ind_id i*');
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{ALTER TABLE tmp_ind_ungrouped ORDER BY submitted_ind_id ASC, ind_id ASC, loc_ind_id ASC, description ASC});
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{
                          CREATE TABLE 
                            tmp_ind
                          SELECT
                            loc_ind_id,
                            description,
                            ind_id
                          FROM
                            tmp_ind_ungrouped
                          GROUP BY
                            ind_id 		
                         });	#table size is small, so no need to change group by
  $self->{'dbVar'}->do(qq{DROP TABLE tmp_ind_ungrouped});
  print $logh Progress::location();

  # load pedigree into seperate tmp table because there are no
  # indexes on it in dbsnp and it makes the left join b/w tables v. slow
  # one individual has 2 (!) pedigree rows, thus the group by

  $stmt = qq{
             SELECT 
               ind_id, 
               pa_ind_id, 
               ma_ind_id, 
               sex
             FROM 
               PedigreeIndividual
            };
  dumpSQL($self->{'dbSNP'}, $stmt);

  create_and_load($self->{'dbVar'}, 'tmp_ped_ungrouped', 'ind_id i*', 'pa_ind_id i', 'ma_ind_id i', 'sex');
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{
                          CREATE TABLE 
                            tmp_ped
                          SELECT
                            ind_id,
                            pa_ind_id,
                            ma_ind_id,
                            sex
                          FROM
                            tmp_ped_ungrouped
                          GROUP BY
                            ind_id;
                         });
  $self->{'dbVar'}->do(qq{DROP TABLE tmp_ped_ungrouped});
  print $logh Progress::location();
                         
  debug(localtime() . "\tLoading individuals into individual table");

  # to make things easier keep dbSNPs individual.ind_id as our individual_id

  #add the individual_id column in the sample table
  $self->{'dbVar'}->do("ALTER TABLE sample ADD column individual_id int, add index ind_idx (individual_id)");   
  print $logh Progress::location();

  #and the individual data in the sample table
  $self->{'dbVar'}->do(qq{INSERT INTO sample (individual_id, name, description)
				  SELECT ti.ind_id, ti.loc_ind_id, ti.description
			          FROM tmp_ind ti
			      });
  print $logh Progress::location();

  #decide which individual_type should this species bem make sure it's correct when adding new speciesz
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

  $self->{'dbVar'}->do(qq{INSERT INTO individual (sample_id, father_individual_sample_id, mother_individual_sample_id, gender, individual_type_id)
				    SELECT s.sample_id,
				    IF(tp.pa_ind_id > 0, tp.pa_ind_id, null),
				     IF(tp.ma_ind_id > 0, tp.ma_ind_id, null),
				     IF(tp.sex = 'M', 'Male',
				        IF(tp.sex = 'F', 'Female', 'Unknown')), $individual_type_id
				    FROM sample s
				      LEFT JOIN tmp_ped tp ON s.individual_id = tp.ind_id
				      WHERE s.individual_id is not null
				});
  print $logh Progress::location();

    $self->{'dbVar'}->do("DROP table tmp_ind");
  print $logh Progress::location();

    #need to convert the father_individual_id and mother_individual_id in father_sample_individual_id and mother_sample_individual_id
    $self->{'dbVar'}->do("UPDATE individual i, sample s, individual i2 set i.father_individual_sample_id = s.sample_id where i.father_individual_sample_id = s.individual_id and i2.sample_id = s.sample_id");
  print $logh Progress::location();
   $self->{'dbVar'}->do("UPDATE individual i, sample s, individual i2 set i.mother_individual_sample_id = s.sample_id where i.mother_individual_sample_id = s.individual_id and i2.sample_id = s.sample_id");
  print $logh Progress::location();

    $self->{'dbVar'}->do("DROP table tmp_ped");
  print $logh Progress::location();
    
    #necessary to fill in the individual_population table with the relation between individual and populations
    $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX uniq_idx ON individual_population (individual_sample_id, population_sample_id)});
  print $logh Progress::location();
    
    $stmt = qq{
               SELECT 
                 si.pop_id, 
                 i.ind_id
	       FROM   
	         SubmittedIndividual si, 
	         Individual i
	       WHERE  
	         si.ind_id = i.ind_id
              };
    dumpSQL($self->{'dbSNP'},$stmt);
    
    create_and_load($self->{'dbVar'}, 'tmp_ind_pop', 'pop_id i*', 'ind_id i*');
  print $logh Progress::location();

    debug(localtime() . "\tLoading individuals_population table");

    $self->{'dbVar'}->do(qq(INSERT IGNORE INTO individual_population (individual_sample_id, population_sample_id)
				      SELECT s1.sample_id, s2.sample_id
				      FROM tmp_ind_pop tip, sample s1, sample s2
				      WHERE tip.pop_id = s2.pop_id
				      AND s1.individual_id = tip.ind_id
				  ));
  print $logh Progress::location();
    $self->{'dbVar'}->do("DROP INDEX uniq_idx ON individual_population");
  print $logh Progress::location();
    $self->{'dbVar'}->do("DROP table tmp_ind_pop");
  print $logh Progress::location();

    #necessary to fill in the sample_synonym table with the relation between individual_id and sample_id
    $self->{'dbVar'}->do(qq{INSERT INTO sample_synonym (sample_id,source_id,name)
 				      SELECT sample_id, 1, individual_id
 				      FROM sample
 				      WHERE individual_id is NOT NULL
 				  });
  print $logh Progress::location();

    return;
}

sub parallelized_allele_table {
  my $self = shift;
  my $load_only = shift;
  
  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  my $stmt;
  
  # The tempfile to be used for loading
  my $file_prefix = $self->{'tmpdir'} . '/allele_table';
  my $loadfile = $file_prefix . '_loadfile.txt';
  # Sample file is used for caching the samples
  my $samplefile = $file_prefix . '_population_samples.txt';
  # Allele file is used for caching the alleles
  my $allelefile = $file_prefix . '_alleles.txt';
  
  # Do the extraction of the alleles unless we're resuming and will only load the alleles
  my $task_manager_file;
  my $jobindex;
  my $jobid;
  my $jobname;
  my $result;
  unless ($load_only) {
    #ÊFirst, get the population_id -> sample_id mapping and write it to the file. The subroutine can also get it but it's faster to get all at once since we expect many to be used. 
    $stmt = qq{
      SELECT DISTINCT
	s.pop_id,
	s.sample_id
      FROM
	sample s
      WHERE
	s.pop_id IS NOT NULL AND
	s.pop_id > 0
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
  
    #ÊProcess the alleles in chunks based on the SubSNP id. This number should be kept at a reasonable level, ideally so that we don't need to request more than 4GB memory on the farm. If so, we'll have access to the most machines.
    #ÊThe limitation is that the results from the dbSNP query is read into memory. If necessary, we can dump it to a file (if the results are sorted)
    my $chunksize = 5e5;
    
    $stmt = qq{
      SELECT
	MIN(ss.subsnp_id),
	MAX(ss.subsnp_id)
      FROM
	SubSNP ss
    };
    my ($min_ss,$max_ss) = @{$self->{'dbSNP'}->db_handle()->selectall_arrayref($stmt)->[0]};
    my $ss_offset = ($min_ss - 1);
    
    $task_manager_file = $file_prefix . '_task_management.txt';
    $jobindex = 0;
    
    #ÊLoop over all subtask intervals and write the parameters to the manager file
    open(MGMT,'>',$task_manager_file);
    while ($ss_offset < $max_ss) {
      
      $jobindex++;
      my $start = ($ss_offset + 1);
      my $end = $ss_offset + $chunksize;
      $ss_offset = $end;
      
      print MGMT qq{$jobindex $loadfile $start $end $allelefile $samplefile\n};
    }
    close(MGMT);
    
    # Run the job on the farm
    $jobname = 'allele_table';
    $result = $self->run_on_farm($jobname,$file_prefix,'allele_table',$task_manager_file,1,$jobindex);
    $jobid = $result->{'jobid'};
    
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
  
  $task_manager_file = $file_prefix . '_task_management_load.txt';
  open(MGMT,'>',$task_manager_file);
  $jobindex = 1;
  print MGMT qq{$jobindex $loadfile allele variation_id subsnp_id sample_id allele frequency count\n};
  close(MGMT);
  
  # Run the job on the farm
  print $logh Progress::location() . "\tSubmitting loading of the allele_table to the farm\n";
  $jobname = 'allele_table_load';
  $result = $self->run_on_farm($jobname,$file_prefix,'load_data_infile',$task_manager_file,1,1);
  $jobid = $result->{'jobid'};
  
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
  die("Loading of the allele_table failed (see log output). This needs to be resolved before proceeding!") unless ($result->{'success'});
  
  #ÊFinally, create the allele_string table needed for variation_feature    
  $stmt = qq{
    SELECT 
      snp.snp_id, 
      a.allele
    FROM   
      SNP snp JOIN 
      $self->{'dbSNP_share_db'}..UniVariAllele uva ON (
	uva.univar_id = snp.univar_id
      ) JOIN
      $self->{'dbSNP_share_db'}..Allele a ON (
	a.allele_id = uva.allele_id
      )
  };
  dumpSQL($self->{'dbSNP'},$stmt);
  print $logh Progress::location();
  create_and_load($self->{'dbVar'},"tmp_allele_string","snp_name *","allele");
  print $logh Progress::location();
    
  $stmt = qq{
    CREATE TABLE
      allele_string
    SELECT
      v.variation_id AS variation_id,
      tas.allele AS allele
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
}

sub allele_table {
  my $self = shift;
  
  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  my $stmt;
  
  #ÊPrepared statement for getting the SubSNPs with or without population frequency data
  $stmt = qq{
    SELECT
      ss.subsnp_id,
      CASE WHEN
	b.pop_id IS NULL
      THEN
	0
      ELSE
	b.pop_id
      END AS pop_id,
      uv.allele_id,
      sssl.substrand_reversed_flag,
      '\\N' AS frequency,
      '\\N' AS count
    FROM
      SNPSubSNPLink sssl JOIN
      SubSNP ss ON (
	ss.subsnp_id = sssl.subsnp_id
      ) JOIN
      Batch b ON (
	b.batch_id = ss.batch_id
      ) JOIN
      $self->{'dbSNP_share_db'}..ObsVariation ov ON (
	ov.var_id = ss.variation_id
      ) JOIN
      $self->{'dbSNP_share_db'}..UniVariAllele uv ON (
	uv.univar_id = ov.univar_id
      )
    WHERE
      ss.subsnp_id BETWEEN ? AND ?
    ORDER BY
      ss.subsnp_id ASC,
      uv.allele_id ASC,
      b.pop_id ASC
  };
  my $ss_sth = $self->{'dbSNP'}->prepare($stmt);
  
  # Prepared statement to get all SubSNPs that have population frequency data
  $stmt = qq{
    SELECT
      ss.subsnp_id,
      afbsp.pop_id,
      afbsp.allele_id,
      sssl.substrand_reversed_flag,
      afbsp.freq,
      afbsp.cnt
    FROM
      SNPSubSNPLink sssl JOIN
      SubSNP ss ON (
	ss.subsnp_id = sssl.subsnp_id
      ) JOIN
      AlleleFreqBySsPop afbsp ON (
	afbsp.subsnp_id = ss.subsnp_id
      )
    WHERE
      ss.subsnp_id BETWEEN ? AND ?
    ORDER BY
      ss.subsnp_id ASC,
      afbsp.allele_id ASC,
      afbsp.pop_id ASC
  };
  my $ss_freq_sth = $self->{'dbSNP'}->prepare($stmt);
  
  #ÊPrepared statement to get the corresponding variation_ids for a range of subsnps from variation_synonym
  $stmt = qq{
    SELECT
      vs.subsnp_id,
      vs.variation_id
    FROM
      variation_synonym vs
    WHERE
      vs.subsnp_id BETWEEN ? AND ?
  };
  my $vs_sth = $self->{'dbVar'}->prepare($stmt);
  
  # Prepared statement to get the alleles
  $stmt = qq{
    SELECT
      a.allele,
      ar.allele
    FROM
      $self->{'dbSNP_share_db'}..Allele a JOIN
      $self->{'dbSNP_share_db'}..Allele ar ON (
	ar.allele_id = a.rev_allele_id
      )
    WHERE
      a.allele_id = ?
  };
  my $allele_sth = $self->{'dbSNP'}->prepare($stmt);
  
  # Prepared statement to get the sample_id from a pop_id
  $stmt = qq{
    SELECT
      s.sample_id
    FROM
      sample s
    WHERE
      s.pop_id = ?
    LIMIT 1
  };
  my $sample_sth = $self->{'dbVar'}->prepare($stmt);
  
  #ÊProcess the alleles in chunks based on the SubSNP id
  my $chunksize = 1e6;
  
  $stmt = qq{
    SELECT
      MIN(ss.subsnp_id),
      MAX(ss.subsnp_id)
    FROM
      SubSNP ss
  };
  my ($min_ss,$max_ss) = @{$self->{'dbSNP'}->db_handle()->selectall_arrayref($stmt)->[0]};
  
  #ÊHash to hold the alleles in memory
  my %alleles;
  #ÊHash to keep sample_ids in memory
  my %samples;
  # The pop_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value
  $samples{0} = '\N';
  
  # The tempfile to be used for loading
  my $tempfile = $self->{'tmpdir'} . '/' . $self->{'tmpfile'};
  # Open a file handle to the temp file that will be used for loading
  open(IMP,'>',$tempfile) or die ("Could not open import file $tempfile for writing");
  
  my $ss_offset = ($min_ss - 1);
  
  while ($ss_offset < $max_ss) {
    
    print $logh Progress::location() . "\tProcessing allele data for SubSNP ids " . ($ss_offset + 1) . "-" . ($ss_offset + $chunksize) . "\n";
    
    # First, get all SubSNPs
    $ss_sth->bind_param(1,($ss_offset + 1),SQL_INTEGER);
    $ss_freq_sth->bind_param(1,($ss_offset + 1),SQL_INTEGER);
    $vs_sth->bind_param(1,($ss_offset + 1),SQL_INTEGER);
    $ss_offset += $chunksize;
    $ss_sth->bind_param(2,$ss_offset,SQL_INTEGER);
    $ss_freq_sth->bind_param(2,$ss_offset,SQL_INTEGER);
    $vs_sth->bind_param(2,$ss_offset,SQL_INTEGER);
    
    #ÊFetch the query result as a hashref
    $ss_sth->execute();
    my $subsnp_alleles = $ss_sth->fetchall_hashref(['subsnp_id','pop_id','allele_id']);
    print $logh Progress::location();
    
    # Fetch the population frequency alleles as an arrayref
    $ss_freq_sth->execute();
    my $subsnp_freq_alleles = $ss_freq_sth->fetchall_arrayref();
    print $logh Progress::location();
    
    # Add new alleles with frequency data and update the ones with data missing in the first hash (where subsnp_id, pop_id and allele_id match).
    #ÊIs this potentially VERY memory intensive?
    map {
    $subsnp_alleles->{$_->[0]}{$_->[1]}{$_->[2]} = {'substrand_reversed_flag' => $_->[3], 'frequency' => $_->[4], 'count' => $_->[5]};
    } @{$subsnp_freq_alleles};
    print $logh Progress::location();
    
    # Fetch the subsnp_id to variation_id mapping and store as hashref
    $vs_sth->execute();
    my $variation_ids = $vs_sth->fetchall_hashref(['subsnp_id']);
    print $logh Progress::location();
    
    #ÊNow, loop over the subsnp hash and print it to the tempfile so we can import the data. Replace the allele_id with the corresponding allele on-the-fly
    while (my ($subsnp_id,$pop_hash) = each(%{$subsnp_alleles})) {
      while (my ($pop_id,$allele_hash) = each(%{$pop_hash})) {
	while (my ($allele_id,$allele_data) = each(%{$allele_hash})) {
	  # Look up the allele in the database if necessary
	  if (!exists($alleles{$allele_id})) {
	    $allele_sth->execute($allele_id);
	    my ($a,$arev);
	    $allele_sth->bind_columns(\$a,\$arev);
	    $allele_sth->fetch();
	    $alleles{$allele_id} = [$a,$arev];
	  }
	  #ÊLook up the sample id in the database if necessary
	  if (!exists($samples{$pop_id})) {
	    $sample_sth->execute($pop_id);
	    my $sample_id;
	    $sample_sth->bind_columns(\$sample_id);
	    $sample_sth->fetch();
	    $samples{$pop_id} = $sample_id;
	  }
	  print IMP join("\t",($variation_ids->{$subsnp_id}{'variation_id'},$subsnp_id,$samples{$pop_id},$alleles{$allele_id}->[$allele_data->{'substrand_reversed_flag'}],$allele_data->{'frequency'},$allele_data->{'count'})) . "\n";
	}
      }
    }
    print $logh Progress::location();
  }
  close(IMP);
  
  #ÊDisable the keys on allele before doing an insert in order to speed up the insert
  $stmt = qq {
    ALTER TABLE
      allele
    DISABLE KEYS
  };
  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();
  
  # Load the data from the infile
  load($self->{'dbVar'},"allele","variation_id","subsnp_id","sample_id","allele","frequency","count");
  print $logh Progress::location();
  
  # Enable the keys on allele after the inserts
  $stmt = qq {
    ALTER TABLE
      allele
    ENABLE KEYS
  };
  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();
  
  #ÊFinally, create the allele_string table needed for variation_feature    
  $stmt = qq{
    SELECT 
      snp.snp_id, 
      a.allele
    FROM   
      SNP snp JOIN 
      $self->{'dbSNP_share_db'}..UniVariAllele uva ON (
	uva.univar_id = snp.univar_id
      ) JOIN
      $self->{'dbSNP_share_db'}..Allele a ON (
	a.allele_id = uva.allele_id
      )
  };
  dumpSQL($self->{'dbSNP'},$stmt);
  create_and_load($self->{'dbVar'},"tmp_allele_string","snp_name *","allele");
  print $logh Progress::location();
    
  $stmt = qq{
    CREATE TABLE
      allele_string
    SELECT
      v.variation_id AS variation_id,
      tas.allele AS allele
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
}

#
# loads the allele table
#
sub old_allele_table {
    my $self = shift;

  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
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
                   $self->{'dbSNP_share_db'}..Allele a1, 
                   $self->{'dbSNP_share_db'}..Allele a2
                 WHERE 
                   a1.rev_allele_id = a2.allele_id
                };
      dumpSQL($self->{'dbSNP'},$stmt);
    
      create_and_load($self->{'dbVar'}, "tmp_rev_allele", "allele_id i*","allele *", "rev_allele");
  print $logh Progress::location();
    }
        
    #first load the allele data for alleles that we have population and
    #frequency data for

    $stmt = "SELECT ";
    if ($self->{'limit'}) {
      $stmt .= "TOP $self->{'limit'} ";
    }
    $stmt .= qq{
                  afsp.subsnp_id AS sorting_id, 
                  afsp.pop_id, 
                  afsp.allele_id, 
                  afsp.freq,
		  afsp.cnt
		FROM   
		  AlleleFreqBySsPop afsp, 
		  SubSNP ss
		WHERE    
		  afsp.subsnp_id = ss.subsnp_id
	       };
     if ($self->{'limit'}) {
       $stmt .= qq{    
		   ORDER BY
		     sorting_id ASC  
	          };
     }
    dumpSQL($self->{'dbSNP'},$stmt);

   debug(localtime() . "\tLoading allele frequency data");

   create_and_load($self->{'dbVar'}, "tmp_allele", "subsnp_id i*", "pop_id i*","allele_id i*", "freq","cnt i");
  print $logh Progress::location();

    debug(localtime() . "\tCreating allele table");
  
  #necessary to create a unique index to simulate the GROUP BY clause
  $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX unique_allele_idx ON allele (variation_id,allele(2),frequency,sample_id)});
  print $logh Progress::location();

  # PL: This is a slow query (~5h), should we try to parallelize it? Perhaps split it in chunks by subsnp_id? Is there any point in doing that when all threads still has to access the same db/table?
   $self->{'dbVar'}->do(qq(INSERT IGNORE INTO allele (variation_id, subsnp_id,allele,frequency, sample_id, count)
               SELECT vs.variation_id,vs.subsnp_id,
                      IF(vs.substrand_reversed_flag,
                         tra.rev_allele,tra.allele) as allele,
                      ta.freq, s.sample_id, ta.cnt
               FROM   tmp_allele ta, tmp_rev_allele tra, variation_synonym vs,
                      sample s
               WHERE  ta.subsnp_id = vs.subsnp_id
	       AND    ta.allele_id = tra.allele_id
	       AND    ta.pop_id = s.pop_id),{mysql_use_result => 1} );
  print $logh Progress::location();

   $self->{'dbVar'}->do("ALTER TABLE allele ENABLE KEYS"); #after ignoring in the insertion, we must enable keys again
  print $logh Progress::location();


   $self->{'dbVar'}->do("DROP TABLE tmp_allele");    
  print $logh Progress::location();
   #going to add the other allele for the variations with 1 allele (have frequency 1 but no frequency for the other allele)
   debug(localtime() . "\tLoading allele data without frequency");
   
   $self->{'dbVar'}->do("CREATE TABLE tmp_allele (variation_id int, subsnp_id int, allele text, primary key (variation_id,allele(10)))");
  print $logh Progress::location();
   $self->{'dbVar'}->do("INSERT IGNORE INTO tmp_allele SELECT variation_id, subsnp_id, allele FROM allele");
  print $logh Progress::location();

   $self->{'dbVar'}->do(qq{CREATE TABLE tmp_unique_allele 
				      SELECT ta.variation_id,  ta.subsnp_id, ta.allele, substring(vs.name,3) as snp_id
				         FROM variation vs,
				                    (SELECT variation_id, subsnp_id, allele 
						     FROM tmp_allele
						     GROUP BY variation_id 
						     HAVING COUNT(*) = 1) as ta
					 WHERE ta.variation_id = vs.variation_id});
  print $logh Progress::location();

   $self->{'dbVar'}->do("CREATE INDEX tmp_unique_allele_idx on tmp_unique_allele (variation_id)");
  print $logh Progress::location();
   $self->{'dbVar'}->do("DROP TABLE tmp_allele");
  print $logh Progress::location();
   #create table with unique alleles from dbSNP very slow with unique index??????
   $self->{'dbVar'}->do("CREATE TABLE tmp_allele (refsnp_id int, subsnp_id int, allele text, primary key (refsnp_id,allele(10)))");
  print $logh Progress::location();

  # PL: Another slow query (~5h), parallelize?
   $self->{'dbVar'}->do(qq{INSERT IGNORE INTO tmp_allele
				      SELECT tva.refsnp_id, tua.subsnp_id, IF (tva.substrand_reversed_flag, tra.rev_allele,tra.allele) as allele
				      FROM tmp_var_allele tva, tmp_rev_allele tra, tmp_unique_allele tua
				      WHERE tva.allele_id = tra.allele_id
				      AND tua.snp_id = tva.refsnp_id
				  });
  print $logh Progress::location();

   $self->{'dbVar'}->do(qq{INSERT IGNORE INTO allele (variation_id, subsnp_id,allele, frequency)
						  SELECT tua.variation_id, ta.subsnp_id, ta.allele,0
						  FROM tmp_unique_allele tua, tmp_allele ta
						  WHERE tua.snp_id = ta.refsnp_id
						  AND tua.allele <> ta.allele
						  });
  print $logh Progress::location();
   #remove the index
   $self->{'dbVar'}->do("DROP INDEX unique_allele_idx ON allele");
  print $logh Progress::location();
   $self->{'dbVar'}->do("DROP TABLE tmp_unique_allele");
  print $logh Progress::location();
   $self->{'dbVar'}->do("DROP TABLE tmp_allele");
  print $logh Progress::location();


   # load remaining allele data which we do not have frequency data for
   # this will not import alleles without frequency for a variation which
   # already has frequency
   
   debug(localtime() . "\tLoading other allele data");
  
  #ÊPL: This query is really slow (~12h), I split it up a bit to speed up but perhaps parallelize as well?
  # PL: First load everything into a new tmp_allele table
  $stmt = qq{
    CREATE TABLE
      tmp_allele
    SELECT
      vs.variation_id AS variation_id,
      vs.subsnp_id AS subsnp_id,
      tva.pop_id,
      IF(
	vs.substrand_reversed_flag,
        tra.rev_allele,
	tra.allele
      ) AS allele
    FROM
      variation_synonym vs,
      tmp_var_allele tva,
      tmp_rev_allele tra
    WHERE
      tva.subsnp_id = vs.subsnp_id AND
      tva.allele_id = tra.allele_id
  };
  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();
  
  # PL: Second, delete alleles that have already been imported into the allele table
  $stmt = qq{
    DELETE FROM
      ta
    USING
      tmp_allele ta JOIN
      allele a ON (
	ta.variation_id = a.variation_id AND
	ta.allele = a.allele
  };
  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();
  
  $self->{'dbVar'}->do("ALTER TABLE tmp_allele ADD INDEX pop_id(pop_id)");
  print $logh Progress::location();

   $self->{'dbVar'}->do(qq{INSERT INTO allele (variation_id, subsnp_id, allele,
                                      frequency, sample_id)
                  SELECT ta.variation_id, ta.subsnp_id, ta.allele, null, s.sample_id
                  FROM   tmp_allele ta
                  LEFT JOIN sample s ON s.pop_id = ta.pop_id
                  GROUP BY ta.variation_id, s.sample_id, ta.allele });
  print $logh Progress::location();

    ###generate allele_string table for variation_feature allele_string column

    debug(localtime() . "\tGenerating tmp_allele_string...");
    
    $stmt = qq{
               SELECT 
                 snp.snp_id, 
                 a.allele
               FROM   
                 SNP snp, 
                 $self->{'dbSNP_share_db'}..UniVariAllele uva, 
                 $self->{'dbSNP_share_db'}..Allele a
               WHERE  
                 snp.univar_id = uva.univar_id AND    
                 uva.allele_id = a.allele_id 
              };
    dumpSQL($self->{'dbSNP'},$stmt);

    create_and_load($self->{'dbVar'},"tmp_allele_string","snp_name *","allele");
  print $logh Progress::location();
    
    $self->{'dbVar'}->do(qq{CREATE TABLE allele_string select v.variation_id as variation_id, tas.allele as allele                                  FROM variation v, tmp_allele_string tas
                                  WHERE v.name = concat("rs",tas.snp_name)});
  print $logh Progress::location();

    $self->{'dbVar'}->do(qq{ALTER TABLE allele_string ADD INDEX variation_idx(variation_id)});
  print $logh Progress::location();


     #$self->{'dbVar'}->do("DROP TABLE tmp_allele_string");
     #$self->{'dbVar'}->do("DROP TABLE tmp_var_allele");
     #$self->{'dbVar'}->do("DROP TABLE tmp_rev_allele");
     #$self->{'dbVar'}->do("DROP TABLE tmp_allele");
}


#
# loads the flanking sequence table
#
sub flanking_sequence_table {
  my $self = shift;

  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  my $stmt;
  $self->{'dbVar'}->do(qq{CREATE TABLE tmp_seq (variation_id int,
						      subsnp_id int,
						      seq_type int,
						      line_num int,
						      type enum ('5','3'),
						      line varchar(255),
						      revcom tinyint)
				MAX_ROWS = 100000000});
  print $logh Progress::location();

# import both the 5prime and 3prime flanking sequence tables
## in human the flanking sequence tables have been partitioned
  if($self->{'dbm'}->dbCore()->species =~ /human|homo/i) {
	foreach my $type ('3','5') {
  
	  foreach my $partition('p1_human','p2_human','p3_human','ins') {
  
		debug("Dumping $type\_$partition flanking sequence");
  
		$stmt = "SELECT ";
		if ($self->{'limit'}) {
		  $stmt .= "TOP $self->{'limit'} ";
		}
		$stmt .= qq{
		    seq.subsnp_id AS sorting_id,
		    seq.type,
		    seq.line_num,
		    seq.line
		  FROM
		    SubSNPSeq$type\_$partition seq,
		    SNP snp
		  WHERE
		    snp.exemplar_subsnp_id = seq.subsnp_id
		};
		if ($self->{'limit'}) {
		  $stmt .= qq{    
			      ORDER BY
				sorting_id ASC  
			     };
		}
		
		dumpSQL($self->{'dbSNP'}, $stmt);
  
  
		$self->{'dbVar'}->do(qq{CREATE TABLE tmp_seq_$type\_$partition (
									  subsnp_id int,
									  seq_type int,
									  line_num int,
									  line varchar(255),
									  KEY subsnp_id_idx(subsnp_id))
					  MAX_ROWS = 100000000 });
  print $logh Progress::location();
  
		load($self->{'dbVar'}, "tmp_seq_$type\_$partition", "subsnp_id", "seq_type", "line_num", "line");
  print $logh Progress::location();
  
		# merge the tables into a single tmp table
		$self->{'dbVar'}->do(qq{INSERT INTO tmp_seq (variation_id, subsnp_id, seq_type,
								   line_num, type, line, revcom)
					  SELECT vs.variation_id, ts.subsnp_id, ts.seq_type, ts.line_num, '$type',
					  ts.line, vs.substrand_reversed_flag
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
		    seq.line
		  FROM 
		    SubSNPSeq$type seq, 
		    SNP snp
		  WHERE 
		    snp.exemplar_subsnp_id = seq.subsnp_id
		 };
       if ($self->{'limit'}) {
	 $stmt .= qq{    
		     ORDER BY
		       sorting_id ASC  
		    };
       }
      dumpSQL($self->{'dbSNP'},$stmt);
      
  
      $self->{'dbVar'}->do(qq{CREATE TABLE tmp_seq_$type (
								subsnp_id int,
								seq_type int,
								line_num int,
								line varchar(255),
								KEY subsnp_id_idx(subsnp_id))
				    MAX_ROWS = 100000000 });
  print $logh Progress::location();
      
      load($self->{'dbVar'}, "tmp_seq_$type", "subsnp_id", "seq_type", "line_num", "line");
  print $logh Progress::location();
  
      # merge the tables into a single tmp table
      $self->{'dbVar'}->do(qq{INSERT INTO tmp_seq (variation_id, subsnp_id, seq_type,
							 line_num, type, line, revcom)
				    SELECT vs.variation_id, ts.subsnp_id, ts.seq_type, ts.line_num, '$type',
				    ts.line, vs.substrand_reversed_flag
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
  $self->{'dbVar'}->do("DROP TABLE tmp_seq");
  print $logh Progress::location();

  debug(localtime() . "\tLoading flanking sequence data");

  # import the generated data
  load($self->{'dbVar'},"flanking_sequence","variation_id","up_seq","down_seq");
  print $logh Progress::location();

  unlink($self->{'tmpdir'} . "/" . $self->{'tmpfile'});

  return;
}



sub variation_feature {
  my $self = shift;

  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  debug(localtime() . "\tDumping seq_region data");

  dumpSQL($self->{'dbCore'}->dbc()->db_handle, qq{SELECT sr.seq_region_id, sr.name
 				      FROM   seq_region sr, coord_system cs
                                      WHERE  sr.coord_system_id = cs.coord_system_id
                                      AND    cs.attrib like "%default_version%"});
    
  debug(localtime() . "\tLoading seq_region data");
  load($self->{'dbVar'}, "seq_region", "seq_region_id", "name");
  print $logh Progress::location();

  debug(localtime() . "\tDumping SNPLoc data");

  my ($tablename1,$tablename2,$row);

  my ($assembly_version) =  $self->{'assembly_version'} =~ /^[a-zA-Z]+(\d+)\.*.*$/;
# override for platypus
  $assembly_version = 1 if $self->{'dbm'}->dbCore()->species =~ /ornith/i;
  $assembly_version = 3 if $self->{'dbm'}->dbCore()->species =~ /rerio/i;

  print $logh "assembly_version again is $assembly_version\n";
  
  my $stmt = qq{
                SELECT 
                  name 
                FROM 
                  $self->{'snp_dbname'}..sysobjects 
                WHERE 
                  name LIKE '$self->{'dbSNP_version'}\_SNPContigLoc\_$assembly_version\__'
             };
  my $sth = $self->{'dbSNP'}->prepare($stmt);
  $sth->execute();

  while($row = $sth->fetchrow_arrayref()) {
    $tablename1 = $row->[0];
  }
  
  $stmt = qq{
             SELECT 
               name 
             FROM 
               $self->{'snp_dbname'}..sysobjects 
             WHERE 
               name LIKE '$self->{'dbSNP_version'}\_ContigInfo\_$assembly_version\__'
          };
  my $sth1 = $self->{'dbSNP'}->prepare($stmt);
  $sth1->execute();

  while($row = $sth1->fetchrow_arrayref()) {
    $tablename2 = $row->[0];
  }
  debug(localtime() . "\ttable_name1 is $tablename1 table_name2 is $tablename2");
  #note the contig based cordinate is 0 based, ie. start at 0, lc_ngbr+2, t1.rc_ngbr
  
  $stmt = "SELECT ";
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
                END
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
      dumpSQL($self->{'dbSNP'},$stmt);
    
    
  debug(localtime() . "\tLoading SNPLoc data");

  create_and_load($self->{'dbVar'}, "tmp_contig_loc", "snp_id i*", "contig *", "start i", 
		  "end i", "strand i");
  print $logh Progress::location();
    
  debug(localtime() . "\tCreating genotyped variations");
    
  #creating the temporary table with the genotyped variations
  $self->{'dbVar'}->do(qq{CREATE TABLE tmp_genotyped_var SELECT DISTINCT variation_id FROM tmp_individual_genotype_single_bp});
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX variation_idx ON tmp_genotyped_var (variation_id)});
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{INSERT IGNORE INTO  tmp_genotyped_var SELECT DISTINCT variation_id FROM individual_genotype_multiple_bp});
  print $logh Progress::location();

    
  debug(localtime() . "\tCreating tmp_variation_feature data");
    
  dumpSQL($self->{'dbVar'},qq{SELECT v.variation_id, ts.seq_region_id, tcl.start, tcl.end, tcl.strand, v.name, v.source_id, v.validation_status
					  FROM variation v, tmp_contig_loc tcl, seq_region ts
					  WHERE v.snp_id = tcl.snp_id
					  AND ts.name = tcl.contig});
    
  create_and_load($self->{'dbVar'},'tmp_variation_feature',"variation_id *","seq_region_id i", "seq_region_start i", "seq_region_end i", "seq_region_strand i", "variation_name", "source_id i", "validation_status i");
  print $logh Progress::location();
    
  debug(localtime() . "\tDumping data into variation_feature table");
  $self->{'dbVar'}->do(qq{INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,variation_name, flags, source_id, validation_status)
				SELECT tvf.variation_id, tvf.seq_region_id, tvf.seq_region_start, tvf.seq_region_end, tvf.seq_region_strand,tvf.variation_name,IF(tgv.variation_id,'genotyped',NULL), tvf.source_id, tvf.validation_status
				FROM tmp_variation_feature tvf LEFT JOIN tmp_genotyped_var tgv ON tvf.variation_id = tgv.variation_id
				 });
  print $logh Progress::location();
    
    #$self->{'dbVar'}->do("DROP TABLE tmp_contig_loc");
    #$self->{'dbVar'}->do("DROP TABLE seq_region");
    #$self->{'dbVar'}->do("DROP TABLE tmp_genotyped_var");
    #$self->{'dbVar'}->do("DROP TABLE tmp_variation_feature");
}


#ÊThe task of getting the genotype will be chunked up and distributed to the farm in small pieces. Each result will be written to a loadfile that in the end
# will be used to populate the tmp_individual.. tables after everything has finished. We will chunk up the results by a) chromosome and b) individual. Because
#Êthe number of genotypes per individual varies _drastically_ with 1000 genomes data having lots of genotypes and most other data only having a few, we need
#Êto determine the best way to chunk up the individuals in order to get an even distribution of the workload
sub parallelized_individual_genotypes {
  my $self = shift;
  my $load_only = shift;
  
  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  my $genotype_table = 'tmp_individual_genotype_single_bp';
  my $multi_bp_gty_table = 'individual_genotype_multiple_bp';
  
  #ÊGet the create statement for tmp_individual_genotype_single_bp from master schema. We will need this to create the individual chromosome tables
  my $ind_gty_stmt = get_create_statement($genotype_table,$self->{'schema_file'});
  print $logh Progress::location();
  
  # Get the SubInd tables that may be split by chromosomes
  my $stmt = qq{
    SELECT 
      name 
    FROM 
      $self->{'snp_dbname'}..sysobjects 
    WHERE 
      name LIKE 'SubInd%'
  };
  my @subind_tables = map {$_->[0]} @{$self->{'dbSNP'}->db_handle()->selectall_arrayref($stmt)};

  # Prepared statement to find out if a table exists
  $stmt = qq{
    SHOW TABLES
      LIKE ?
  };
  my $table_sth = $self->{'dbVar'}->prepare($stmt);
  
  # Use one common file for alleles and one for samples.
  my $file_prefix = $self->{'tmpdir'} . '/individual_genotypes';
  #ÊMulti-bp genotypes will be written to a separate loadfile
  my $multi_bp_gty_file = $file_prefix . '_multi_bp_gty';
  # Allele file is used for caching the alleles
  my $allele_file = $file_prefix . '_alleles.txt';
  # Sample file is used for caching the samples
  my $sample_file = $file_prefix . '_individual_samples.txt';
  
  # For each SubInd_ch.. table, determine the best chunking of the individuals. We aim to get $target_rows number of rows to work on for each subtask but this will just be an approximate number.
  my $target_rows = 10e6;
  # Store the data for the farm submission in a hash with the data as key. This way, the jobs will be "randomized" w.r.t. chromosomes and chunks when submitted so all processes shouldn't work on the same tables/files at the same time.
  my %job_data;
  my $genotype_counts;
  my %gty_tables;
  my @skip_tables;
  my $sth;
  print $logh Progress::location() . "\tDividing the import task into suitable chunks\n";
  foreach my $subind_table (@subind_tables) {
    
    print $logh Progress::location() . "\t\tProcessing $subind_table\n";
    #ÊThe subtable to store the data for this subind table and the loadfile to use. The mapping file is used to temporarily store the subsnp_id -> variation_id mapping
    my $dst_table = "tmp_individual_genotype_single_bp\_$subind_table";
    my $loadfile = $file_prefix . '_' . $dst_table . '.txt';
    my $mapping_file = $file_prefix . '_subsnp_mapping_' . $subind_table . '.txt';
    $gty_tables{$subind_table} = [$dst_table,$loadfile,$mapping_file];
    
    # If the subtable already exists, warn about this but skip the iteration (perhaps we are resuming after a crash)
    $table_sth->execute($dst_table);
    print $logh Progress::location();
    if (defined($table_sth->fetchrow_arrayref())) {
      warn("Table $dst_table already exists in destination database. Will skip importing from $subind_table");
      push(@skip_tables,$subind_table);
      next;
    }
    
    #ÊIf we are skipping ahead directly to the loading of the genotypes, skip to the next iteration here
    next if ($load_only);
    
    $stmt = qq{
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
    $sth = $self->{'dbSNP'}->prepare($stmt);
    $sth->execute();
    print $logh Progress::location();
    my ($submitted_ind_id,$genotype_count);
    $sth->bind_columns(\$submitted_ind_id,\$genotype_count);
    
    #ÊLoop over the counts and when a large enough number of genotypes have been counted, split that off as a chunk to be submitted to the farm
    my $total_count = 0;
    my $start_id = -1;
    while ($sth->fetch()) {
      $total_count += $genotype_count;
      # Set the start id if it's not specified
      $start_id = $submitted_ind_id if ($start_id < 0);
      #ÊBreak off the chunk if it's large enough
      if ($total_count >= $target_rows) {
	$job_data{qq{$subind_table $start_id $submitted_ind_id}}++;
	$total_count = 0;
	$start_id = -1;
      }
    }
    #ÊAdd the rest of the individual_ids if a new chunk has been started
    $job_data{qq{$subind_table $start_id $submitted_ind_id}}++ if ($start_id >= 0);
    print $logh Progress::location();
  }
  
  #ÊIf we are only doing the load (when resuming after a crash), do not run any job on the farm
  unless ($load_only) {
  
    # If there are no genotypes at all to be imported, just return here
    return unless (scalar(keys(%job_data)));
    
    # Print the job parameters to a file
    my $task_manager_file = $file_prefix . '_task_management.txt';
    my $jobindex = 0;
    open(MGMT,'>',$task_manager_file);
    foreach my $params (keys(%job_data)) {
      $jobindex++;
      my @arr = split(/\s+/,$params);
      print MGMT qq{$jobindex $arr[0] $gty_tables{$arr[0]}->[1] $multi_bp_gty_file $arr[1] $arr[2] $gty_tables{$arr[0]}->[2]\n};
    }
    close(MGMT);
    
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
  
  # Loop over the subtables and load each of them
  my @subtables;
  my $jobindex = 0;
  my $task_manager_file = $file_prefix . '_task_management_load.txt';
  open(MGMT,'>',$task_manager_file);
  foreach my $subind_table (keys(%gty_tables)) {
    
    my $dst_table = $gty_tables{$subind_table}->[0];
    my $loadfile = $gty_tables{$subind_table}->[1];
    
    # Skip if necessary
    next if (grep(/^$subind_table$/,@skip_tables));
    
    # Create the sub table
    $stmt = $ind_gty_stmt;
    $stmt =~ s/$genotype_table/$dst_table/;
    $self->{'dbVar'}->do($stmt);
    print $logh Progress::location();
    
    # If the loadfile doesn't exist, we can't load anything
    next unless (-e $loadfile);
    
    $jobindex++;
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
  
  #ÊDrop the tmp_individual.. table if it exists
  $stmt = qq{
    DROP TABLE IF EXISTS
      $genotype_table
  };
  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();
  
  # Merge all the subtables into the final table, or if there is just one subtable, rename it to the final table name
  print $logh Progress::location() . "\tMerging the genotype subtables into a big $genotype_table table\n";
  my $merge_subtables = join(",",map {$gty_tables{$_}->[0]} keys(%gty_tables));
  if (scalar(keys(%gty_tables)) > 1) {
    
    #ÊAdd an empty table where any subsequent inserts will end up
    my $extra_table = $genotype_table . '_extra';
    $stmt = $ind_gty_stmt;
    $stmt =~ s/$genotype_table/$extra_table/;
    $self->{'dbVar'}->do($stmt);
    print $logh Progress::location();
    $merge_subtables .= ",$extra_table";
    
    $stmt = $ind_gty_stmt;
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
  
  #ÊMove multiple bp genotypes into the multiple bp table
  $stmt = qq{
      SELECT DISTINCT
          variation_id
      FROM
          $multi_bp_gty_table
  };
  dumpSQL($self->{'dbVar'},$stmt);
  print $logh Progress::location();
  
  create_and_load($self->{'dbVar'},"tmp_multiple_bp_gty_variations","variation_id");
  print $logh Progress::location();
  
  $stmt = qq{
      INSERT INTO
          $multi_bp_gty_table (
              variation_id,
              subsnp_id,
              sample_id,
              allele_1,
              allele_2
          )
      SELECT
          s.variation_id,
          s.subsnp_id,
          s.sample_id,
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
  
}

#
# loads individual genotypes into the individual_genotype table,need variation_synonym and sample tables
#
sub old_individual_genotypes {
   my $self = shift;
   
  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
   #
   # load SubInd individual genotypes into genotype table
   #

   my $TMP_DIR = $self->{'tmpdir'};
   my $TMP_FILE = $self->{'tmpfile'};
   my @gty_tables;
   my @tmp2_gty_tables;
   my @subind_tables;
   my %rec_pid1;
   my $num_process = 3;

   debug(localtime() . "\tDumping SubInd and ObsGenotype data");

   debug(localtime() . "\tTime starting to dump subind_tmp_gty table: ",scalar(localtime(time)));
   
   my $stmt = qq{
                 SELECT 
                   name 
                 FROM 
                   $self->{'snp_dbname'}..sysobjects 
                 WHERE 
                   name LIKE 'SubInd_ch%'
                };
   @subind_tables = map{$_->[0]} @{$self->{'dbSNP'}->db_handle()->selectall_arrayref($stmt)};

   if ($self->{'dbm'}->dbCore()->species !~ /homo|hum/i) {
     $num_process =1;
     @subind_tables = ("SubInd");
   }
   #@subind_tables = ("SubInd_ch3","SubInd_ch5");##only do ch3 and ch5
   foreach my $subind (@subind_tables) {
     debug(localtime() . "\tStarting with dumping $subind...");
     my $kid;
     my $pid = fork;
     if (! defined $pid){
       throw("Not possible to fork: $!\n");
     }
     elsif ($pid == 0){
       #you are the child.....
       my $dbh_child = $self->{'dbSNP'}->clone;
       $self->{'dbSNP'}->{'InactiveDestroy'} = 1;
       undef $self->{'dbSNP'};

       $ImportUtils::TMP_FILE = $subind;

       debug(localtime() . "\tTime starting to dump $subind\_tmp_gty table: ",scalar(localtime(time)));

# Watch out with this statement, if the obs column in ObsGenotype is not properly formatted (with a '/' separating alleles), the entire contents of the column will be returned (similar to the MySQL behaviour).
# The column obs_upp_fix seems to hold corrected values for most of these cases
       $stmt = "SELECT ";
       if ($self->{'limit'}) {
         $stmt .= "TOP $self->{'limit'} ";
       }
       $stmt .= qq{
                     si.subsnp_id AS sorting_id, 
                     sind.ind_id, 
                     LEN(og.obs_upp_fix) AS length_pat,
                     CASE WHEN
                       CHARINDEX('/',og.obs_upp_fix) > 0
                     THEN
                       (SELECT LEFT(og.obs_upp_fix,CHARINDEX('/',og.obs_upp_fix)-1) AS allele_1)
                     ELSE
                       (SELECT og.obs_upp_fix AS allele_1)
                     END,
                     CASE WHEN
                       CHARINDEX('/',og.obs_upp_fix) > 0
                     THEN
                       (SELECT RIGHT(og.obs_upp_fix,LEN(og.obs_upp_fix)-CHARINDEX('/',og.obs_upp_fix)) AS allele_2)
                     ELSE
                       (SELECT og.obs_upp_fix AS allele_2)
                     END,
		     CASE WHEN
		      si.submitted_strand_code IS NOT NULL
		     THEN
		      (SELECT si.submitted_strand_code)
		     ELSE
		      (SELECT 0)
		     END
		   FROM   
		     $subind si, 
		     $self->{'dbSNP_share_db'}..ObsGenotype og, 
		     SubmittedIndividual sind
		   WHERE  
		     og.gty_id = si.gty_id AND    
		     sind.submitted_ind_id = si.submitted_ind_id AND    
		     og.obs_upp_fix != 'N/N'
	          };
     if ($self->{'limit'}) {
       $stmt .= qq{    
		   ORDER BY
		     sorting_id ASC  
	          };
     }
       dumpSQL($dbh_child,$stmt);

       create_and_load($self->{'dbVar'}, "$subind\_tmp_gty", 'subsnp_id i*', 'ind_id i', 'length_pat i','allele_1', 'allele_2','submitted_strand i');
  print $logh Progress::location();
       debug(localtime() . "\tTime finishing to dump $subind\_tmp_gty table: ",scalar(localtime(time)));

     POSIX:_exit(0);
     }
     $rec_pid1{$pid}++;
     if (keys %rec_pid1 == $num_process){
       do{
	 #the parent waits for the child to finish; foreach my $p (keys
	 #%rec_pid) {
	  print $logh Progress::location();
	  debug(localtime() . "\tWaiting for a free slot to fork another process");
	 $kid = waitpid(-1,WNOHANG); if ($kid > 0){
	   delete $rec_pid1{$kid};
	 }
       } until keys %rec_pid1 < $num_process;
     }
   }
   do{
     #the parent waits for the child to finish for the last three processes;
     #foreach my $p (keys %rec_pid) {
      print $logh Progress::location();
      debug(localtime() . "\tWaiting for the last " . scalar(keys(%rec_pid1)) . " forked processes to finish");
     my $kid = waitpid(-1,WNOHANG);  
     if ($kid > 0){
       delete $rec_pid1{$kid};
     }
   } until keys %rec_pid1 <1;

    print $logh Progress::location();
    debug(localtime() . "\tAll forked processes have finished, proceeding");
    
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
                   $self->{'dbSNP_share_db'}..Allele a1, 
                   $self->{'dbSNP_share_db'}..Allele a2
                 WHERE 
                   a1.rev_allele_id = a2.allele_id
                };
      dumpSQL($self->{'dbSNP'},$stmt);

      create_and_load($self->{'dbVar'}, "tmp_rev_allele", "allele_id i*","allele *", "rev_allele");
      print $logh Progress::location();
    }

#     #we have truncated the individual_genotype table, one contains the genotypes single bp, and the other, the rest
#     #necessary to create a unique index to remove duplicated genotypes


  #ÊGet the create statement for tmp_individual_genotype_single_bp from master schema
  my $ind_gty_stmt = get_create_statement('tmp_individual_genotype_single_bp',$self->{'schema_file'});
  print $logh Progress::location();

  #ÊDrop the tmp_individual.. table if it exists
  $stmt = qq{
    DROP TABLE IF EXISTS
      tmp_individual_genotype_single_bp
  };
  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();
  
   my $insert = 'INSERT';
   if ($self->{'dbm'}->dbCore()->species !~ /homo|hum/i) {
    
    #ÊCreate the tmp_individual_genotype_single_bp table
    $self->{'dbVar'}->do($ind_gty_stmt);
    print $logh Progress::location();
    
     ###only non human one needs unique index
     $insert = "INSERT IGNORE";

     $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX ind_genotype_idx ON tmp_individual_genotype_single_bp(variation_id,subsnp_id,sample_id,allele_1,allele_2)});

     old_calculate_gtype($self,$self->{'dbVar'},"SubInd_tmp_gty","tmp_individual_genotype_single_bp",$insert);
     $self->{'dbVar'}->do(qq{RENAME TABLE SubInd_tmp_gty TO tmp2_gty});

    print $logh Progress::location();
#     $self->{'dbVar'}->do(qq{DROP INDEX ind_genotype_idx ON tmp_individual_genotype_single_bp});

   }
   else {
     my %rec_pid;
     my $kid;

     foreach my $subind (@subind_tables) {
       my $table1 = "$subind\_tmp_gty"; #CHANGED
       my $table2 = "tmp_individual_genotype_single_bp\_$subind";
       push @tmp2_gty_tables, $table1;
       push @gty_tables, $table2;

       #next unless $subind eq "SubInd_ch1";
       debug(localtime() . "\ttable1 is $table1 and table2 is $table2");
       
      #ÊCreate a sub table for the genotypes
      $stmt = $ind_gty_stmt;
      $stmt =~ s/tmp_individual_genotype_single_bp/$table2/;
      $self->{'dbVar'}->do($stmt);
      print $logh Progress::location();
       
      $self->{'dbVar'}->do(qq{ALTER TABLE $table2 engine = innodb});
      print $logh Progress::location();

       #need to fork to processce gtype data
       my $pid = fork;
       if (! defined $pid){
	 throw("Not possible to fork: $!\n");
       }
       elsif ($pid == 0){
	 #you are the child.....
	 old_calculate_gtype($self,$self->{'dbVar'},$table1,$table2,$insert);
       POSIX:_exit(0);
       }
       $rec_pid{$pid}++;
       if (keys %rec_pid == $num_process){
	 do{
	 #the parent waits for the child to finish;
	 #foreach my $p (keys %rec_pid) {
	  print $logh Progress::location();
	  debug(localtime() . "\tWaiting for a free slot to fork another process");
	   $kid = waitpid(-1,WNOHANG);  
	   if ($kid > 0){
	     delete $rec_pid{$kid};
	   }
	 } until keys %rec_pid < $num_process;
       }
     }

     do{
       #the parent waits for the child to finish for the last three processes;
	print $logh Progress::location();
	debug(localtime() . "\tWaiting for the last " . scalar(keys(%rec_pid)) . " forked processes to finish");
       $kid = waitpid(-1,WNOHANG);  
       if ($kid > 0){
	 delete $rec_pid{$kid};
       }
     } until keys %rec_pid < 1;

    print $logh Progress::location();
    debug(localtime() . "\tAll forked processes have finished, proceeding with merging the tmp_individual_single_bp table");     
     
     my $merge_gty_tables = join ',',@gty_tables;
     
     #ÊI (Pontus) think this next statement messes up the merge statement because it refers to a non-existing table (MySQL crashes when trying to access the merged table).
     # It is fine if creating the merged table manually from the existing tmp_individual...SubInd.. tables. At the moment (2010-06-01), I cannot see any use for this next statement so I comment it out.
     # $merge_gty_tables .= ",tmp_individual_genotype_single_bp_last_insert"; #to hold new insert from ensembl called SNPs
     
     my $merge_tmp2_gty_tables = join ',',@tmp2_gty_tables;
     #make a merge table for tmp_individual_genotype_single_bp
    
    # Merge all the sub tables into a big tmp_individual_single_bp table
    $stmt = $ind_gty_stmt;
    $stmt =~ s/MyISAM/MERGE INSERT_METHOD=LAST UNION=($merge_gty_tables)/;
    $self->{'dbVar'}->do($stmt);
    print $logh Progress::location();
    
    $stmt = qq{
      CREATE TABLE
	tmp2_gty (
	  subsnp_id int(11) DEFAULT NULL,
	  ind_id int(11) DEFAULT NULL,
	  length_pat int(11) DEFAULT NULL,
	  allele_1 varchar(255) DEFAULT NULL,
	  allele_2 varchar(255) DEFAULT NULL,
	  submitted_strand int(11) DEFAULT NULL,
	  KEY subsnp_id_idx (subsnp_id)
	)
	ENGINE = MERGE,
	INSERT_METHOD = LAST,
	UNION = ($merge_tmp2_gty_tables)
    };
    $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();
  }
  
  debug(localtime() . "\tTime finishing to insert to tmp_individual_genotype_single_bp table: ",scalar(localtime(time)));


  debug(localtime() . "\tTime start to insert to individual_genotype_multiple_bp table: ",scalar(localtime(time)));


  $self->{'dbVar'}->do(qq{INSERT INTO individual_genotype_multiple_bp (variation_id, subsnp_id, sample_id, allele_1, allele_2) 
			     SELECT vs.variation_id, vs.subsnp_id, s.sample_id,
                             IF(vs.substrand_reversed_flag,tra1.allele,tra1.rev_allele) as allele_1,
                             IF(vs.substrand_reversed_flag,tra2.allele,tra2.rev_allele) as allele_2
                             FROM tmp2_gty tg, variation_synonym vs, sample s, tmp_rev_allele tra1, tmp_rev_allele tra2
			     WHERE tg.subsnp_id = vs.subsnp_id
                             AND   tg.submitted_strand IN (1,3,5)
                             AND   tra1.allele = tg.allele_1
                             AND   tra2.allele = tg.allele_2
			     AND   tg.length_pat > 3
			     AND   tg.ind_id = s.individual_id});
  print $logh Progress::location();

  $self->{'dbVar'}->do(qq{INSERT INTO individual_genotype_multiple_bp (variation_id,subsnp_id, sample_id, allele_1, allele_2) 
			     SELECT vs.variation_id, vs.subsnp_id, s.sample_id,
                             IF(vs.substrand_reversed_flag,tra1.rev_allele,tra1.allele) as allele_1,
                             IF(vs.substrand_reversed_flag,tra2.rev_allele,tra2.allele) as allele_2
	                     FROM tmp2_gty tg, variation_synonym vs, sample s, tmp_rev_allele tra1, tmp_rev_allele tra2
			     WHERE tg.subsnp_id = vs.subsnp_id
                             AND   tg.submitted_strand IN (0,2,4)
                             AND   tra1.allele = tg.allele_1
                             AND   tra2.allele = tg.allele_2
			     AND   tg.length_pat > 3
			     AND   tg.ind_id = s.individual_id});
  print $logh Progress::location();


  #$self->{'dbVar'}->do("DROP TABLE tmp2_gty");
  debug(localtime() . "\tTime finishing to insert to individual_genotype_multiple_bp table: ",scalar(localtime(time)));

  debug(localtime() . "\tTime to delete from tmp_individual_genotype_sangle_bp where length(allele_1)>1 and insert the data into individual_genotype_multiple_bp");

  $self->{'dbVar'}->do(qq{INSERT INTO individual_genotype_multiple_bp select * 
                       FROM tmp_individual_genotype_single_bp
                       WHERE length(allele_1)>1 or length(allele_2)>1
                       OR allele_1 = '-' or allele_2 = '-'});
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{DELETE FROM tmp_individual_genotype_single_bp
                          WHERE length(allele_1)>1 or length(allele_2)>1
                          OR allele_1 = '-' or allele_2 = '-'});
  print $logh Progress::location();
   $self->{'dbVar'}->do(qq{INSERT INTO tmp_individual_genotype_single_bp select * 
                       FROM individual_genotype_multiple_bp
                       WHERE length(allele_1)=1 and length(allele_2=1)
                       AND allele_1 != '-' and allele_2 != '-'});
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{DELETE FROM individual_genotype_multiple_bp
                          WHERE length(allele_1)=1 and length(allele_2)=1
                          AND allele_1 != '-' and allele_2 != '-'});
  print $logh Progress::location();
  $self->{'dbVar'}->do(qq{DELETE FROM individual_genotype_multiple_bp
                          WHERE allele_1 like "%indeterminate%"});
  print $logh Progress::location();
   
}

sub calculate_gtype {
  my $self = shift;
  my $dbVar = shift;
  my $dbSNP = shift;
  my $subind_table = shift;
  my $dst_table = shift;
  my $mlt_file = shift;
  my $allele_file = shift;
  my $sample_file = shift;
  
  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  print $logh Progress::location() . "\tImporting from $subind_table to $dst_table\n";
  
  #ÊProcess the alleles in chunks based on the SubSNP id
  my $chunksize = 1e6;
  
  my $stmt = qq{
    SELECT
      MIN(subsnp_id),
      MAX(subsnp_id)
    FROM
      $subind_table
  };
  my ($min_ss,$max_ss) = @{$dbSNP->selectall_arrayref($stmt)->[0]};
  print $logh Progress::location();
  
  # If there were no subsnp_ids to be found in the table, return
  if (!defined($min_ss) || !defined($max_ss)) {
    print $logh Progress::location() . "\tNo SubSNP ids available in $subind_table, skipping this table\n";
    return 1;
  }
  
  # The tempfile to be used for loading
  my $tempfile = $self->{'tmpdir'} . '/' . $dst_table . '_' . $self->{'tmpfile'};
  
  #ÊA prepared statement for getting the genotype data from the dbSNP mirror
  $stmt = qq{
    SELECT 
      si.subsnp_id,
      sind.ind_id, 
      ug.allele_id_1,
      ug.allele_id_2,
      CASE WHEN
	si.submitted_strand_code IS NOT NULL
      THEN
	si.submitted_strand_code
      ELSE
	0
      END,
      LEN(ug.gty_str) AS pattern_length
    FROM   
      $subind_table si JOIN 
      $self->{'dbSNP_share_db'}..GtyAllele ga ON (
	ga.gty_id = si.gty_id
      ) JOIN
      $self->{'dbSNP_share_db'}..UniGty ug ON (
	ug.unigty_id = ga.unigty_id
      ) JOIN
      SubmittedIndividual sind ON (
	sind.submitted_ind_id = si.submitted_ind_id
      )
    WHERE
      si.subsnp_id BETWEEN ? AND ? AND
      si.submitted_strand_code <= 5
  };
  my $subind_sth = $dbSNP->prepare($stmt);
  
  #ÊPrepared statement to get the corresponding variation_ids for a range of subsnps from variation_synonym
  $stmt = qq{
    SELECT
      vs.subsnp_id,
      vs.variation_id,
      vs.substrand_reversed_flag
    FROM
      variation_synonym vs
    WHERE
      vs.subsnp_id BETWEEN ? AND ?
  };
  my $vs_sth = $dbVar->prepare($stmt);
  
  # Prepared statement to get the sample_id from a ind_id
  $stmt = qq{
    SELECT
      s.sample_id
    FROM
      sample s
    WHERE
      s.individual_id = ?
    LIMIT 1
  };
  my $sample_sth = $dbVar->prepare($stmt);
  
  # Prepared statement to get the alleles. We get each one of these when needed.
  $stmt = qq{
    SELECT
      a.allele,
      ar.allele
    FROM
      $self->{'dbSNP_share_db'}..Allele a JOIN
      $self->{'dbSNP_share_db'}..Allele ar ON (
	ar.allele_id = a.rev_allele_id
      )
    WHERE
      a.allele_id = ?
  };
  my $allele_sth = $dbSNP->prepare($stmt);
  
  #ÊHash to hold the alleles in memory
  my %alleles;
  #ÊIf the allele file exist, we'll read alleles from it
  if (defined($allele_file) && -e $allele_file) {
    open(ALL,'<',$allele_file);
    # Lock the filehandle
    flock(ALL,LOCK_SH);
    # Parse the file
    while (<ALL>) {
      my ($id,$a,$arev) = split;
      $alleles{$id} = [$a,$arev];
    }
    close(ALL);
  }
  # Keep track if we did any new lookups
  my $new_alleles = 0;
  
  #ÊHash to keep sample_ids in memory
  my %samples;
  #ÊIf the sample file exist, we'll read alleles from it
  if (defined($sample_file) && -e $sample_file) {
    open(SPL,'<',$sample_file);
    # Lock the filehandle
    flock(SPL,LOCK_SH);
    # Parse the file
    while (<SPL>) {
      my ($ind_id,$sample_id) = split;
      $samples{$ind_id} = $sample_id;
    }
    close(SPL);
  }
  # The individual_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value
  $samples{0} = '\N';
  # Keep track if we did any new lookups
  my $new_samples = 0;
  
  # Open a file handle to the temp file that will be used for loading
  open(SGL,'>',$tempfile) or die ("Could not open import file $tempfile for writing");
  #ÊGet a lock on the file
  flock(SGL,LOCK_EX);
  
  # Open a file handle to the file that will be used for loading the multi-bp genotypes
  open(MLT,'>>',$mlt_file) or die ("Could not open import file $mlt_file for writing");
  #ÊGet a lock on the file
  flock(MLT,LOCK_EX);
  
  my $ss_offset = ($min_ss - 1);
  while ($ss_offset < $max_ss) {
    
    print $logh Progress::location() . "\tProcessing genotype data for SubSNP ids " . ($ss_offset + 1) . "-" . ($ss_offset + $chunksize) . "\n";
    
    # For debugging
    $ss_offset = ($max_ss - 10000);
    
    # First, get all SubSNPs
    $subind_sth->bind_param(1,($ss_offset + 1),SQL_INTEGER);
    $vs_sth->bind_param(1,($ss_offset + 1),SQL_INTEGER);
    $ss_offset += $chunksize;
    $subind_sth->bind_param(2,$ss_offset,SQL_INTEGER);
    $vs_sth->bind_param(2,$ss_offset,SQL_INTEGER);
    
    # Fetch the import data as an arrayref
    $subind_sth->execute();
    my $import_data = $subind_sth->fetchall_arrayref();
    print $logh Progress::location();
    
    # Fetch the subsnp_id to variation_id mapping and store as hashref
    $vs_sth->execute();
    my $variation_ids = $vs_sth->fetchall_hashref(['subsnp_id']);
    print $logh Progress::location();
    
    #ÊNow, loop over the import data and print it to the tempfile so we can import the data. Replace the allele_id with the corresponding allele on-the-fly
    while (my $data = shift(@{$import_data})) {
      my $subsnp_id = $data->[0];
      my $ind_id = $data->[1];
      my $allele_1 = $data->[2];
      my $allele_2 = $data->[3];
      my $sub_strand = $data->[4];
      my $pattern_length = $data->[5];
      
      # If pattern length is less than 3, skip this genotype
      next if ($pattern_length < 3);
      
      #ÊShould the alleles be flipped?
      my $reverse = (($sub_strand + $variation_ids->{$subsnp_id}{'substrand_reversed_flag'})%2 != 0);
      
      # Look up the alleles if necessary
      foreach my $allele_id ($allele_1,$allele_2) {
	if (!exists($alleles{$allele_id})) {
	  $allele_sth->execute($allele_id);
	  my ($a,$arev);
	  $allele_sth->bind_columns(\$a,\$arev);
	  $allele_sth->fetch();
	  $alleles{$allele_id} = [$a,$arev];
	  $new_alleles = 1;
	}
      }
      #ÊLook up the sample id in the database if necessary
      if (!exists($samples{$ind_id})) {
        $sample_sth->execute($ind_id);
        my $sample_id;
        $sample_sth->bind_columns(\$sample_id);
        $sample_sth->fetch();
        $samples{$ind_id} = $sample_id;
        $new_samples = 1;
      }
      
      # Print to the import file
      $allele_1 = $alleles{$allele_1}->[$reverse];
      $allele_2 = $alleles{$allele_2}->[$reverse];
      
      #ÊSkip this genotype if the alleles are N
      next if ($allele_1 eq 'N' && $allele_2 eq 'N');
      
      my $row = join("\t",($variation_ids->{$subsnp_id}{'variation_id'},$subsnp_id,$samples{$ind_id},$allele_1,$allele_2)) . "\n";
      if (length($allele_1) == 1 && length($allele_2) == 1 && $allele_1 ne '-' && $allele_2 ne '-') {
	print SGL $row;
      }
      else {
	#ÊSkip this genotype if the first allele contains the string 'indeterminate'
	next if ($allele_1 =~ m/indeterminate/i);
	print MLT $row;
      }
    }
    print $logh Progress::location();
    
    # Just for debugging!
    last;
  }
  
  close(SGL);
  close(MLT);
  
  #ÊDisable the keys on the destination table before loading in order to speed up the insert
  $stmt = qq {
    ALTER TABLE
      $dst_table
    DISABLE KEYS
  };
  $dbVar->do($stmt);
  print $logh Progress::location();
  
  # Load the data from the infile
  print $logh Progress::location() . "\tLoads the infile $tempfile to $dst_table\n";
  loadfile($tempfile,$dbVar,$dst_table,"variation_id","subsnp_id","sample_id","allele_1","allele_2");
  print $logh Progress::location();
  
  # Enable the keys after the inserts
  $stmt = qq {
    ALTER TABLE
      $dst_table
    ENABLE KEYS
  };
  $dbVar->do($stmt);
  print $logh Progress::location();
  
  # Remove the import file
  unlink($tempfile);
  
  # If we had an allele file and we need to update it, do that
  if (defined($allele_file) && -e $allele_file && $new_alleles) {
    open(ALL,'>',$allele_file);
    # Lock the file
    flock(ALL,LOCK_EX);
    while (my ($id,$a) = each(%alleles)) {
      print ALL join("\t",($id,@{$a})) . "\n";
    }
    close(ALL);
  }
  
  # If we had a sample file and we need to update it, do that
  if (defined($sample_file) && -e $sample_file && $new_samples) {
    open(SPL,'>',$sample_file);
    # Lock the file
    flock(SPL,LOCK_EX);
    while (my @row = each(%samples)) {
      print SPL join("\t",@row) . "\n";
    }
    close(SPL);
  }
  
  print $logh Progress::location() . "\tFinished importing from $subind_table to $dst_table\n";
  return 1;
}

sub old_calculate_gtype {
  
  my ($self,$dbVariation,$table1,$table2,$insert) = @_;

  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};

  #ÊThe different cases for reversed alleles were split into 4 different queries. As a first step towards optimizing,
  #ÊI combine these into one query with a conditional based on the value of substrand_reversed_flag + submitted_strand
  # An even value means that we'll get the forward allele. An odd value means the reverse allele
  
  # I haven't tested this code so for now, I'll leave it commented out and use the old one
=head  
  my $stmt = qq{
    $insert INTO
      $table2 (
	variation_id,
	subsnp_id,
	sample_id,
	allele_1,
	allele_2
      )
      SELECT
	vs.variation_id,
	vs.subsnp_id,
	s.sample_id,
	IF (
	  MOD(tg.submitted_strand+vs.substrand_reversed_flag,2) = 0,
	  tra1.allele,
	  tra1.rev_allele
	),
	IF (
	  MOD(tg.submitted_strand+vs.substrand_reversed_flag,2) = 0,
	  tra2.allele,
	  tra2.rev_allele
	)
      FROM
	$table1 tg,
	variation_synonym vs,
	tmp_rev_allele tra1,
	tmp_rev_allele tra2,
	sample s
      WHERE
	tg.subsnp_id = vs.subsnp_id AND
	tra1.allele = tg.allele_1 AND
	tra2.allele = tg.allele_2 AND
	tg.ind_id = s.individual_id AND
	tg.length_pat = 3 AND
	tg.submitted_strand <= 5
  };
  my $sth = $dbVariation->prepare($stmt);

  debug(localtime() . "\tTime starting to insert alleles from $table1 into $table2 : ",scalar(localtime(time)));
  $sth->execute();
  print $logh Progress::location();
=cut

  debug(localtime() . "\tTime starting to insert into $table2 : ",scalar(localtime(time)));
  $dbVariation->do(qq{$insert INTO $table2 (variation_id, subsnp_id, sample_id, allele_1, allele_2) 
				 SELECT vs.variation_id, vs.subsnp_id, s.sample_id, tra1.allele, tra2.allele
				 FROM $table1 tg, variation_synonym vs, tmp_rev_allele tra1, tmp_rev_allele tra2, sample s
				 WHERE tg.submitted_strand IN (1,3,5)
				 AND vs.substrand_reversed_flag =1
				 AND tg.subsnp_id = vs.subsnp_id
				 AND tra1.allele = tg.allele_1
				 AND tra2.allele = tg.allele_2
				 AND tg.length_pat = 3
				 AND tg.ind_id = s.individual_id});
  print $logh Progress::location();
  debug(localtime() . "\tTime starting to insert into $table2 table 2: ",scalar(localtime(time)));

  $dbVariation->do(qq{$insert INTO $table2 (variation_id, subsnp_id, sample_id, allele_1, allele_2)
				 SELECT vs.variation_id, vs.subsnp_id, s.sample_id, tra1.rev_allele, tra2.rev_allele
				 FROM $table1 tg, variation_synonym vs, tmp_rev_allele tra1, tmp_rev_allele tra2, sample s
				 WHERE tg.submitted_strand IN (1,3,5)
				 AND vs.substrand_reversed_flag =0
				 AND tg.subsnp_id = vs.subsnp_id
				 AND tra1.allele = tg.allele_1
				 AND tra2.allele = tg.allele_2
				 AND tg.length_pat = 3
 				 AND tg.ind_id = s.individual_id});
  print $logh Progress::location();

  debug(localtime() . "\tTime starting to insert into $table2 table 3: ",scalar(localtime(time)));

  $dbVariation->do(qq{$insert INTO $table2 (variation_id, subsnp_id, sample_id, allele_1, allele_2)

				 SELECT vs.variation_id, vs.subsnp_id, s.sample_id, tra1.rev_allele, tra2.rev_allele
				 FROM $table1 tg, variation_synonym vs, tmp_rev_allele tra1, tmp_rev_allele tra2, sample s
				 WHERE tg.submitted_strand IN (0,2,4)
				 AND vs.substrand_reversed_flag =1
				 AND tg.subsnp_id = vs.subsnp_id
				 AND tra1.allele = tg.allele_1
				 AND tra2.allele = tg.allele_2
				 AND tg.length_pat = 3
				 AND tg.ind_id = s.individual_id});
  print $logh Progress::location();

  debug(localtime() . "\tTime starting to insert into $table2 table 4: ",scalar(localtime(time)));

  $dbVariation->do(qq{$insert INTO $table2 (variation_id, subsnp_id, sample_id, allele_1, allele_2)

				 SELECT vs.variation_id, vs.subsnp_id, s.sample_id, tra1.allele, tra2.allele
				 FROM $table1 tg, variation_synonym vs, tmp_rev_allele tra1, tmp_rev_allele tra2, sample s
				 WHERE tg.submitted_strand IN (0,2,4)
				 AND vs.substrand_reversed_flag =0
				 AND tg.subsnp_id = vs.subsnp_id
				 AND tra1.allele = tg.allele_1
				 AND tra2.allele = tg.allele_2
				 AND tg.length_pat = 3
				 AND tg.ind_id = s.individual_id});
  print $logh Progress::location();

  $dbVariation->do(qq{alter TABLE $table2 engine = MyISAM}) if ($self->{'dbm'}->dbCore()->species =~ /homo|hum/i);
  print $logh Progress::location();
  
}

#
# loads population genotypes into the population_genotype table
#
sub population_genotypes {
    my $self = shift;

  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
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
		   $self->{'dbSNP_share_db'}..Allele a1, 
		   $self->{'dbSNP_share_db'}..Allele a2
                 WHERE 
                   a1.rev_allele_id = a2.allele_id
                };
      dumpSQL($self->{'dbSNP'},$stmt);

      create_and_load($self->{'dbVar'}, "tmp_rev_allele", "allele_id i*","allele *", "rev_allele");
      print $logh Progress::location();
      
    }
    debug(localtime() . "\tDumping GtyFreqBySsPop and UniGty data");
 
    $stmt = "SELECT ";
    if ($self->{'limit'}) {
      $stmt .= "TOP $self->{'limit'} ";
    }
    $stmt .= qq{
                  gtfsp.subsnp_id AS sorting_id, 
                  gtfsp.pop_id, 
                  gtfsp.freq,
                  ROUND(gtfsp.cnt,0),
                  a1.allele, 
                  a2.allele
	       	FROM   
	       	  GtyFreqBySsPop gtfsp, 
	       	  $self->{'dbSNP_share_db'}..UniGty ug, 
	       	  $self->{'dbSNP_share_db'}..Allele a1, 
	       	  $self->{'dbSNP_share_db'}..Allele a2
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
    dumpSQL($self->{'dbSNP'},$stmt);

     debug(localtime() . "\tloading population_genotype data");

     create_and_load($self->{'dbVar'}, "tmp_pop_gty", 'subsnp_id i*', 'pop_id i*', 'freq','count','allele_1', 'allele_2');
  print $logh Progress::location();

   $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX pop_genotype_idx ON population_genotype(variation_id,subsnp_id,frequency,sample_id,allele_1(5),allele_2(5))});


    $self->{'dbVar'}->do(qq{INSERT IGNORE INTO population_genotype (variation_id,subsnp_id,allele_1, allele_2, frequency, sample_id, count)
			                  SELECT vs.variation_id,vs.subsnp_id,tra1.rev_allele as allele_1,tra2.rev_allele as allele_2,tg.freq,s.sample_id, tg.count
					  FROM   variation_synonym vs, tmp_pop_gty tg,tmp_rev_allele tra1,tmp_rev_allele tra2, sample s
					  WHERE  vs.subsnp_id = tg.subsnp_id
                                          AND    tg.allele_1 = tra1.allele
                                          AND    tg.allele_2 = tra2.allele
                                          AND    vs.substrand_reversed_flag = 1
					  AND    s.pop_id = tg.pop_id});
  print $logh Progress::location();

    $self->{'dbVar'}->do(qq{INSERT IGNORE INTO population_genotype (variation_id,subsnp_id,allele_1, allele_2, frequency, sample_id, count)
					  SELECT vs.variation_id,vs.subsnp_id,tg.allele_1,tg.allele_2,tg.freq,s.sample_id, tg.count
					  FROM   variation_synonym vs, tmp_pop_gty tg,sample s
					  WHERE  vs.subsnp_id = tg.subsnp_id
                                          AND    vs.substrand_reversed_flag = 0
					  AND    s.pop_id = tg.pop_id});
  print $logh Progress::location();

#    $self->{'dbVar'}->do(qq{DROP INDEX pop_genotype_idx ON population_genotype});
    #$self->{'dbVar'}->do("DROP TABLE tmp_pop_gty");
}

# cleans up some of the necessary temporary data structures after the
# import is complete
sub cleanup {
    my $self = shift;

  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
    debug(localtime() . "\tIn cleanup...");
    #remove populations that are not present in the Individual or Allele table for the specie
    $self->{'dbVar'}->do('CREATE TABLE tmp_pop (sample_id int PRIMARY KEY)'); #create a temporary table with unique populations
  print $logh Progress::location();
    $self->{'dbVar'}->do('INSERT IGNORE INTO tmp_pop SELECT distinct(sample_id) FROM allele'); #add the populations from the alleles
  print $logh Progress::location();
    $self->{'dbVar'}->do('INSERT IGNORE INTO tmp_pop SELECT distinct(sample_id) FROM population_genotype'); #add the populations from the population_genotype
  print $logh Progress::location();
    $self->{'dbVar'}->do('INSERT IGNORE INTO tmp_pop SELECT population_sample_id FROM individual_population'); #add the populations from the individuals
  print $logh Progress::location();
    $self->{'dbVar'}->do(qq{INSERT IGNORE INTO tmp_pop SELECT super_population_sample_id 
 				      FROM population_structure ps, tmp_pop tp 
 				      WHERE tp.sample_id = ps.sub_population_sample_id}); #add the populations from the super-populations
  print $logh Progress::location();
    
    #necessary to difference between MySQL 4.0 and MySQL 4.1
    my $sql;
    my $sql_2;
    my $sql_3;
    my $sql_4;
    my $sth = $self->{'dbVar'}->prepare(qq{SHOW VARIABLES LIKE 'version'});
    $sth->execute();
    my $row_ref = $sth->fetch();
    $sth->finish();
    #check if the value in the version contains the 4.1
    $sql = qq{DELETE FROM p USING population p
		  LEFT JOIN tmp_pop tp ON p.sample_id = tp.sample_id
		  LEFT JOIN sample s on p.sample_id = s.sample_id
		  WHERE tp.sample_id is null AND s.individual_id is null};

    $sql_2 = qq{DELETE FROM ss USING sample_synonym ss
		    LEFT JOIN tmp_pop tp ON ss.sample_id = tp.sample_id
		    LEFT JOIN sample s on s.sample_id = ss.sample_id
		    WHERE tp.sample_id is null AND s.individual_id is null};

    $sql_3 = qq{DELETE FROM ps USING population_structure ps 
		    LEFT JOIN tmp_pop tp ON ps.sub_population_sample_id = tp.sample_id 
		    WHERE tp.sample_id is null};

    $sql_4 = qq{DELETE FROM s USING sample s
		    LEFT JOIN population p ON s.sample_id = p.sample_id
		    WHERE p.sample_id is null
		    AND s.individual_id is null
		    };

    $self->{'dbVar'}->do($sql); #delete from population
  print $logh Progress::location();
    # populations not present
    $self->{'dbVar'}->do($sql_2); #delete from sample_synonym
  print $logh Progress::location();
    $self->{'dbVar'}->do($sql_3); #and delete from the population_structure table
  print $logh Progress::location();
    $self->{'dbVar'}->do($sql_4); #and delete from Sample table the ones that are not in population
  print $logh Progress::location();

    $self->{'dbVar'}->do('DROP TABLE tmp_pop'); #and finally remove the temporary table
  print $logh Progress::location();

    $self->{'dbVar'}->do('ALTER TABLE variation  DROP COLUMN snp_id');
  print $logh Progress::location();  
    $self->{'dbVar'}->do('ALTER TABLE variation_synonym DROP COLUMN substrand_reversed_flag');
  print $logh Progress::location();
    $self->{'dbVar'}->do('ALTER TABLE sample DROP COLUMN pop_class_id, DROP COLUMN pop_id, DROP COLUMN individual_id');
  print $logh Progress::location();
    

}

1;

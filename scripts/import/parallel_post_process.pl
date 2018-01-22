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

use POSIX;
use Getopt::Long;
use ImportUtils qw(dumpSQL load create_and_load debug);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
#use DBH;
use DBI qw(:sql_types);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use FindBin qw( $Bin );
use Data::Dumper;
use constant MAX_GENOTYPES => 6_000_000; #max number of genotypes per file. When more, split the file into regions
use constant REGIONS => 15; #number of chunks the file will be splited when more than MAX_GENOTYPES

my ($TMP_DIR, $TMP_FILE, $LIMIT, $PERLBIN, $SCHEDULER);


my ($vhost, $vport, $vdbname, $vuser, $vpass,
    $chost, $cport, $cdbname, $cuser, $cpass,
    $limit, $num_processes, $top_level, $species,
    $variation_feature, $flanking_sequence,
    $transcript_variation, $ld_populations, $reverse_things, $make_allele_string_table,$merge_rs_features,
    $merge_ensembl_snps,$new_source_id,$remove_wrong_variations,$remove_multi_tables,$check_seq_region_id, 
    $read_updated_files,$merge_multi_tables,
    $registry_file, $bsub_queue_name,@into_source_id);

$variation_feature = $flanking_sequence = $transcript_variation = $ld_populations = $reverse_things = $merge_rs_features = $merge_ensembl_snps = $new_source_id = $remove_wrong_variations = $remove_multi_tables = $make_allele_string_table = $check_seq_region_id = $read_updated_files = $merge_multi_tables = '';

GetOptions(
  'tmpdir=s'  => \$ImportUtils::TMP_DIR,
  'tmpfile=s' => \$ImportUtils::TMP_FILE,
  'species=s' => \$species,
  'limit=i'   => \$limit,
  'perlbin=s' => \$PERLBIN,
  'scheduler=s' => \$SCHEDULER,
  'top_level=i' => \$top_level,
  'num_processes=i' => \$num_processes,
  'variation_feature' => \$variation_feature,
  'flanking_sequence' => \$flanking_sequence,
  'transcript_variation' => \$transcript_variation,
  'reverse_things'       => \$reverse_things,
  'merge_rs_features'         => \$merge_rs_features,
  'remove_wrong_variations' => \$remove_wrong_variations,
  'remove_multi_tables'  => \$remove_multi_tables,
  'merge_ensembl_snps'   => \$merge_ensembl_snps,
  'new_source_id=n'      => \$new_source_id,
  'into_source_id=i'      => \@into_source_id,
  'check_seq_region_id'  => \$check_seq_region_id,
  'read_updated_files'   => \$read_updated_files,
  'make_allele_string_table' => \$make_allele_string_table,
  'merge_multi_tables' => \$merge_multi_tables,
  'ld_populations' => \$ld_populations,
  'registry_file=s' => \$registry_file,
  'bsub_queue_name=s' => \$bsub_queue_name,
);

$num_processes ||= 1;

$bsub_queue_name ||= 'long';

$LIMIT = ($limit) ? " $limit " : ''; #will refer to position in a slice
$PERLBIN   ||= '/software/bin/env perl'; # Could use $^X?
-x $PERLBIN or usage("Perl interpretor at $PERLBIN is not suitable. " .
                     "See the -perlbin argument" );
$SCHEDULER ||= 'LSF';
$SCHEDULER = uc($SCHEDULER);

my %valid_schedulers = map{$_=>1} qw( LSF PBS );
$valid_schedulers{$SCHEDULER} || usage("-scheduler $SCHEDULER is invalid");

usage('-num_processes must at least be 1') if ($num_processes == 0);
usage('-species argument required') if(!$species);

warn("Make sure you have an updated ensembl.registry file!\n");

$registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core')
    || usage( "Cannot find core db for $species in $registry_file" );
my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation')
    || usage( "Cannot find variation db for $species in $registry_file" );


my $dbVar = $vdba->dbc->db_handle;
my $dbCore = $cdba;

#added default options
$chost = $cdba->dbc->host; 
$cuser = $cdba->dbc->username;
$cport = $cdba->dbc->port;
$cdbname = $cdba->dbc->dbname;

$vhost = $vdba->dbc->host;
$vport = $vdba->dbc->port;
$vuser =  $vdba->dbc->username;
$vdbname = $vdba->dbc->dbname;
$vpass = $vdba->dbc->password;

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;


##Apart from mosquitos, we directly import top_level coordinates from dbSNP
if (! defined $top_level && $species !~ /mosquito|anoph/i) {
  $top_level=1;
}
#find out seq_region_id for haplotype chromosomes    
my $sr_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{select sra.seq_region_id from seq_region_attrib sra, attrib_type at where sra.attrib_type_id=at.attrib_type_id and at.name="Non Reference"});

my %hap_seq_id = map {$_->[0],1} @$sr_ref;
my $hap_id_string = "(". join ",", keys %hap_seq_id, .")";

$hap_id_string = '(0)' if $hap_id_string eq '()';

#we need to create a tmp table to store variations that we filter out
#create_failed_variation_table($dbVar);
parallel_variation_feature($dbVar, $top_level, $Bin . '/parallel_variation_feature.pl', $bsub_queue_name) if ($variation_feature);
parallel_flanking_sequence($dbVar, $Bin . '/parallel_flanking_sequence.pl', $bsub_queue_name) if ($flanking_sequence);
parallel_transcript_variation($dbVar, $Bin . '/parallel_transcript_variation.pl', $registry_file, $bsub_queue_name) if ($transcript_variation);
parallel_ld_populations($dbVar, $bsub_queue_name) if ($ld_populations);
reverse_things($dbVar) if ($reverse_things);
merge_rs_feature($dbVar) if ($merge_rs_features);
remove_wrong_variations($dbVar, $Bin . '/../../sql') if ($remove_wrong_variations);
remove_multi_tables($dbVar) if ($remove_multi_tables);
merge_ensembl_snps($dbVar,$new_source_id,\@into_source_id) if ($merge_ensembl_snps and $new_source_id);
make_allele_string_table($dbVar) if $make_allele_string_table;
read_updated_files($dbVar) if $read_updated_files;
merge_multi_tables($dbVar) if $merge_multi_tables;

#will take the number of processes, and divide the total number of entries in the variation_feature table by the number of processes
sub parallel_variation_feature{
    my $dbVar = shift;
    my $top_level = shift;
    my $script = shift;
    my $bsub_queue_name = shift;
    
    #before do anything, create a copy of old table
    $dbVar->do(qq{CREATE TABLE variation_feature_before_pp like variation_feature});
    $dbVar->do(qq{INSERT INTO variation_feature_before_pp select * from variation_feature});

    my $min_variation; #minim variation_feature_id
    my $max_variation; #maximum variation_feature_id
    my $call;
    my $variation_status_file = "status_file_variation_feature_$$\.log";
    #first, create the log file for the variation_feature
    open STATUS, ">$TMP_DIR/$variation_status_file"
	or throw("Could not open tmp file: $TMP_DIR/$variation_status_file\n"); 
    close STATUS;
    #then, calculate the rows for each subprocess
    my $sth = $dbVar->prepare(qq{SELECT min(variation_feature_id),max(variation_feature_id)
			   FROM variation_feature
    });
    $sth->execute();
    ($min_variation, $max_variation) = $sth->fetchrow_array();
    $sth->finish();
    my $dbname = $vdbname; #get the name of the database to create the file    

    if ($hap_id_string !~ /\([\,\d]+\)/) {
      $hap_id_string = "";
    }

    #create a temporary table to store the map_weight, that will be deleted by the last process
    $dbVar->do(qq{CREATE TABLE tmp_map_weight
                SELECT variation_id, count(*) as count
                FROM   variation_feature
                WHERE  seq_region_id not IN $hap_id_string
                GROUP BY variation_id}
	       );
    $dbVar->do(qq{ALTER TABLE tmp_map_weight 
		      ADD UNIQUE INDEX variation_idx(variation_id)});
    #add additional variation_ids only appear in haplotype chromosomes
    $dbVar->do(qq{INSERT IGNORE INTO tmp_map_weight
		  SELECT variation_id, count(*) as count
		  FROM   variation_feature
		  WHERE  seq_region_id IN $hap_id_string
		  GROUP BY variation_id});
    


    $dbVar->do(qq{DELETE FROM variation_feature WHERE seq_region_start=1});
    $dbVar->do(qq{DELETE FROM variation_feature WHERE seq_region_end=1});

    #if no allele_string table, create one from variation_feature table
    my $as_ref = $dbVar->selectall_arrayref(qq{show tables like 'allele_string'});
    my $as = $as_ref->[0][0];
    if (! $as_ref->[0][0]) {
      debug("There is no allele_string table, create one from variation_feature table ?");
      $dbVar->do(qq{CREATE TABLE allele_string (
        variation_id int(10) unsigned NOT NULL,
        allele longtext,
        UNIQUE KEY variation_allele_idx(variation_id,allele(4)))
        });
        $dbVar->do(qq{INSERT IGNORE INTO allele_string SELECT variation_id,substring_index(allele_string,'/',1) as allele
                      FROM variation_feature
                     });
        $dbVar->do(qq{INSERT IGNORE INTO allele_string SELECT variation_id,substring_index(allele_string,'/',-1) as allele
                      FROM variation_feature
                     });
	$dbVar->do(qq{INSERT IGNORE INTO allele_string SELECT variation_id,substring_index(substring_index(allele_string,'/',2),'/',-1) as allele
                      FROM variation_feature
		      });
        $dbVar->do(qq{INSERT IGNORE INTO allele_string SELECT variation_id,substring_index(substring_index(allele_string,'/',3),'/',-1) as allele
                      FROM variation_feature
                     });
    }

    my $sub_variation = int(($max_variation - $min_variation)/ $num_processes);

    $bsub_queue_name ||= 'normal';

    #the limit will be (AND variation_feature_id > min and variation_feature_id < max)
    for (my $i = 0; $i < $num_processes ; $i++){
	$limit = "AND variation_feature_id <= " . (($i+1) * $sub_variation + $min_variation-1) . " AND variation_feature_id >= " . ($i*$sub_variation + $min_variation) if ($i+1 < $num_processes);
	$limit =  "AND variation_feature_id <= " .  $max_variation . " AND variation_feature_id >= " . ($i*$sub_variation + $min_variation) if ($i + 1 == $num_processes); #the last one takes the left rows

	$call = "bsub -q $bsub_queue_name -J $dbname\_variation_job_$i -o $TMP_DIR/output_variation_feature_$i\_$$.txt $PERLBIN $script -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -limit '$limit' -tmpdir $TMP_DIR -tmpfile $TMP_FILE -num_processes $num_processes -status_file $variation_status_file ";
	$call .= "-cpass $cpass " if ($cpass);
	$call .= "-cport $cport " if ($cport);
	$call .= "-vpass $vpass " if ($vpass);
	$call .= "-toplevel $top_level " if ($top_level);
	#print $call,"\n";
	system($call);      
    }
    $call = "bsub -q $bsub_queue_name -K -w 'done($dbname\_variation_job*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
    system($call);

}

#Will take the number of processes, and divide the total number of entries in the flanking_sequence table by the number of processes
#has to wait until the variation_feature table has been filled
sub parallel_flanking_sequence{
  my $dbVar = shift;
  my $script = shift;
  my $bsub_queue_name = shift;

  my $pid = $$;
  my $flanking_status_file = "flanking_status_file_$pid\.log";

  # To avoid any confusion regarding the haplotypes and reversing etc, only attempt to map to non-haplotype seq regions 
  my $hap_line;
  if ($hap_id_string =~ /\([\,\d]+\)/) {
    $hap_line = "vf.seq_region_id NOT IN $hap_id_string";
  }
  else {
    $hap_line ='';
  }
  
  #before do anything, create a copy of old table, this is done in post_flanking_sequence.pl
  $dbVar->do(qq{CREATE TABLE flanking_sequence_before_pp like flanking_sequence});
  $dbVar->do(qq{INSERT INTO flanking_sequence_before_pp select * from flanking_sequence});

  #first, create the log file for the variation_feature
  open STATUS, ">$TMP_DIR/$flanking_status_file"
    or throw("Could not open tmp file: $TMP_DIR/$flanking_status_file\n"); 
  close STATUS;

  my $dbname = $vdbname; #get the name of the database to create the file

  #find out the total number of variations to split the into the files
  my $sequences; #total number of flanking sequences in table
  my $sth_variations = $dbVar->prepare(qq{SELECT COUNT(*) from flanking_sequence});
  $sth_variations->execute();
  ($sequences) = $sth_variations->fetch();
  $sth_variations->finish();

  my $sth = $dbVar->prepare(qq{
	SELECT fs.variation_id, fs.up_seq, fs.down_seq,
	vf.seq_region_id, vf.seq_region_start,
	vf.seq_region_end, vf.seq_region_strand
	FROM flanking_sequence fs FORCE INDEX (PRIMARY) LEFT JOIN variation_feature vf
	ON vf.variation_id = fs.variation_id
	WHERE $hap_line
	ORDER BY fs.variation_id
	$LIMIT
  }, {mysql_use_result => 1});

  $sth->execute();
  my $count = 0; #to know the number of variations
  my $process = 1; #number of file to write the process to
  my $sub_sequences = int($sequences->[0] / $num_processes);
  my $previous_variation_id = 0;
  my $curr_variation_id;
  my $buffer = {}; #hash containing all the files to be parallelized
  
  # Open a new file
  my $file = $TMP_DIR . "/" . $dbname . ".flanking_sequence." . $pid . "." . $process . ".txt";
  open(SEQ,">",$file) or die ("Could not open sequence dumpfile for writing");
  
  #create the files to send to parallelize
  while (my $row = $sth->fetch()){
      $count++; #counting the total number of entries in the table
      $curr_variation_id = $row->[0]; #get the current variation_id

      if ($curr_variation_id ne $previous_variation_id){
		if((($sub_sequences * $process) <= $count) && ($process < $num_processes)){ #need to write in a new file
		  
		  # Close the current file 
		  close(SEQ);
		  
		  #ÊSend this job off to the farm
		  my $log_out = $TMP_DIR . "/output_flanking." . $$ . "." . $process . ".txt";		  
		  my $call = qq{bsub -q $bsub_queue_name -o $log_out $PERLBIN $script -tmpdir $TMP_DIR -tmpfile $TMP_FILE.$process -num_processes $num_processes -status_file $flanking_status_file -file $process -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname };
          $call .= "-cpass $cpass " if ($cpass);
          $call .= "-cport $cport " if ($cport);
          $call .= "-vpass $vpass " if ($vpass);
          $call .= "-pid $pid ";
          system($call);
      
		  $process++;
		  
		  # Open a new file
		  $file = $TMP_DIR . "/" . $dbname . ".flanking_sequence." . $pid . "." . $process . ".txt";
		  open(SEQ,">",$file) or die ("Could not open sequence dumpfile for writing");
		}
      }
	  
      my @a = map {defined($_) ? $_ : '\N'} @$row;
      
      print SEQ join("\t",@a) . "\n";
      
      $previous_variation_id = $curr_variation_id;
  }
  
	  # Close the current file 
	  close(SEQ);
	  
	  #ÊSend this job off to the farm
	  my $log_out = $TMP_DIR . "/output_flanking." . $pid . "." . $process . ".txt";		  
	  my $call = qq{bsub -q $bsub_queue_name -o $log_out $PERLBIN $script -tmpdir $TMP_DIR -tmpfile $TMP_FILE.$process -num_processes $num_processes -status_file $flanking_status_file -file $process -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname };
      $call .= "-cpass $cpass " if ($cpass);
      $call .= "-cport $cport " if ($cport);
      $call .= "-vpass $vpass " if ($vpass);
      $call .= "-pid $pid ";
      system($call);
  
  $sth->finish();  
  
}


#will fill in the transcriptoutput_transcript_0_6030 variation table. Has to wait until the variation feature table has been filled up. Then, divide the number of entries
# by the number of processes to make the subprocesses
sub parallel_transcript_variation{
  my $dbVar = shift;
  my $script = shift;
  my $registry_file = shift;
  my $bsub_queue_name = shift;
  
  my $total_process = 0;

  my $length_slices = 0; #number of entries in the variation_feature table
  my $call;
  my $transcript_status_file = "transcript_status_file_$$\.log";
  
  #first, create the log file for the transcript_variation table
  open STATUS, ">$TMP_DIR/$transcript_status_file"
  or throw("Could not open tmp file: $TMP_DIR/$transcript_status_file\n"); 
  close STATUS;
  
  #then, calculate the rows for each subprocess
  my $sa = $dbCore->get_SliceAdaptor();
  my $inc_non_ref = 1;
  my $slices = $sa->fetch_all('toplevel', undef, $inc_non_ref);
  
  #order the slices by name
  my @slices_ordered = sort {$a->seq_region_name cmp $b->seq_region_name }  @{$slices};
  
  # assumes that variation features have already been pushed to toplevel
  foreach my $slice (@slices_ordered) {
	$length_slices += $slice->length;
  }

  my $batchid;
  if( $SCHEDULER eq 'PBS' ){ # For PBS/qsub;
	# Set a job that depends on the completion of all other jobs
	$batchid = submit_qsub(
	  q( echo 'sleep 1' ),
	  [ -W => "depend=on:$num_processes", -N => "ptv_waiter" ]
	);
  }
  
  #I must add up the length of all the slices to find the limit for each 
  my $sub_slice = int($length_slices / $num_processes);
  my $slice_max; #the number of slices in the chunk
  my $slice_min = 0; #first slice in the chunk
  
  $bsub_queue_name ||= 'long';
  
  for (my $i = 0; $i < $num_processes ; $i++){
	$length_slices = 0;
	$slice_max = 0;
	
	foreach my $slice (@slices_ordered){
	  if (($length_slices + $slice->length )<= ($i+1)*$sub_slice){
		$length_slices += $slice->length;
		$slice_max++;	       
	  }
	  
	  else{last;}
	}
	
	$limit = $slice_min . "," . ($slice_max-$slice_min) if ($i+1 < $num_processes or $num_processes==1);
	$limit = $slice_min . "," . (scalar(@slices_ordered)-$slice_min) 
	
	if ($i+1 == $num_processes and $num_processes != 1); #the last slice, get the left slices
	
	my $command = (
	  "$PERLBIN $script ".
	  "-registry_file $registry_file ".
	  "-species $species ".
	  "-limit $limit ".
	  "-tmpdir $TMP_DIR ".
	  "-tmpfile $TMP_FILE ".
	  "-num_processes $num_processes ".
	  "-status_file $transcript_status_file "
	);
	
	my $outfile = "$TMP_DIR/output_transcript_$i\_$$.txt";
	my $errfile = "$TMP_DIR/error_transcript_$i\_$$.txt";
	my $jobname = "ptv_job_$i";
	
	if( $SCHEDULER eq 'PBS' ){ # Use specified PBS/qsub
	  my $jobid = submit_qsub(
		$command,
		[
		  '-W' => "depend=beforeany:$batchid",
		  '-o' => $outfile,
		  '-e' => $errfile,
		  '-v' => 'PERL5LIB',
		  '-N' => $jobname
		]
	  );
	}
	
	else{ # Use default LSF/bsub
	  $call =
		"bsub -a normal -o $outfile -q $bsub_queue_name ".
		($species =~ /human|homo/ ? "-M8000000 -R'select[mem>8000] rusage[mem=8000]'" : "").
		" $command";
	  system($call);
	}
	
	$slice_min = $slice_max;
  }
}

##use genotype method to change alleles in genotype table to forward strand, this is for mouse only in order to merge with sanger called SNPs

sub reverse_things {

  my $dbVar = shift;
  
  # create a table of things to be reversed
  my $hap_line;
  if ($hap_id_string =~ /\([\,\d]+\)/) {
    $hap_line = "AND seq_region_id not in $hap_id_string";
  }
  else {
    $hap_line ='';
  }
  
  $dbVar->do(qq{
	CREATE TABLE IF NOT EXISTS variations_to_reverse
	SELECT vf.variation_id FROM variation_feature vf LEFT JOIN failed_variation fv
	ON fv.variation_id = vf.variation_id
	WHERE map_weight = 1 AND seq_region_strand = -1
	AND fv.variation_id IS NULL
	$hap_line;
  });
  
  $dbVar->do(qq{create index variation_idx on variations_to_reverse(variation_id);});
  
  reverse_genotype($dbVar);#run genotype/allele first, then run variation_feature and flanking_sequence
  reverse_variation_feature($dbVar);
  reverse_flanking_sequence($dbVar);
  
  $dbVar->do(qq{update variation v, variations_to_reverse vr set v.flipped = 1 where vr.variation_id = v.variation_id;});
  $dbVar->do(qq{update variation set flipped = 0 where flipped is null;});
}

sub reverse_genotype{
  my $dbVar = shift;
  my $buffer = {}; #hash containing all the files to be parallelized
  my $dbsnp = 0; #for non-dbsnp data
  
  my ($source_id,$name,$validation_status,$ancestral_allele,$flipped,$count,$population_genotype_id,$variation_id,$ss_id,
      $allele_id,$allele,$allele_1,$allele_2,$frequency,$sample_id,$allele_string,$seq_region_id,$class_so_id,$is_somatic,$minor_allele,$minor_allele_freq,$minor_allele_count,$clin_sig);
  
  my (%rec_seq_region_ids,%rec_table_name,$table_name);

  my $hap_line;
  if ($hap_id_string =~ /\([\,\d]+\)/) {
    $hap_line = "AND seq_region_id not in $hap_id_string";
  }
  else {
    $hap_line ='';
  }
  
  my @table_list = ("variation","individual_genotype_multiple_bp","population_genotype","allele");

## this section adds the genotype tables
## we are not doing this anymore as compress_genotypes.pl will cope with strands
#  if($species =~ /hum|homo/i or $vdbname =~/hum|homo/i) {
#	my $sth1 = $dbVar->prepare(qq{SHOW TABLES LIKE 'tmp_individual_genotype_single_bp_SubInd%'});
#	$sth1->execute;
#	
#	while(my ($table_name) = $sth1->fetchrow_array()) {
#	  push @table_list, $table_name;
#	}
#  }
#  
#  else {
#	push @table_list, "tmp_individual_genotype_single_bp";
#  }

  foreach my $table (@table_list) {
  #we only change things that have map_weight=1 that is vf.map_weight=1 and seq_region_id not in (hap_id_string)
    debug("processing table $table");
    if ($table !~ /single/) {
      open OUT, ">$TMP_DIR/$table\.$vdbname" or die "can't open file $table\_out : $!";
    }
    elsif ($table =~ /single/) {
      system("rm -f $TMP_DIR/$table\_del");
      system("rm -f $TMP_DIR/$table\.import");
    }

    my $sth = $dbVar->prepare(qq{
	  SELECT tg.*,vf.allele_string,vf.seq_region_id 
	  FROM $table tg, variation_feature vf, variations_to_reverse vtr
	  WHERE vtr.variation_id=tg.variation_id
      AND vtr.variation_id=vf.variation_id
	  AND vf.map_weight=1
	  $hap_line 
	}, {mysql_use_result=>1} );  

    $sth->execute();
    
	# bind columns
	if ($table =~ /pop/i) {
      $sth->bind_columns(\$population_genotype_id,\$variation_id,\$ss_id,\$allele_1,\$allele_2,\$frequency,\$sample_id,\$count,\$allele_string,\$seq_region_id);
    }
    elsif($table =~ /allele/i) {
      $sth->bind_columns(\$allele_id,\$variation_id,\$ss_id,\$allele,\$frequency,\$sample_id,\$count,\$allele_string,\$seq_region_id);
    }
    elsif($table =~ /variation/) {
      $sth->bind_columns(\$allele_id,\$source_id,\$name,\$validation_status,\$allele,\$flipped,\$class_so_id,\$is_somatic,\$minor_allele,\$minor_allele_freq,\$minor_allele_count,\$clin_sig,\$allele_string,\$seq_region_id);
    }
    else {
      $sth->bind_columns(\$variation_id,\$ss_id,\$allele_1,\$allele_2,\$sample_id,\$allele_string,\$seq_region_id);
    }
	
    LINE : while ($sth->fetch()){
      next if (!$allele and $table =~/variation|allele/);
      my ($old_allele,$old_allele_1,$old_allele_2);
      my @alleles = split /\//, $allele_string;
      foreach my $a (@alleles) {
        if ($a =~ /\(|\)/g) {
          $a =~ /\((\S+)\)(\d+)/;
          if ($1 and $2) {
            $a = $1;
            reverse_comp(\$a);
            $a = "($a)$2";
          }
        }
        else {
          reverse_comp(\$a);
        }
      }  
      my $new_allele_string = join "/", @alleles;
      
      if ($table =~ /allele|variation/) {
		next if ($allele =~ /\-|\+/);
		next if ($allele =~ /N/);
		next if ($allele =~ /^\(.*\)$/ig);
		$old_allele = $allele;
		my ($al,$nu,$allele_new);
		
		if ($allele =~ /\(|\)/g) {            
		  $allele =~ /\((\S+)\)(\d+)/;
		  if ($1 and $2) {
			$al = $1;
			$nu = $2;
		  }
		}
		
		next LINE if ($al and $nu and $new_allele_string =~ /$al/); #not change it if new_allele_string already match allele
		
		if ($new_allele_string !~ /$allele/ or ($new_allele_string =~ /$allele/ and $allele_string =~ /$allele/)) {
		  if ($allele =~ /\(|\)/g) {
			$allele =~ /\((\S+)\)(\d+)/; 
			
			if ($1 and $2) {
			  $allele = $1;
			  reverse_comp(\$allele);
			  $allele_new = "($allele)$2";
			}
		  }
		  else {
			reverse_comp(\$allele);
			$allele_new = $allele;
		  }
		  
		  if ($new_allele_string =~ /$allele/) {#use allele here to avoid case like ($allele)num
			if ($table =~ /allele/) {
			  print OUT "update ignore $table set allele = \"$allele_new\" where allele_id=$allele_id;\n";
			}
			elsif ($table =~ /variation/) {
			  print OUT "update ignore $table set ancestral_allele = \"$allele_new\" where variation_id=$allele_id;\n";
			}
		  }
		}
      }
	  
      else { 
		next if ($allele_1 =~ /\-/ and $allele_2 =~ /\-/);
		next if ($allele_1 =~/\+/ or $allele_2 =~/\+/);
		next if ($allele_1 =~ /^\(.*\)$/ or $allele_2 =~ /^\(.*\)$/);
		next if ($allele_1 =~/[^-acgtnACGTN]/ or $allele_2 =~/[^-acgtnACGTN]/);
		next if ($allele_1 =~/homozygous|indeterminate|SEQUENCE/ig or $allele_2 =~/homozygous|indeterminate|SEQUENCE/ig);
		$old_allele_1 = $allele_1;
		$old_allele_2 = $allele_2;
		my ($al,$nu,$new_allele_1,$new_allele_2);
		
		foreach my $allele ($allele_1,$allele_2) {
		  if ($allele =~ /\(|\)/g) {
			$allele =~ /\((\S+)\)(\d+)/;
			if ($1 and $2) {
			  $al = $1;
			  $nu = $2;
			}
		  } 
		  next LINE if ($al and $nu and $new_allele_string =~ /$al/);
		}
		
		if (($new_allele_string !~ /$allele_1|$allele_2/) or ($new_allele_string =~ /$allele_1|$allele_2/ and $allele_string =~ /$allele_1|$allele_2/)) {
		  if ($allele_1 =~ /\(|\)/g) {
			$allele_1 =~ /\((\S+)\)(\d+)/; 
			if ($1 and $2) {
			  $allele_1 = $1;
			  reverse_comp(\$allele_1);
			  $new_allele_1 = "($allele_1)$2";
			}
		  }
		  else {
			reverse_comp(\$allele_1);
			$new_allele_1 = $allele_1;
		  }
		  
		  if ($allele_2 =~ /\(|\)/g) {
			$allele_2 =~ /\((\S+)\)(\d+)/; 
			if ($1 and $2) {
			  $allele_2 = $1;
			  reverse_comp(\$allele_2);
			  $new_allele_2 = "($allele_2)$2";
			}
		  }
		  else {
			reverse_comp(\$allele_2) if $allele_2 !~/r|y/i;
			$new_allele_2 = $allele_2;
		  }
		  
		  if (($new_allele_string =~ /$allele_1/ and $new_allele_string =~ /$allele_2/) or ($new_allele_string =~ /$allele_1|$allele_2/ and ($allele_1 =~/r|y|-/i or $allele_2 =~ /r|y|-/i))) {
			if ($table =~ /pop/i) {
			  print OUT
				"update ignore $table set allele_1 = \"$new_allele_1\", ".
				"allele_2 = \"$new_allele_2\" ".
				"where population_genotype_id = $population_genotype_id;\n";
			}
			elsif ($table =~ /single/) {
			  $table_name = $table;
			  print_buffered($buffer, "$TMP_DIR/$table_name\_del", 
			  "DELETE FROM $table_name WHERE variation_id=$variation_id and subsnp_id=$ss_id and sample_id=$sample_id;\n");
			  print_buffered( $buffer, "$TMP_DIR/$table_name\.import",
			  join ("\t",$variation_id,$ss_id,$new_allele_1,$new_allele_2,$sample_id)."\n");
			}
			else {
			  print OUT
				"update ignore $table set allele_1 = \"$new_allele_1\", ".
				"allele_2 = \"$new_allele_2\" where variation_id=$variation_id ".
				"and subsnp_id=$ss_id and  sample_id=$sample_id ".
				"and allele_1 = \"$old_allele_1\" ".
				"and allele_2 = \"$old_allele_2\";\n";
			}
		  }
		}
      }
    }
    print_buffered( $buffer);
    #system("mysql  -u$vuser -p$vpass  -h $vhost $vdbname <$TMP_DIR/$table\_out");
  }
}

sub reverse_variation_feature{
  my $dbVar = shift;

  debug("processing reverse variation_feature strand");

  my ($variation_feature_id,$allele_string,$strand);

  my $hap_line;
  if ($hap_id_string =~ /\([\,\d]+\)/) {
    $hap_line = "AND seq_region_id not in $hap_id_string";
  }
  else {
    $hap_line ='';
  }

  foreach my $table ("variation_feature") {
    $dbVar->do(qq{CREATE TABLE variation_feature_before_re LIKE variation_feature});
    $dbVar->do(qq{INSERT INTO variation_feature_before_re SELECT * FROM variation_feature});
    my $sth;
    $sth = $dbVar->prepare(qq{
	  select vf.variation_feature_id,vf.allele_string,vf.seq_region_strand
	  from $table vf, variations_to_reverse vtr
      where vf.variation_id = vtr.variation_id
	});
  

    $sth->execute();
    $sth->bind_columns(\$variation_feature_id,\$allele_string,\$strand);

    while ($sth->fetch()){
      my @alleles = split /\//,$allele_string;
      foreach my $allele (@alleles) {
		next if $allele =~ /^\(.*\)$/;
		next if $allele =~ /CNV|INS|DEL/;
		if ($allele =~ /\(|\)/g) {
		  $allele =~ /\((\S+)\)(\d+)/; 
		  if ($1 and $2) {
			$allele = $1;
			reverse_comp(\$allele);
			$allele = "($allele)$2";
		  }
		}
		else {
		  reverse_comp(\$allele);
		}
      }
      
      # flip strand
      $strand *= -1;
      
      my $new_allele_string = join "/",@alleles;
      $dbVar->do(qq{update $table set allele_string = "$new_allele_string", seq_region_strand =$strand where variation_feature_id=$variation_feature_id});
    }
  }
}

sub reverse_flanking_sequence {

  my $dbVar = shift;

  debug("Processing flanking_sequence reverse strand");

  my $hap_line;
  if ($hap_id_string =~ /\([\,\d]+\)/) {
    $hap_line = "AND vf.seq_region_id not in $hap_id_string";
  }
  else {
    $hap_line ='';
  }

  $dbVar->do(qq{CREATE TABLE flanking_sequence_before_re LIKE flanking_sequence});
  $dbVar->do(qq{INSERT INTO flanking_sequence_before_re SELECT * FROM flanking_sequence});

  foreach my $table ("flanking_sequence") {
    my $sth=$dbVar->prepare(qq{
	  select fl.* from flanking_sequence fl, variations_to_reverse vtr
	  where fl.variation_id = vtr.variation_id
	});
	
    $sth->execute();
    my ($variation_id,$up_seq,$down_seq,$up_seq_start,$up_seq_end,$down_seq_start,$down_seq_end,$seq_region_id,$seq_region_strand);

    $sth->bind_columns(\$variation_id,\$up_seq,\$down_seq,\$up_seq_start,\$up_seq_end,\$down_seq_start,\$down_seq_end,\$seq_region_id,\$seq_region_strand);
    while($sth->fetch()) {
	  #print "$variation_id,$up_seq,$down_seq,$up_seq_start,$up_seq_end,$down_seq_start,$down_seq_end,$seq_region_id,$seq_region_strand\n";
	  if ($up_seq and $down_seq) {
		($up_seq, $down_seq) = ($down_seq, $up_seq);
		reverse_comp(\$up_seq);
		reverse_comp(\$down_seq);
		$up_seq_start = $up_seq_end = $down_seq_start = $down_seq_end = '\N';
	  }
	  elsif (! $up_seq and ! $down_seq) {
		my $tmp_seq_start = $up_seq_start;
		my $tmp_seq_end = $up_seq_end;
		($up_seq_start, $up_seq_end) = ($down_seq_start, $down_seq_end);
		($down_seq_start, $down_seq_end) = ($tmp_seq_start, $tmp_seq_end);
		$up_seq = $down_seq = '\N';
	  }
	  elsif ($up_seq and ! $down_seq) {
		$down_seq = $up_seq;
		reverse_comp(\$down_seq);
		$up_seq = '\N';
		($up_seq_start, $up_seq_end) = ($down_seq_start, $down_seq_end);
		$down_seq_start = '\N';
		$down_seq_end = '\N';
	  }
	  elsif (! $up_seq and $down_seq) {
		$up_seq = $down_seq;
		reverse_comp(\$up_seq);
		$down_seq = '\N';
		($down_seq_start, $down_seq_end) = ($up_seq_start, $up_seq_end);
		$up_seq_start = '\N';
		$up_seq_end = '\N';
	  }
	  $seq_region_strand = 1;
	  
	  $dbVar->do(qq{update $table set up_seq = "$up_seq", down_seq = "$down_seq", up_seq_region_start = $up_seq_start, up_seq_region_end = $up_seq_end, down_seq_region_start=$down_seq_start, down_seq_region_end=$down_seq_end, seq_region_strand=1  where variation_id = $variation_id and seq_region_id = $seq_region_id and seq_region_strand = -1});
    }
	
    $dbVar->do(qq{update $table set up_seq = null where up_seq = 'N'});
    $dbVar->do(qq{update $table set down_seq = null where down_seq = 'N'});
  }
}

sub merge_ensembl_snps {

  #this one is used to merge ENSEMBL SNPs with dbSNP SNPs and other SNPs
  my $dbVar = shift;
  my $new_source_id = shift ;
  
  #ÊOptional listref to a list of source_ids that the new source_id will be merged with 
  my $into_source_id = shift;
  $into_source_id ||= [];
  
  my $hap_line;
  if ($hap_id_string =~ /\([\,\d]+\)/) {
      $hap_line = "AND vf1.seq_region_id not in $hap_id_string
                   AND vf2.seq_region_id not in $hap_id_string";
  }
  else {
    $hap_line ='';
  }

  debug("Note there are venter and celera variation_name are shared, so if celera variation name is in first, then venter name can't be imported into variation,i.e variation_id for those variation with same name are not in variation table, but in all other tables, should put those names in variation_synonym table with source_id=venter, then delete/change those variation_id in other tables see table vf_share_name_same_vf in both celera and venter databases");
  die ("You need to provide the new source_id that you want the data to merge") if !$new_source_id;

  if (! $new_source_id) {
    my $new_source_id_ref = $dbVar->selectall_arrayref(qq{select max(source_id) from source});
    $new_source_id = $new_source_id_ref->[0][0];
  }
  # The default behaviour is to merge with source_ids less than the new source id, so get a list of all source_ids less than the new id from db.
  unless (scalar(@{$into_source_id})) {
      $into_source_id = $dbVar->selectcol_arrayref(qq{SELECT source_id FROM source WHERE source_id < $new_source_id});
  }
  
  debug("new source_id is $new_source_id and it will be merged with source_ids " . join(",",@{$into_source_id}) . ", is this correct???\n");
  sleep(20);

  #ÊCreate a condition on the source_ids to merge with
  my $source_id_constraint = " (vf1.source_id = $new_source_id AND vf2.source_id IN (" . join(",",@{$into_source_id}) . ")) ";
  
  debug("Creat table tmp_ids...");#check for same position,different alleles as well
  $dbVar->do(qq{
	CREATE table tmp_ids_$new_source_id
	SELECT vf1.variation_id as variation_id1, vf1.variation_name as name1,
	vf2.variation_id as variation_id2,vf2.variation_name as name2 
	FROM (variation_feature vf1, variation_feature vf2)
	WHERE vf1.seq_region_id=vf2.seq_region_id 
	AND vf1.seq_region_start = vf2.seq_region_start 
	AND vf1.seq_region_end = vf2.seq_region_end
	AND vf1.seq_region_strand = vf2.seq_region_strand
	AND vf1.map_weight=1 and vf2.map_weight=1
	$hap_line
	AND vf1.variation_id > vf2.variation_id
	AND $source_id_constraint 
	AND vf1.somatic = vf2.somatic
 });

  $dbVar->do(qq{alter table tmp_ids_$new_source_id add index variation_idx1(variation_id1), add index variation_idx2(variation_id2)});
  
  #die "FINITO\n";

  debug("Change allele_string if different...");
  my ($variation_id1,$name1,$variation_id2,$name2,$allele_string1,$allele_string2,$variation_feature_id2);
  $dbVar->do(qq{
	CREATE table allele_string_diff_$new_source_id
	SELECT t.variation_id1, t.name1, SUBSTR(vf1.allele_string,1,25000) as allele_string1,
	t.variation_id2, t.name2, SUBSTR(vf2.allele_string,1,25000) as allele_string2,
	vf2.variation_feature_id as variation_feature_id2
	FROM tmp_ids_$new_source_id t, variation_feature vf1, variation_feature vf2
	WHERE t.variation_id1 =vf1.variation_id 
	AND t.variation_id2 =vf2.variation_id 
	AND vf2.map_weight=1 ##vf2 is dbSNP
	AND vf1.map_weight=1
	AND vf1.allele_string != 'N/A'#cnv has 'N' in allele_string, not merge this
	AND vf2.allele_string != 'N/A'#cnv has 'N' in allele_string, not merge this
	$hap_line
	AND vf1.allele_string != vf2.allele_string
 });

  my $sth = $dbVar->prepare(qq{
	SELECT variation_id1,name1,allele_string1,variation_id2,name2,
	allele_string2,variation_feature_id2
	FROM allele_string_diff_$new_source_id
  });
  $sth->execute();
  $sth->bind_columns(\$variation_id1,\$name1,\$allele_string1,\$variation_id2,\$name2,\$allele_string2,\$variation_feature_id2);

  while($sth->fetch()) {
    my @alleles1 = split /\//, $allele_string1; #ensembl allele_string
    my @alleles2 = split /\//, $allele_string2; #dbSNP allele_string

    my %rec_allele = map {$_,1} @alleles2 ;
    my $n=0;
    while ( $n < @alleles1) {
      if (! $rec_allele{$alleles1[$n]} and $alleles1[$n] !~ /MUTATION|CNV/) {#allele not seen and is a real allele
		$allele_string2 .= "/$alleles1[$n]"; #in this way to make sure the last one is from deleted allele_string if it's not in allele_string already
      }
      $n++;
    }
    my $allele_string = $allele_string2;
    $dbVar->do(qq{
	  UPDATE variation_feature set allele_string = \"$allele_string\" 
	  WHERE variation_id=$variation_id2
	  AND variation_feature_id=$variation_feature_id2
	});
  }

  debug("Inserting into variation_synonym table...");
  $dbVar->do(qq{
	INSERT IGNORE INTO variation_synonym (variation_id,source_id,name) 
	SELECT t.variation_id2,vf.source_id,vf.variation_name #dbSNP varid with ensembl-snp source and name
	FROM  variation_feature vf, tmp_ids_$new_source_id t 
	WHERE t.variation_id1=vf.variation_id
  });

  debug("Deleting from tables ...");
  foreach my $table(qw(variation variation_feature flanking_sequence)) {
	debug("Deleting from $table");
	$dbVar->do(qq{DELETE from v using $table v, tmp_ids_$new_source_id t where t.variation_id1=v.variation_id});
  }

  $dbVar->do(qq{
	DELETE from tv using transcript_variation tv 
    LEFT JOIN variation_feature vf on tv.variation_feature_id=vf.variation_feature_id 
    WHERE vf.variation_feature_id is null
  });
  
  

  # get rid of duplicated middle variation_id/names
  $dbVar->do(qq{
	CREATE TABLE IF NOT EXISTS tmp_ids_$new_source_id\_final
	SELECT variation_id1,name1,min(variation_id2) as variation_id2, 'placeholder' as name2
	FROM tmp_ids_$new_source_id
	GROUP BY variation_id1
  });
  $dbVar->do(qq{alter table tmp_ids_$new_source_id\_final add index name2(name2)});
  $dbVar->do(qq{alter table tmp_ids_$new_source_id\_final add index variation_idx1(variation_id1), add index variation_idx2(variation_id2)});
  $dbVar->do(qq{
	UPDATE tmp_ids_$new_source_id\_final f, tmp_ids_$new_source_id o
	SET f.name2 = o.name2
	WHERE f.variation_id2=o.variation_id2
  });

  debug("Updating tables using tmp_ids_$new_source_id\_final ...");
  $dbVar->do(qq{UPDATE variation_synonym vs, tmp_ids_$new_source_id\_final t set vs.variation_id=t.variation_id2 
      	  WHERE vs.variation_id = t.variation_id1});
  #if t.variation_id2 is already exist in table vsv, then the update will be ignored, these ignored variation_id1 needs to be deleted
  $dbVar->do(qq{DELETE vsv from variation_set_variation vsv
	        LEFT JOIN variation v ON vsv.variation_id=v.variation_id
		WHERE v.variation_id is null});

  debug("Updating tables...");
  
  foreach my $table(qw(allele population_genotype variation_set_variation variation_annotation individual_genotype_multiple_bp tmp_individual_genotype_single_bp variation_synonym failed_variation)) {
	debug("Updating $table");
	$dbVar->do(qq{UPDATE IGNORE $table a, tmp_ids_$new_source_id\_final t SET a.variation_id = t.variation_id2 where t.variation_id1 = a.variation_id});
  }
  
  # when insert ignore, if the entry is already exist, it will do nothing,
  # i.e leave the one unchanged, so need to delete them
  foreach my $table (qw(tmp_individual_genotype_single_bp failed_variation)) {
    $dbVar->do(qq{
	   DELETE FROM vs USING $table vs, tmp_ids_$new_source_id\_final t 
	   WHERE vs.variation_id=t.variation_id1
    });
  }
  
  #$dbVar->do(qq{DROP table tmp_ids_$source_id});
}

sub read_updated_files {
  my $dbVar = shift;

  #read in updated files generated from reverse_things:
  
  my $command = "mysql -u$vuser -p$vpass -h$vhost $vdbname";
  my $command1 = "mysqlimport -L -u$vuser -p$vpass -h$vhost $vdbname";
  
  print "command is $command command1 is $command1\n";
  
  my %rec_pid1;

  opendir DIRENTRY, "$TMP_DIR" || die "Failed to open dir : $!";
  my @files = grep /\_del$/, readdir(DIRENTRY);

  my $num_process = 4;
  
  foreach my $file (@files) {
    debug("Starting read file $file...");
    my $kid;
    my $pid = fork;
    if (! defined $pid) {
      throw("Not possible to fork: $!\n");
    }
    elsif ($pid ==0) {
	  print "$command \< $TMP_DIR\/$file\n";
      system("$command < $TMP_DIR/$file");
      my $file1 = $file;
      $file1 =~ s/\_del/\.import/;
      print "file $file has done\n";
      system("mv $TMP_DIR/$file $TMP_DIR/$file\_done");
	  print "$command1 $TMP_DIR\/$file1\n";
      system("$command1 $TMP_DIR/$file1");
      system("mv $TMP_DIR/$file1 $TMP_DIR/$file1\_done");
      
      POSIX:_exit(0);
    }
    $rec_pid1{$pid}++;
    if (keys %rec_pid1 == $num_process){
      do{
        #the parent waits for the child to finish;
        $kid = waitpid(-1, WNOHANG);
        if ($kid > 0){
          delete $rec_pid1{$kid};
        }
      } until keys %rec_pid1 < $num_process;
    }
  }
  do{
    #the parent waits for the child to finish for the last three processes;
    my $kid = waitpid(-1, WNOHANG);  
    if ($kid > 0){
      delete $rec_pid1{$kid};
    }
  } until keys %rec_pid1 <1;

}

sub merge_multi_tables {

  my $dbVar = shift;
   
  foreach my $table("combine_feature","failed_flanking_qual","failed_gtype","flanking_qual","gtype_allele","snp_pos","tmp_individual_genotype_single_bp") {
    my $single_table_ref = $dbVar->selectall_arrayref(qq{SHOW tables like "$table%"});
    
    my @tables = map {$_->[0] } @$single_table_ref;
    my $table_names = join ',',@tables;
    #my $final_names = $table_names . ",tmp_individual_genotype_single_bp_last_insert";
    my $final_names = $table_names;
    print "dbname is $vdbname and table_names is $final_names\n";

    #$dbVar->do(qq{CREATE TABLE $table\_last_insert like $tables[0]});
    $dbVar->do(qq{DROP TABLE IF EXISTS $table});
    $dbVar->do(qq{
	  CREATE TABLE $table (
	  variation_id int not null, subsnp_id int(15) unsigned, 
	  allele_1 varchar(255),allele_2 varchar(255),sample_id int,
	  key variation_idx(variation_id),
	  key subsnp_idx(subsnp_id),
	  key sample_idx(sample_id)
	  ) ENGINE=MERGE UNION=($final_names) INSERT_METHOD=LAST
	});
  }
}

sub merge_rs_feature{
  ###this method is used to merge dbSNP SNPs to merge SNPs with same position and allele, also with map_weight=1, put the one with longest flanking_sequence into variation table, the other one into variation_synonym table
  ###not longest sequence any more, but keep smaller rs number
  my $dbVar = shift;

  my $hap_line;
  if ($hap_id_string =~ /\([\,\d]+\)/) {
    $hap_line = "AND vf1.seq_region_id not in $hap_id_string
                 AND vf2.seq_region_id not in $hap_id_string";
  }
  else {
    $hap_line ='';
  }

  debug("Creat table tmp_ids...");
  $dbVar->do(qq{
	CREATE table tmp_ids_rs
	SELECT vf1.variation_id as variation_id1, vf1.variation_name as name1,
	vf2.variation_id as variation_id2,vf2.variation_name as name2
	FROM variation_feature vf1, variation_feature vf2
	WHERE vf1.seq_region_id=vf2.seq_region_id 
	AND vf1.seq_region_start = vf2.seq_region_start 
	AND vf1.seq_region_end = vf2.seq_region_end
	AND vf1.seq_region_strand = vf2.seq_region_strand  #only merge same strand
	And vf1.map_weight=1 and vf2.map_weight=1
	AND vf1.source_id=1 and vf2.source_id=1
	AND vf1.allele_string not like '%-%'
	AND vf2.allele_string not like '%-%'
	AND vf1.allele_string like '%/%'
	AND vf2.allele_string like '%/%'
	$hap_line
	AND round(substring(vf1.variation_name,3)) > round(substring(vf2.variation_name,3))
	AND vf1.somatic = vf2.somatic
  });

  $dbVar->do(qq{alter table tmp_ids_rs add index variation_idx1(variation_id1), add index variation_idx2(variation_id2)});
  $dbVar->do(qq{alter table tmp_ids_rs add index name2(name2)}); #for later use
  
  # added to merge allele_string for dbSNP variation, if not merging, two sets of genotype data associated to different variation_ids but with same location, this will confuse web display using table compressed_genotype_single_bp 
  debug("Change allele_string if different..."); #added for merging allele_string for dbSNP variation
  
  # Modified the query to only select a substring of each allele_string since the total row length will otherwise exceed the maximum
  $dbVar->do(qq{
	CREATE table allele_string_diff_rs
	SELECT t.variation_id1, t.name1, SUBSTR(vf1.allele_string,1,25000) as allele_string1,
	t.variation_id2,t.name2, SUBSTR(vf2.allele_string,1,25000) as allele_string2,
	vf2.variation_feature_id as variation_feature_id2
	FROM tmp_ids_rs t, variation_feature vf1, variation_feature vf2
	WHERE t.variation_id1 =vf1.variation_id 
	AND t.variation_id2 =vf2.variation_id 
	AND vf1.allele_string != vf2.allele_string
 });

  my ($variation_id1,$name1,$variation_id2,$name2,$allele_string1,$allele_string2,$variation_feature_id2);
  my $sth = $dbVar->prepare(qq{SELECT variation_id1,name1,allele_string1,variation_id2,name2,
                                      allele_string2,variation_feature_id2
                               FROM allele_string_diff_rs});
  $sth->execute();
  $sth->bind_columns(\$variation_id1,\$name1,\$allele_string1,\$variation_id2,\$name2,\$allele_string2,\$variation_feature_id2);
  
  my %rec_allele_string2; #to hold cumulated allele_string
  
  while($sth->fetch()) {
    #check if the allele_string already changed, take the changed one if it was changed
    $allele_string1 = ($rec_allele_string2{$variation_id1}) ? $rec_allele_string2{$variation_id1} : $allele_string1;
    $allele_string2 = ($rec_allele_string2{$variation_id2}) ? $rec_allele_string2{$variation_id2} : $allele_string2;
    
    my @alleles1 = split /\//, $allele_string1; #dbSNP allele_string with bigger rsnumber
    my @alleles2 = split /\//, $allele_string2; #dbSNP allele_string with small rsnumber

    my %rec_allele = map {$_,1} @alleles2 ;
    my $n=0;
    while ( $n < @alleles1) {
      if (! $rec_allele{$alleles1[$n]}) {
		$allele_string2 .= "/$alleles1[$n]"; #in this way to make sure the last one is from deleted allele_string if it's not in allele_string already
      }
      $n++;
    }
    my $allele_string = $allele_string2;
    
    $rec_allele_string2{$variation_id2} = $allele_string;

    #update allele_string both in variation_feature and allele_string_diff_rs
    $dbVar->do(qq{
	  UPDATE variation_feature set allele_string = \"$allele_string\" 
	  WHERE variation_feature_id = $variation_feature_id2
	});
  }
  #upto here -- added to merge allele_string for dbSNP variation
  


  debug("Inserting into variation_synonym table...");
  $dbVar->do(qq{
	INSERT IGNORE INTO variation_synonym (variation_id,source_id,name) 
	SELECT t.variation_id2,v.source_id,v.name 
	FROM  variation v, tmp_ids_rs t
	WHERE t.variation_id1=v.variation_id
  });

  debug("Deleting from tables ...");
  
  foreach my $table(qw(variation variation_feature flanking_sequence)) {
	debug("Deleting from $table");
	$dbVar->do(qq{DELETE from v using $table v, tmp_ids_rs t where t.variation_id1=v.variation_id});
  }

  $dbVar->do(qq{
	DELETE from tv using transcript_variation tv 
    LEFT JOIN variation_feature vf on tv.variation_feature_id=vf.variation_feature_id 
    WHERE vf.variation_feature_id is null
  });

  #get rid of duplicated middle variation_id/names
  $dbVar->do(qq{
	CREATE TABLE IF NOT EXISTS tmp_ids_rs_final
	SELECT variation_id1,name1,variation_id2,concat('rs',min(substring(name2,3))) as name2
	FROM tmp_ids_rs
	GROUP BY variation_id1
  });
  $dbVar->do(qq{alter table tmp_ids_rs_final add index name2(name2)});
  $dbVar->do(qq{alter table tmp_ids_rs_final add index variation_idx1(variation_id1), add index variation_idx2(variation_id2)});
  $dbVar->do(qq{
	UPDATE tmp_ids_rs_final f, tmp_ids_rs o
	SET f.variation_id2 = o.variation_id2
	WHERE f.name2=o.name2
  });

  debug("Updating tables using tmp_ids_rs_final ...");
  $dbVar->do(qq{UPDATE variation_synonym vs, tmp_ids_rs_final t set vs.variation_id=t.variation_id2 
      	  WHERE vs.variation_id = t.variation_id1});
  #if t.variation_id2 is already exist in table vsv, then the update will be ignored, these ignored variation_id1 needs to be deleted
  $dbVar->do(qq{DELETE vsv from variation_set_variation vsv
	        LEFT JOIN variation v ON vsv.variation_id=v.variation_id
		WHERE v.variation_id is null});
  
  ## we need to deal with the situation where one of the variants has been flipped and the other hasn't
  
  # first create a flipper table to flip the genotypes
  $dbVar->do(qq{CREATE TABLE flipper (a char(1) NOT NULL, b char(1) NOT NULL);});
  $dbVar->do(qq{INSERT INTO flipper VALUES('A','T');});
  $dbVar->do(qq{INSERT INTO flipper VALUES('C','G');});
  $dbVar->do(qq{INSERT INTO flipper VALUES('G','C');});
  $dbVar->do(qq{INSERT INTO flipper VALUES('T','A');});
  
  debug("Flipping genotypes in cases where variation to be merged has different flipped state");
  
  # flip the genotypes where the flipped column does not match
  # otherwise we will get some genotypes on + and some on -
  # when we update variation_id below
  $dbVar->do(qq{
	UPDATE tmp_individual_genotype_single_bp tg, tmp_ids_rs_final t, variation v1, variation v2, flipper f1, flipper f2
	SET tg.allele_1 = f1.b, tg.allele_2 = f2.b
	WHERE tg.variation_id = t.variation_id1
	AND t.variation_id1 = v1.variation_id
	AND t.variation_id2 = v2.variation_id
	AND v1.flipped != v2.flipped
	AND tg.allele_1 = f1.a
	AND tg.allele_2 = f2.a
  });
  
  # drop the flipper table
  $dbVar->do(qq{DROP TABLE flipper;});
  
  #ÊFor the tables that need updating and have primary keys, create temp tables rather than join directly
  foreach my $table (qw(variation_annotation allele population_genotype failed_variation)) {
    debug("Updating $table");
	
	#ÊDrop the temp table if it already exists
	my $stmt = qq{
	    DROP TABLE IF EXISTS
	       tmp_merge
	};
	$dbVar->do($stmt);
	
	#ÊCreate the temp table to map the entries that needs updating
	$stmt = qq{
	    CREATE TABLE
	       tmp_merge
	    SELECT
	       tbl.$table\_id
	       m.variation_id2 AS variation_id
	    FROM
	       tmp_ids_rs_final m JOIN
	       $table tbl ON (
	           tbl.variation_id = m.variation_id1
	       )
	};
	$dbVar->do($stmt);
	
	#ÊUpdate the table
	$stmt = qq{
	    UPDATE IGNORE
	       $table tbl,
	       tmp_merge m
	    SET
	       tbl.variation_id = m.variation_id
	    WHERE
	       tbl.$table\_id = m.$table\_id
	};
	$dbVar->do($stmt);
	
	#ÊDrop the temp table
	$stmt = qq{
	    DROP TABLE
	       tmp_merge
	};
	$dbVar->do($stmt);
  }
	
  #ÊFor the tables that don't have a straightforward primary key, join on variation_id and update
  foreach my $table(qw(variation_set_variation individual_genotype_multiple_bp tmp_individual_genotype_single_bp)) {
	debug("Updating $table");
	$dbVar->do(qq{
	  UPDATE IGNORE $table vs, tmp_ids_rs_final t
	  SET vs.variation_id=t.variation_id2 
	  WHERE vs.variation_id = t.variation_id1
	});
  }

  #ÊDelete rows from failed_variation that were ignored in the update statement due to duplication of keys
  $dbVar->do(qq{DELETE FROM fv USING tmp_ids_rs_final t JOIN failed_variation fv ON (fv.variation_id = t.variation_id1)});
   
  #$dbVar->do(qq{DROP table tmp_ids_rs}); #keep this table for future checking
}

sub remove_wrong_variations {

  my $dbVar = shift;
  
  # flag no variation mappings
  $dbVar->do(qq{
	INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) 
	SELECT v.variation_id, fd.failed_description_id
	FROM (failed_description fd, variation v) LEFT JOIN variation_feature vf
	ON v.variation_id=vf.variation_id
	WHERE vf.variation_id IS NULL
	AND fd.description = 'Variation does not map to the genome';
  });
  
  # flag non-variant
  $dbVar->do(qq{
	INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) 
	SELECT a.variation_id, fd.failed_description_id
	FROM  allele a, failed_description fd
	WHERE a.allele = 'NoVariation'
	AND fd.description = 'Loci with no observed variant alleles in dbSNP';
  });
  
  # flag N allele subsnps
  $dbVar->do(qq{
	INSERT INTO failed_allele (allele_id, failed_description_id)
	SELECT a.allele_id, fd.failed_description_id
	FROM allele a, failed_description fd
	WHERE a.allele = 'N'
	AND fd.description = 'Alleles contain ambiguity codes'
  });
    
  # now check for variations with >3 alleles
  my $sth = $dbVar->prepare(qq{
	SELECT variation_id, allele_string
	FROM  variation_feature
	WHERE length(allele_string) - length(REPLACE(allele_string,'/','')) > 2
	AND allele_string not like '%-%';
  });
  
  my ($variation_id, $allele_string);
  $sth->execute();
  $sth->bind_columns(\$variation_id, \$allele_string);
  
  my %failed;
  
  # check them all
  while($sth->fetch()) {
	my $ok = 0;
	
	# we only want to fail e.g. A/C/T/G, not e.g. A/CC/TT/G or A/-/T/G
	foreach my $a(split /\//, $allele_string) {
	  $ok = 1 if length($a) > 1 or $a eq "-";
	}
	
	$failed{$variation_id} = 1 unless ($ok);
  }
  
  $sth->finish();
  
  my $sth2 = $dbVar->prepare(qq{
	INSERT IGNORE INTO failed_variation (variation_id, failed_description_id)
	SELECT ?, failed_description_id
	FROM failed_description
	WHERE description = 'Variation has more than 3 different alleles';
  });
  
  foreach my $v_id(keys %failed) {
	$sth2->execute($v_id);
  }
}

sub check_seq_region_id {

  #This needs to be run when gene builder changed their gene set to check they still use same seq_region_id
  my $dbCore = shift;
  my $dbVar = shift;

  my $sth1 = $dbCore->dbc->db_handle->prepare(qq{
	SELECT seq_region_id 
	FROM seq_region_attrib sa, attrib_type at 
	WHERE sa.attrib_type_id=at.attrib_type_id 
	AND at.code='toplevel'
  });
  $sth1->execute();
  my %rec_seq_region_ids;
  while (my ($seq_region_id) = $sth1->fetchrow_array()) {
    $rec_seq_region_ids{$seq_region_id}=1;
  }

  my $id_res = $dbVar->selectall_arrayref(qq{SELECT distinct seq_region_id from variation_feature});
  my @seq_region_ids = map {$_->[0] } @$id_res;

  my $core_dbname = $dbCore->dbname();
  foreach my $seq_reg_id (@seq_region_ids) {
    if (!$rec_seq_region_ids{$seq_reg_id}) {
      debug("seq_region_id $seq_reg_id is not exist in $core_dbname");
    }
  }
}

sub make_allele_string_table {
#after mapping, needs post processing which needs allele_string table
#now in post-processing, it will check allele_string table exist? if not, will create one before run post_variation_feature check
  my $dbVar = shift;
  dumpSQL($dbVar,qq{SELECT variation_id,allele_string from variation_feature});
  open IN, "$TMP_DIR/$TMP_FILE";
  open OUT, ">$TMP_DIR/allele_string.txt";

  while(<IN>) {
    my ($variation_id,$allele_string) = split;
    my @alleles = split /\//,$allele_string;
    foreach my $allele (@alleles) {
      print OUT "$variation_id\t$allele\n";
    }
  }

  system("mv $TMP_DIR/allele_string.txt $TMP_DIR/$TMP_FILE");
  create_and_load($dbVar,"allele_string","variation_id i*","allele");
}

sub find_old_new_ids {

  my $old_new_ids = shift;
  my $table = shift;
  my $table_id = shift;
  my @colnames = @_;
  my $cols = join( ",", @colnames );
  
  my (@old_ids,$max_var_id,$new_id);
  
  if ($table =~ /^variation$/) {
    my $max_var_id_ref = $dbVar->selectall_arrayref(qq{select max(variation_id) from variation});
    $max_var_id = $max_var_id_ref->[0][0];
    $new_id = $max_var_id;
  }
  
  $dbVar->do(qq{CREATE TABLE old_new\_$table\_id (old\_$table_id int, new\_$table_id int, key old_idx(old\_$table_id), key new_idx(new\_$table_id))});
  my $sth = $dbVar->prepare(qq{select $table_id from $table\_merge});
  $sth->execute();
  while (my ($id) = $sth->fetchrow_array()) {
    push @old_ids,$id;
  } 
  foreach my $old_id (@old_ids) {
    #print "old_id is $old_id\n";
    if ($old_id) {
      if ($table =~ /^individual$/) {
        my $sample_id_ref = $dbVar->selectall_arrayref(qq{select new_sample_id from old_new_sample_id where old_sample_id=$old_id});
        $new_id = $sample_id_ref->[0][0];
        $dbVar->do(qq{INSERT INTO $table ($table_id,$cols) SELECT $new_id as sample_id,$cols from $table\_merge where $table_id = $old_id});
      }
      elsif ($table =~ /^sample$/) {
        $dbVar->do(qq{INSERT INTO $table ($cols) SELECT $cols from $table\_merge where $table_id = $old_id});
        $new_id = $dbVar->{'mysql_insertid'};
      }
      elsif ($table =~ /^variation$/) {
        $new_id++;
        $dbVar->do(qq{INSERT INTO $table ($table_id,$cols) SELECT $new_id as variation_id,$cols from $table\_merge where $table_id = $old_id});
      }
      $dbVar->do(qq{INSERT INTO old_new\_$table\_id (old\_$table_id,new\_$table_id) values ($old_id,$new_id)});
    }
  }
}

#will have to wait until the variation_feature has finished. Then, select all the genotype, and split the data into files (1 per population)
sub parallel_ld_populations {
    my $dbVar = shift;    
    my $bsub_queue_name = shift;

    my $call;

    my %seq_region; #hash containing the mapping between seq_region_id->name region
    my %alleles_variation = (); #will contain a record of the alleles in the variation. A will be the major, and a the minor. When more than 2 alleles
    my %genotype_information; #will contain all the genotype information to write in the file, if necessary
    #, the genotypes for that variation will be discarded
    my %regions; #will contain all the regions in the population and the number of genotypes in each one
    my $previous_variation_id = ''; #to know if it is a new variation and we can get the new alleles
    my $buffer = {}; #will contain a buffer where will be written all the LD information
    my %populations;
    my %genotypes_file; #foreach file, will contain the number of genotypes, so we can split it later

    #going to get the population_id for the HapMap and PerlEgene populations and a hash with the individuals that shouldn't
    #be present in the LD calculation
    my $pop_id;
    my $population_name;
    #get all populations to be tagged (HapMap and PerlEgen)
    my $sth = $dbVar->prepare(qq{SELECT s.sample_id, s.name
				     FROM population p, sample s
				     WHERE (s.name like 'PERLEGEN:AFD%'
				     OR s.name like 'CSHL-HAPMAP%')
				     AND s.sample_id = p.sample_id
				 });
    

    $sth->execute();
    $sth->bind_columns(\$pop_id,\$population_name);
    #get all the children that we do not want in the genotypes
    my $siblings = {}; # hash {$individual_id} ,where the individual is sibling of another one
    my @pops;
    while($sth->fetch){
	if($population_name =~ /CEU|YRI|MEX/){
	    &get_siblings($dbVar,$pop_id,$siblings);
	}
	push @pops, $pop_id;
    }
 
    my $in_str = " IN (" . join(',', @pops). ")";

    #necessary the order to know when we change variation. Not get genotypes with a NULL variation or map_weight > 1

    $sth = $dbVar->prepare
	(qq{
	    SELECT  STRAIGHT_JOIN ig.variation_id, vf.variation_feature_id, vf.seq_region_id, vf.seq_region_start, 
                          ig.sample_id, ig.allele_1, ig.allele_2, vf.seq_region_end, ip.population_sample_id
		    FROM  variation_feature vf FORCE INDEX(pos_idx), individual_genotype_single_bp ig, individual_population ip
		   WHERE  ig.variation_id = vf.variation_id

		    AND   ig.allele_2 IS NOT NULL
		    AND   vf.map_weight = 1
		    AND   ip.individual_sample_id = ig.sample_id
		    AND   ip.population_sample_id $in_str
		    ORDER BY  vf.seq_region_id,vf.seq_region_start}, {mysql_use_result => 1} );


    print "Time starting to dump data from database: ",scalar(localtime(time)),"\n";
    $sth->execute();

    my $dbname = $vdbname; #get the name of the database to create the file
    my ($variation_id, $variation_feature_id, $seq_region_id, $seq_region_start, 
	$individual_id, $allele_1,$allele_2,$seq_region_end,$population_id);


    $sth->bind_columns(\$variation_id, \$variation_feature_id, \$seq_region_id, \$seq_region_start, 
		       \$individual_id, \$allele_1,\$allele_2,\$seq_region_end,\$population_id); 
    while ($sth->fetch()){
	if ($previous_variation_id eq ''){
	    $previous_variation_id = $variation_id;
	}
	#only print genotypes without parents genotyped
	if (!exists $siblings->{$population_id . '-' . $individual_id}){ #necessary to use the population_id
	    #if it is a new variation, write to the file (if necessary) and empty the hash
	    if ($previous_variation_id ne $variation_id){
		foreach my $population (keys %alleles_variation){
		    #if the variation has 2 alleles, print all the genotypes to the file
		    if (keys %{$alleles_variation{$population}} == 2){		
			&convert_genotype($alleles_variation{$population},$genotype_information{$population});
			foreach my $individual_id (keys %{$genotype_information{$population}}){
			    &print_individual_file($buffer,$population, 
						   $previous_variation_id, $individual_id,
						   \%genotype_information,$dbname,\%genotypes_file,\%regions);
			}
		    }
		}
		$previous_variation_id = $variation_id;
		%alleles_variation = (); #new variation, flush the hash
		%genotype_information = (); #new variation, flush the hash
	    }
	    #we store the genotype information for the variation
	    if ($allele_1 ne 'N' and $allele_2 ne 'N'){
	      $genotype_information{$population_id}{$individual_id}{variation_feature_id} = $variation_feature_id;
	      $genotype_information{$population_id}{$individual_id}{seq_region_start} = $seq_region_start;
	      $genotype_information{$population_id}{$individual_id}{allele_1} = $allele_1;
	      $genotype_information{$population_id}{$individual_id}{allele_2} = $allele_2;
	      $genotype_information{$population_id}{$individual_id}{seq_region_end} = $seq_region_end;
	      $genotype_information{$population_id}{$individual_id}{seq_region_id} = $seq_region_id;
	      
	      #and the alleles
	      $alleles_variation{$population_id}{$allele_1}++;
	      $alleles_variation{$population_id}{$allele_2}++;
	    
	      $populations{$population_id}++;
	    }
	}
    }
    $sth->finish();
    #we have to print the last variation
    foreach my $population (keys %alleles_variation){
	#if the variation has 2 alleles, print all the genotypes to the file
	if (keys %{$alleles_variation{$population}} == 2){		
	    &convert_genotype($alleles_variation{$population},$genotype_information{$population});
	    foreach my $individual_id (keys %{$genotype_information{$population}}){
		&print_individual_file($buffer,$population, 
				       $previous_variation_id, $individual_id,
				       \%genotype_information,$dbname,\%genotypes_file,\%regions);
	    }
	}
    }

    print_buffered( $buffer );
    print "Time starting to submit jobs to queues: ",scalar(localtime(time)),"\n";

    $bsub_queue_name ||= 'normal';
    #let's run a job array
    foreach my $file (keys %genotypes_file){
      $call = "bsub -q $bsub_queue_name -J '$dbname.pairwise_ld' -m 'bc_hosts' ./ld_wrapper.sh $file $file\_out";	    
      system($call);
    }
    $call = "bsub -q $bsub_queue_name -w 'done($dbname.pairwise_ld)' -m 'ecs4_hosts' -o $TMP_DIR/output_ld_populations_import.txt $PERLBIN parallel_ld_populations.pl -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -tmpdir $TMP_DIR -tmpfile $TMP_FILE ";
    $call .= "-cpass $cpass " if ($cpass);
    $call .= "-cport $cport " if ($cport);
    $call .= "-vpass $vpass " if ($vpass);
#    system($call); #send the last job, that will wait for the rest to finish and upload everything to the table 

}


#prints to the buffer the different individuals
sub print_individual_file{
    my $buffer = shift;
    my $population = shift;
    my $previous_variation_id = shift;
    my $individual_id = shift;
    my $genotype_information = shift;
    my $dbname = shift;
    my $genotypes_file = shift;
    my $regions = shift;


    $regions->{"$TMP_DIR/$dbname.pairwise_ld_$population\.txt"}->{$genotype_information->{$population}->{$individual_id}->{seq_region_id}}++; #add one more genotype in the region for the population

    print_buffered( $buffer, "$TMP_DIR/$dbname.pairwise_ld_$population\.txt", 
		    join("\t",$previous_variation_id,$genotype_information->{$population}->{$individual_id}->{seq_region_start},
			 $individual_id,
			 $genotype_information->{$population}->{$individual_id}->{genotype},
			 $genotype_information->{$population}->{$individual_id}->{variation_feature_id},
			 $genotype_information->{$population}->{$individual_id}->{seq_region_id},
			 $genotype_information->{$population}->{$individual_id}->{seq_region_end},$population)."\n" );
    $genotypes_file->{"$TMP_DIR/$dbname.pairwise_ld_$population\.txt"}++;

}

sub print_buffered {
    my $buffer = shift;
    my $filename = shift;
    my $text = shift;

    local *FH;

    if( ! $filename ) {
	# flush the buffer
	foreach my $file (keys %{$buffer}){
	    open( FH, ">>$file" ) or die "Could not print to file $! \n";
	    print FH $buffer->{ $file };
	    close FH;
	}
	%{$buffer}=();

    } else {
	$buffer->{ $filename } .= $text;
	if( length( $buffer->{ $filename } ) > 10_000 ) {
	    open( FH, ">>$filename" ) or die;
	    print FH $buffer->{ $filename };
	    close FH;
	    $buffer->{ $filename } = '';
	}
    }
}


#
# Converts the genotype into the required format for the calculation of the pairwise_ld value: AA, Aa or aa
# From the Allele table, will select the alleles and compare to the alleles in the genotype
#
sub convert_genotype{
    my $alleles_variation = shift; #reference to the hash containing the alleles for the variation present in the genotypes
    my $genotype_information = shift; #reference to a hash containing the values to be written to the file
    my @alleles_ordered; #the array will contain the alleles ordered by apparitions in the genotypes (only 2 values possible)
    
    @alleles_ordered = sort({$alleles_variation->{$b} <=> $alleles_variation->{$a}} keys %{$alleles_variation});

    #let's convert the allele_1 allele_2 to a genotype in the AA, Aa or aa format, where A corresponds to the major allele and a to the minor
    foreach my $individual_id (keys %{$genotype_information}){
	#if both alleles are different, this is the Aa genotype
	if ($genotype_information->{$individual_id}{allele_1} ne $genotype_information->{$individual_id}{allele_2}){
	    $genotype_information->{$individual_id}{genotype} = 'Aa';
	}
	#when they are the same, must find out which is the major
	else{	    
	    if ($alleles_ordered[0] eq $genotype_information->{$individual_id}{allele_1}){
		#it is the major allele
		$genotype_information->{$individual_id}{genotype} = 'AA';
	    }
	    else{
		$genotype_information->{$individual_id}{genotype} = 'aa';
	    }
	    
	}
    }
}

#given a file with more than MAX_GENOTYPES, splits it in REGIONS different files
sub split_file{
    my $file = shift;
    my $regions = shift; #hash with all the regions in the population and the number of genotypes in each region
    my $genotypes = shift;
    my $dbname = shift;
    my @regions_ordered = sort {$a<=>$b} keys %{$regions->{$file}};
    my $sub_genotypes = int($genotypes/REGIONS()); #find minimum number of genotypes in each file
    #create the groups of regions for each file: the array will contain the number of genotypes (lines) that the group must contain
    my @groups;
    my $lines = 0;
    my $index = 1; #position in the group. The first position will contain n lines, the second $i*$n,...
    foreach my $region (@regions_ordered){
	$lines += $regions->{$file}->{$region};
	if ($lines > $sub_genotypes * $index){
	    push @groups,$lines;	    
	    $index++;
	}
    }
    open INFILE, "< $file" or die "Could not open file $file: $!";
    my $chunk = 1;
    until(eof INFILE) {
	open OUTFILE, "> $file\_chunk.$chunk"
	    or die "Could not open file $TMP_DIR/$dbname.pairwise_ld_chunk_$$\_$chunk\.txt: $!";	
	while(<INFILE>) {
	    print OUTFILE;
	    if ($chunk != REGIONS()){
		last unless $. % $groups[$chunk-1]; #last chunk will contain all remaining lines in the file
	    }
	}
	++$chunk;
	close OUTFILE;
    }
    close INFILE;
    unlink $file; #remove the original file after splitting it
}

#for a given population, gets all individuals that are children (have father or mother)
sub get_siblings{
    my $dbVariation = shift;
    my $population_id = shift;
    my $siblings = shift;

    my $sth_individual = $dbVariation->prepare(qq{SELECT i.sample_id
							     FROM individual i, individual_population ip
							     WHERE ip.individual_sample_id = i.sample_id
							     AND ip.population_sample_id = ? 
							     AND i.father_individual_sample_id IS NOT NULL
							     AND i.mother_individual_sample_id IS NOT NULL
							 });
    my ($individual_id);
    $sth_individual->execute($population_id);
    $sth_individual->bind_columns(\$individual_id);
    while ($sth_individual->fetch){
	$siblings->{$population_id.'-'.$individual_id}++; #necessary to have in the key the population, since some individuals are shared between
	                                                   #populations
    }
    return $siblings;
}


#method to crete a tmp table to store failed variations
sub create_failed_variation_table{
    my $dbVar = shift;

    $dbVar->do(qq{CREATE TABLE IF NOT EXISTS failed_variation(
			  variation_id int(10) unsigned not null,
			  failed_description_id int(10) unsigned not null,

			   PRIMARY KEY(variation_id))
		  }
	       );
}

#----------------------------------------------------------------------
# A wrapper for submitting PBS qsub jobs;
# First the READER (this process) opens a pipe (-|) on the WRITER.
# The WRITER then opens a pipe (|-) on the RUNNER.
# The RUNNER then execs qsub with command line options,
# the WRITER writes the script to the RUNNER, 
# and any output is collected by the READER.
# Returns the qsub job ID.
sub submit_qsub{
  my $script = shift || die( "Need a script to submit to qsub!" );
  my @qsub_args = @{ shift || [] };

  local *QSUB;
  local *QSUB_READER;
  
  my $jobid;
  
  my $writer_pid;
  if( $writer_pid = open( QSUB_READER, '-|' ) ){
    # READER - Reads stdout from RUNNER process
    while( <QSUB_READER> ){
      if( /^(\d+)/ ){
        $jobid = $1;
        print( "Job ID $1 submitted to qsub\n" );
      }
    }
    close( QSUB_READER );
  }
  else{
    unless( defined($writer_pid) ){ die("Could not fork : $!" ) }
    my $runner_pid;
    if( $runner_pid = open(QSUB, '|-') ){
      # WRITER - Writes command to RUNNER process
      print QSUB $script;
      close QSUB;
      unless( $? == 0 ){ die( "qsub failed; non-zero status" ) }
      exit(0);
    }
    else{
      # RUNNER - Runs the command; STDIN,STDOUT attached to WRITER,READER.
      unless( defined($runner_pid) ){ die("Could not fork : $!" ) }
      #warn join( " ", 'qsub', @qsub_args );
      exec( 'qsub', @qsub_args );
      die("Could not exec qsub : $!");
    }
  }
  return $jobid;
}

#----------------------------------------------------------------------

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl parallel_post_process.pl <options>

options:
    -species <string>       species in ensembl.registry
    -limit <number>         limit the number of rows for testing
    -perlbin <path>         path to perl interpreter (def/usr/bin/env perl)
    -scheduler <string>     job submission system to use (def LSF)
    -tmpdir <dir>           temp directory to use (with lots of space!)
    -tmpfile <filename>     name of temp file to use
    -num_processes <number> number of processes that are running (default = disabled)
    -variation_feature  fill in the Variation_feature table (default = disabled)
    -flanking_sequence  fill in the flanking sequence tables (default = disabled)
    -transcript_variation  fill in the Transcript_variation table (default = disabled)
    -ld_populations  fill in the Pairwise_ld table (default = disabled)
EOF

  die("\n$msg\n\n");
}

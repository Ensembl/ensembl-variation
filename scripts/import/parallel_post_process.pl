use strict;
use warnings;

use Getopt::Long;
use ImportUtils qw(debug);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use DBH;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($TMP_DIR, $TMP_FILE, $LIMIT);


my ($vhost, $vport, $vdbname, $vuser, $vpass,
    $chost, $cport, $cdbname, $cuser, $cpass,
    $limit, $num_processes, $top_level,
    $variation_feature, $flanking_sequence, $variation_group_feature,
    $transcript_variation, $ld_populations);

$variation_feature = $flanking_sequence = $variation_group_feature = $transcript_variation = $ld_populations = '';

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
	   'limit=i'   => \$limit,
	   'num_processes=i' => \$num_processes,
	   'variation_feature' => \$variation_feature,
	   'flanking_sequence' => \$flanking_sequence,
	   'variation_group_feature' => \$variation_group_feature,
	   'transcript_variation' => \$transcript_variation,
	   'ld_populations' => \$ld_populations);

#added default options
$chost    ||= 'ecs2';
$cuser    ||= 'ensro';
$cport    ||= 3365;

$vhost    ||='ecs2';
$vport    ||= 3361;
$vuser    ||= 'ensadmin';

$num_processes ||= 1;

$LIMIT = ($limit) ? " $limit " : ''; #will refer to position in a slice

usage('-vdbname argument is required') if(!$vdbname);
usage('-cdbname argument is required') if(!$cdbname);

usage('-num_processes must at least be 1') if ($num_processes == 0);

my $dbVar = DBH->connect
    ("DBI:mysql:host=$vhost;dbname=$vdbname;port=$vport",$vuser, $vpass );
die("Could not connect to variation database: $!") if(!$dbVar);

 my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host   => $chost,
     -user   => $cuser,
     -pass   => $cpass,
     -port   => $cport,
     -dbname => $cdbname);

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

##Apart from human and mouse, we directly import top_level coordinates from dbSNP
if ($cdbname !~ /homo|mus/i) {
  $top_level=1;
}

parallel_variation_feature($dbVar, $top_level) if ($variation_feature);
parallel_flanking_sequence($dbVar) if ($flanking_sequence);
parallel_variation_group_feature($dbVar) if ($variation_group_feature);
parallel_transcript_variation($dbVar) if ($transcript_variation);
parallel_ld_populations($dbVar) if ($ld_populations);


#will take the number of processes, and divide the total number of entries in the variation_feature table by the number of processes
sub parallel_variation_feature{
    my $dbVar = shift;
    my $top_level = shift;

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

    #create a temporary table to store the map_weight, that will be deleted by the last process
    $dbVar->do(qq{CREATE TABLE tmp_map_weight
                SELECT variation_id, count(*) as count
                FROM   variation_feature
                GROUP BY variation_id}
	       );
    $dbVar->do(qq{ALTER TABLE tmp_map_weight 
		      ADD INDEX variation_idx(variation_id,count)});

    my $sub_variation = int(($max_variation - $min_variation)/ $num_processes);
    #the limit will be (AND variation_feature_id > min and variation_feature_id < max)
    for (my $i = 0; $i < $num_processes ; $i++){
	$limit = "AND variation_feature_id <= " . (($i+1) * $sub_variation + $min_variation-1) . " AND variation_feature_id >= " . ($i*$sub_variation + $min_variation) if ($i+1 < $num_processes);
	$limit =  "AND variation_feature_id <= " .  $max_variation . " AND variation_feature_id >= " . ($i*$sub_variation + $min_variation) if ($i + 1 == $num_processes); #the last one takes the left rows
	$call = "bsub -J variation_job_$i -o $TMP_DIR/output_variation_feature_$i\_$$.txt /usr/local/ensembl/bin/perl parallel_variation_feature.pl -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -limit '$limit' -tmpdir $TMP_DIR -tmpfile $TMP_FILE -num_processes $num_processes -status_file $variation_status_file ";
	$call .= "-cpass $cpass " if ($cpass);
	$call .= "-cport $cport " if ($cport);
	$call .= "-vpass $vpass " if ($vpass);
	$call .= "-toplevel $top_level " if ($top_level);
	system($call);      
    }
    $call = "bsub -K -w 'done(variation_job*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
    system($call);
}

#will take the number of processes, and divide the total number of entries in the flanking_sequence table by the number of processes
#has to wait until the variation_feature table has been filled
sub parallel_flanking_sequence{
  my $dbVar = shift;
  
  my $total_process = 0;
  my $sequences; #number of entries in the flanking_sequence table
  my $call;
  my $flanking_status_file = "flanking_status_file_$$\.log";
  #first, create the log file for the variation_feature
  open STATUS, ">$TMP_DIR/$flanking_status_file"
    or throw("Could not open tmp file: $TMP_DIR/$flanking_status_file\n"); 
  close STATUS;
  #then, calculate the rows for each subprocess
  my $sth = $dbVar->prepare(qq{SELECT count(*)
			       FROM flanking_sequence
			      });
  $sth->execute();
  ($sequences) = $sth->fetchrow_array();
  $sth->finish();
  my $sub_sequences = int($sequences / $num_processes);
  
  for (my $i = 0; $i < $num_processes ; $i++){
   ###num_processes has to be >1, otherwise $i=0, cause limit 0,0
					      
   $limit = $i*$sub_sequences . "," . $sub_sequences  if ($i+1 < $num_processes or $num_processes ==1);
   $limit = $i*$sub_sequences . "," . $sub_sequences*$i if ($i+1 == $num_processes and $num_processes !=1); #the last one will select the left rows
   $call = "bsub -o $TMP_DIR/output_flanking_$i\_$$.txt /usr/local/ensembl/bin/perl parallel_flanking_sequence.pl -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -limit $limit -tmpdir $TMP_DIR -tmpfile $TMP_FILE -num_processes $num_processes -status_file $flanking_status_file ";
   $call .= "-cpass $cpass " if ($cpass);
   $call .= "-cport $cport " if ($cport);
   $call .= "-vpass $vpass " if ($vpass);
   system($call);
  }
  
}

#when the variation_feature table has been filled up, run the variation_group_feature. Not necessary to parallelize as fas as I know....
sub parallel_variation_group_feature{
    my $dbVar = shift;

    my $total_process = 0;
    my $call = "bsub -o $TMP_DIR/output_group_feature_$$\.txt /usr/local/ensembl/bin/perl parallel_variation_group_feature.pl -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -tmpdir $TMP_DIR -tmpfile $TMP_FILE ";
    $call .= "-cpass $cpass " if ($cpass);
    $call .= "-cport $cport " if ($cport);
    $call .= "-vpass $vpass " if ($vpass);
    $call .= "-limit $limit" if ($limit);
    system($call);
}

#will fill in the transcript variation table. Has to wait until the variation feature table has been filled up. Then, divide the number of entries
# by the number of processes to make the subprocesses
sub parallel_transcript_variation{
    my $dbVar = shift;

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
    
    #I must add up the length of all the slices to find the limit for each 
    my $sub_slice = int($length_slices / $num_processes);
    my $slice_max; #the number of slices in the chunk
    my $slice_min = 0; #first slice in the chunk
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
	$limit = $slice_min . "," . (scalar(@slices_ordered)-$slice_min-1) 
	  if ($i+1 == $num_processes and $num_processes != 1); #the last slice, get the left slices

	$call = "bsub -o $TMP_DIR/output_transcript_$i\_$$.txt /usr/local/ensembl/bin/perl parallel_transcript_variation.pl -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -limit $limit -tmpdir $TMP_DIR -tmpfile $TMP_FILE -num_processes $num_processes -status_file $transcript_status_file ";
	$call .= "-cpass $cpass " if ($cpass);
	$call .= "-cport $cport " if ($cport);
	$call .= "-vpass $vpass " if ($vpass);
	$slice_min = $slice_max;
	system($call);
    }
}

#will have to wait until the variation_feature has finished. Then, divide the total number of populations by the processes, and run it
sub parallel_ld_populations{
    my $dbVar = shift;

    my $total_process = 0;

    my $populations; #number of populations with genotypes
    my $call;
    my $ld_populations_status_file = "ld_populations_status_file_$$\.log";
#first, create the log file for the variation_feature
    open STATUS, ">$TMP_DIR/$ld_populations_status_file"
	or throw("Could not open tmp file: $TMP_DIR/$ld_populations_status_file\n"); 
    close STATUS;
    #then, calculate the populations
    my $sth = $dbVar->prepare(qq{SELECT count(distinct i.population_id)
				     FROM individual i, individual_genotype ig 
				     WHERE ig.individual_id = i.individual_id				     
				 });
    $sth->execute();
    ($populations) = $sth->fetchrow_array();
    $sth->finish();
    my $sub_populations = int($populations / $num_processes);
    for (my $i = 0; $i < $num_processes ; $i++){

	$limit = $i*$sub_populations . "," . $sub_populations if ($i+1 < $num_processes or $num_processes==1);
	$limit = $i*$sub_populations . "," . $sub_populations*$i 
	  if ($i+1 == $num_processes and $num_processes !=1); #the last one will select the left rows

	$call = "bsub -o $TMP_DIR/output_ld_populations_$i\_$$.txt /usr/local/ensembl/bin/perl parallel_ld_populations.pl -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -limit $limit -tmpdir $TMP_DIR -tmpfile $TMP_FILE -num_processes $num_processes -status_file $ld_populations_status_file ";
	$call .= "-cpass $cpass " if ($cpass);
	$call .= "-cport $cport " if ($cport);
	$call .= "-vpass $vpass " if ($vpass);
	system($call);
    }

}
sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl parallel_post_process.pl <options>

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
    -num_processes <number> number of processes that are running (default = disabled)
    -variation_feature  fill in the Variation_feature table (default = disabled)
    -flanking_sequence  fill in the flanking sequence tables (default = disabled)
    -variation_group_feature fill in the Variation_group_feature table (default = disabled)
    -transcript_variation  fill in the Transcript_variation table (default = disabled)
    -ld_populations  fill in the Pairwise_ld table (default = disabled)
EOF

  die("\n$msg\n\n");
}

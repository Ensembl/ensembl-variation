use strict;
use warnings;

use Getopt::Long;
use Fcntl ':flock';
use DBI;
use DBH;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);

use ImportUtils qw(debug load);

use constant MAX_GENOTYPES => 1_000_000; #max number of genotypes per file. When more, split the file into regions
use constant REGIONS => 5; #number of chunks the file will be splited when more than MAX_GENOTYPES

my ($TMP_DIR, $TMP_FILE, $LIMIT,$status_file);

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
             'limit=s'   => \$limit,
	     'num_processes=i' => \$num_processes,
	     'status_file=s' => \$status_file);

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

  usage('-num_processes must at least be 1') if ($num_processes == 0);
  usage('-status_file argument is required') if (!$status_file);

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

  ld_populations($dbCore,$dbVar);
  open STATUS, ">>$TMP_DIR/$status_file"
    or throw("Could not open tmp file: $TMP_DIR/$status_file\n"); 
  flock(STATUS,LOCK_EX);
  seek(STATUS, 0, 2); #move to the end of the file
  print STATUS "process finished\n";
  flock(STATUS,LOCK_UN);
  close STATUS;
  #check if it is the last process
  my $processes = `cat $TMP_DIR/$status_file | wc -l`;
  if ($processes == $num_processes){
      #if is the last process, delete the variation table and upload with the new coordinate information from the different tables
      last_process($dbCore,$dbVar);
  }
}

#for a list of populations, creates the pairwise_ld information
sub ld_populations{
    my $dbCore = shift;
    my $dbVar = shift;
    my $population_id;
    my $sth = $dbVar->prepare(qq{SELECT distinct(i.population_id)
				     FROM individual i, individual_genotype ig 
				     WHERE ig.individual_id = i.individual_id
				     ORDER BY i.population_id
				     $LIMIT
    });
    $sth->execute();
    $sth->bind_columns(\$population_id);
    debug("calculating ld distance");
    #for all the populations, we calculate the pairwise_ld informaiton
    while ($sth->fetch()){
	&pairwise_ld($dbCore,$dbVar,$population_id);
    }
    $sth->finish;
}
#
# populates the table pairwise_ld with the information for the genotypes for a certain population
#
sub pairwise_ld{
    my $dbCore = shift;
    my $dbVar = shift;
    my $population_id = shift;

    my %seq_region; #hash containing the mapping between seq_region_id->name region
    my %alleles_variation = (); #will contain a record of the alleles in the variation. A will be the major, and a the minor. When more than 2 alleles
    my %genotype_information; #will contain all the genotype information to write in the file, if necessary
    #, the genotypes for that variation will be discarded
    my %regions; #will contain all the regions in the population and the number of genotypes in each one
    my $previous_variation_id = ''; #to know if it is a new variation and we can get the new alleles
    my $genotype; #genotype in the AA, Aa or aa format to write to the file
#    debug("Loading pairwise_ld table for population $population_id\n");

    #necessary the order to know when we change variation. Not get genotypes with a NULL variation or map_weight > 1
    my $sth = $dbVar->prepare
	(qq{SELECT ig.variation_id, vf.variation_feature_id, vf.seq_region_id, vf.seq_region_start, ig.individual_id, ig.allele_1, ig.allele_2, vf.seq_region_end
		FROM   variation_feature vf, individual_genotype ig, individual i
		WHERE  ig.variation_id = vf.variation_id
		AND i.individual_id = ig.individual_id
		AND i.population_id = ?
		AND ig.allele_2 IS NOT NULL
		AND vf.map_weight = 1
		ORDER BY vf.seq_region_id,vf.seq_region_start
		});
    $sth->execute($population_id);
    my $dbname = $dbVar->dbname(); #get the name of the database to create the file
    open ( FH, ">$TMP_DIR/$dbname.pairwise_ld_$$\.txt");

    #retrieve all the genotypes and format them into the calculation_genotype required format:
    # snp_id chr position individual_id genotype
    my ($variation_id, $variation_feature_id, $seq_region_id, $seq_region_start, $individual_id, $allele_1,$allele_2,$seq_region_end);
    my $seq_region_name;
    $sth->bind_columns(\$variation_id, \$variation_feature_id, \$seq_region_id, \$seq_region_start, \$individual_id, \$allele_1,\$allele_2,\$seq_region_end);
    my $count_not_discarded = 0;
    my $count_discarded = 0;
    while ($sth->fetch()){
	if ($previous_variation_id eq ''){
	    $previous_variation_id = $variation_id;
	}
	#if it is a new variation, write to the file (if necessary) and empty the hash
	if ($previous_variation_id ne $variation_id){
	    #if the variation has 2 alleles, print all the genotypes to the file
	    if (keys %alleles_variation == 2){		
		$count_not_discarded++;
		&convert_genotype(\%alleles_variation,\%genotype_information);
		foreach my $individual_id (keys %genotype_information){
		    $regions{$genotype_information{$individual_id}{seq_region_id}}++; #add one more genotype in the region
		    print FH join("\t",$previous_variation_id,$genotype_information{$individual_id}{seq_region_start},$individual_id,$genotype_information{$individual_id}{genotype},$genotype_information{$individual_id}{variation_feature_id},$genotype_information{$individual_id}{seq_region_id},$genotype_information{$individual_id}{seq_region_end},$population_id),"\n";
		}
	    }
	    else{
		$count_discarded++;
	    }
	    $previous_variation_id = $variation_id;
	    %alleles_variation = (); #new variation, flush the hash
	    %genotype_information = (); #new variation, flush the hash
	}
	#we store the genotype information for the variation
	$genotype_information{$individual_id}{variation_feature_id} = $variation_feature_id;
	$genotype_information{$individual_id}{seq_region_start} = $seq_region_start;
	$genotype_information{$individual_id}{allele_1} = $allele_1;
	$genotype_information{$individual_id}{allele_2} = $allele_2;
	$genotype_information{$individual_id}{seq_region_end} = $seq_region_end;
	$genotype_information{$individual_id}{seq_region_id} = $seq_region_id;

	#and the alleles
	$alleles_variation{$allele_1}++;
	$alleles_variation{$allele_2}++;
    }
#    print "Number of variations discarded for population $population_id are: $count_discarded and used $count_not_discarded\n\n";
    $sth->finish();
    close FH;
    #to make a better use of the farm, break the population into smaller segments when the population contains more than 2_000_000
    #check if it is the last process
    my $lines = `cat $TMP_DIR/$dbname.pairwise_ld_$$\.txt | wc -l`;
    #when there are more genotypes than MAX, split the file into regions to better use the farm and concatenate them at the end
    if ($lines > MAX_GENOTYPES()){
	&split_file("$TMP_DIR/$dbname.pairwise_ld_$$\.txt",\%regions,$lines,$dbname);
	foreach my $i (1..REGIONS){
	    system("bsub -J $dbname\_ld_job_$$\_$i -o $TMP_DIR/output_ld_populations_region_$i\_$$.txt /usr/local/ensembl/bin/perl calc_genotypes.pl $TMP_DIR/$dbname.pairwise_ld_chunk_$i\_$$\.txt $TMP_DIR/$dbname.pairwise_ld_chunk_out_$i\_$$\.txt");
	}
	system("bsub -K -w 'done($dbname\_ld_job_$$\_*)' -J waiting_process sleep 1"); #wait for the calculation to finish
	system("cat $TMP_DIR/$dbname.pairwise_ld_chunk_out_?_$$\.txt > $TMP_DIR/$dbname.pairwise_ld_out_$$\.txt"); #concat all the pieces into the final file	
    }
    else{
	#and, finally, create the table with the information, calling the calc_genotypes.pl script
	system("perl /nfs/acari/dr2/projects/src/ensembl/ensembl-variation/scripts/import/calc_genotypes.pl $TMP_DIR/$dbname.pairwise_ld_$$\.txt $TMP_DIR/$dbname.pairwise_ld_out_$$\.txt");
    }
    #delete the file with the original file
    unlink("$TMP_DIR/$dbname.pairwise_ld_$$\.txt");
    return;
}

#given a file with more than MAX_GENOTYPES, splits it in REGIONS different files
sub split_file{
    my $file = shift;
    my $regions = shift; #hash with all the regions in the population and the number of genotypes in each region
    my $genotypes = shift;
    my $dbname = shift;
    my @regions_ordered = sort {$a<=>$b} keys %{$regions};
    my $sub_genotypes = int($genotypes/REGIONS()); #find minimum number of genotypes in each file
    #create the groups of regions for each file: the array will contain the number of genotypes (lines) that the group must contain
    my @groups;
    my $lines = 0;
    my $index = 1; #position in the group. The first position will contain n lines, the second $i*$n,...
    foreach my $region (@regions_ordered){
	$lines += $regions->{$region};
	if ($lines > $sub_genotypes * $index){
	    push @groups,$lines;	    
	    $index++;
	}
    }
    open INFILE, "< $file" or die "Could not open file $file: $!";
    my $chunk = 1;
    until(eof INFILE) {
	open OUTFILE, "> $TMP_DIR/$dbname.pairwise_ld_chunk_$chunk\_$$\.txt"
	    or die "Could not open file $TMP_DIR/$dbname.pairwise_ld_chunk_$chunk\_$$\.txt: $!";	
	while(<INFILE>) {
	    print OUTFILE;
	    last unless $. % $groups[$chunk-1]; #will return true when the line
	}
	++$chunk;
	close OUTFILE;
    }
    close INFILE;
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

#will know if it is the last process running counting the lines of the status_file.If so, load all data
sub last_process{
    my $dbCore = shift;
    my $dbVar = shift;



    debug("Importing pairwise data");
   #group all the fragments in 1 file
    my $dbname = $dbVar->dbname(); #get the name of the database to create the file
    my $call = "cat $TMP_DIR/$dbname.pairwise_ld_out*.txt > $TMP_DIR/$TMP_FILE";
    system($call);
    unlink(<$TMP_DIR/$dbname.pairwise_ld_out*>);    
    #remove chunk files created
    unlink(<$TMP_DIR/$dbname.pairwise_ld_chunk_*>);

    #and import the data in the database
    load($dbVar, qw(pairwise_ld variation_feature_id_1 variation_feature_id_2 population_id seq_region_id seq_region_start seq_region_end snp_distance_count r2 d_prime sample_count));


    update_meta_coord($dbCore, $dbVar, 'pairwise_ld');
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

#  $sth->execute($table_name, $cs->dbID());

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

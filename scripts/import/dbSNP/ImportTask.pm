use strict;
use warnings;

#generic object for the dbSNP data. Contains the general methods to dump the data into the new Variation database. Any change in the methods
# will need to overload the correspondent method in the subclass for the species

package dbSNP::ImportTask;

use ImportUtils qw(dumpSQL debug create_and_load load loadfile get_create_statement);
use Progress;
use DBI qw(:sql_types);
use Fcntl qw( LOCK_SH LOCK_EX );
use Digest::MD5 qw ( md5_hex );
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp );
use Bio::EnsEMBL::Variation::Utils::dbSNP qw(get_alleles_from_pattern);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(revcomp_tandem);


our %QUICK_COMP = ( "A" => "T",
                    "T" => "A",
                    "C" => "G",
                    "G" => "C"
    ); 

sub new {
  my $class = shift;
  my $dbm = shift;
  my $task = shift;
  my $file_prefix = shift;
  my $unique_suffix = shift;
  my $log = shift;
  
  my $ref;
  $ref->{'db_manager'} = $dbm;
  $ref->{'task'} = $task;
  $ref->{'file_prefix'} = $file_prefix;
  $ref->{'unique_suffix'} = $unique_suffix;
  $ref->{'log'} = $log;
  
  return bless($ref,$class);
}


######## export allele data (with frequencies where possible) from  local dbSNP database 
sub allele_table {

    
  my $self = shift;
  my $loadfile = shift;
  my $task_start = shift;
  my $task_end = shift;

  my $logh = $self->{'log'};
 
  my $stmt;
  my $dbm = $self->{'db_manager'};
  my $shared_db = $dbm->dbSNP_shared();
  
  print  Progress::location() . "\tImporting for SubSNP ids from $task_start to $task_end\n";
  
 # Prepared statement to get the alleles
  $stmt = qq{
    SELECT
      a.allele,
      ar.allele
    FROM
      $shared_db..Allele a JOIN
      $shared_db..Allele ar ON (
      ar.allele_id = a.rev_allele_id
      )
    WHERE
      a.allele_id = ?
  };
  my $allele_sth = $dbm->dbSNP()->dbc->prepare($stmt)||die "ERROR preparing allele_sth: $DBI::errstr\n";
  


  #ÊPrepared statement for getting the SubSNPs with or without population frequency data
  # 2012/09/24 switch from UniVariAllele at dbSNP request
   my $ss_no_freq_stmt = qq{
     SELECT
        ss.subsnp_id,     
        b.pop_id,
        ov.pattern,
        sssl.substrand_reversed_flag
    FROM   SNPSubSNPLink sssl,
           SubSNP ss,
           Batch b,
           $shared_db..ObsVariation  ov
        WHERE ss.subsnp_id BETWEEN ? AND ?
        AND   ss.subsnp_id = sssl.subsnp_id
        AND   b.batch_id = ss.batch_id       
        AND   ov.var_id = ss.variation_id
    ORDER BY
      ss.subsnp_id ASC
  };

## Can't filter s with freq from those without easily due to circular schena

  my $ss_no_freq_sth = $dbm->dbSNP()->dbc->prepare($ss_no_freq_stmt) ||die "ERROR preparing ss_no_freq_sth: $DBI::errstr\n";
  
  # Prepared statement to get all SubSNPs that have population frequency data
  my  $ss_with_freq_stmt = qq{
    SELECT
      ss.subsnp_id,
      afbsp.pop_id,
      afbsp.allele_id ,
      sssl.substrand_reversed_flag,
      afbsp.freq,
      afbsp.cnt,
      pop.handle 
    FROM
      SNPSubSNPLink sssl,
      SubSNP ss,
      AlleleFreqBySsPop afbsp, 
      Population pop      
    WHERE ss.subsnp_id BETWEEN ? AND ?
    AND   ss.subsnp_id = sssl.subsnp_id  
    AND   afbsp.subsnp_id = ss.subsnp_id
    AND   afbsp.pop_id = pop.pop_id
    ORDER BY
      ss.subsnp_id ASC,   
      afbsp.pop_id ASC
  };
  my $ss_freq_sth = $dbm->dbSNP()->dbc->prepare($ss_with_freq_stmt)  ||die "ERROR preparing ss_freq_sth: $DBI::errstr\n";
  
  #ÊPrepared statement to get the corresponding variation_ids for a range of subsnps from variation_synonym
  my $vs_stmt = qq{
    SELECT
      vs.subsnp_id,
      vs.variation_id
    FROM
      variation_synonym vs
    WHERE
      vs.subsnp_id BETWEEN ? AND ?
  };
  my $vs_sth = $dbm->dbVar()->dbc->prepare($vs_stmt) ||die "ERROR preparing vs_sth: $DBI::errstr\n";

  my %handle;
  $handle{0} = '\N';
  my $handle_ext_sth  = $dbm->dbVar()->dbc->prepare(qq[select handle_id from submitter_handle where handle =?])||die "ERROR preparing handle_ins_sth: $DBI::errstr\n";
  my $handle_ins_sth  = $dbm->dbVar()->dbc->prepare(qq[insert into submitter_handle (handle) values (?)])||die "ERROR preparing handle_ins_sth: $DBI::errstr\n";


  # Prepared statement to get the sample_id from a pop_id
  my $sample_stmt = qq{
    SELECT
      s.sample_id
    FROM
      sample s
    WHERE
      s.pop_id = ?
    LIMIT 1
  };
  my $sample_sth = $dbm->dbVar()->dbc->prepare($sample_stmt)||die "ERROR preparing sample_sth: $DBI::errstr\n";
    #ÊHash to keep sample_ids in memory
  my %samples;
  # The pop_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value
  $samples{0} = '\N';


  #ÊHash to hold the alleles in memory
  my %alleles;
  # The allele_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value
  $alleles{0} = ['\N','\N'];
  

  
  # Hash to hold the output until it's time to write to disk
  my @output;
    
 
 
  # Fetch the subsnp_id to variation_id mapping and store as hashref
  $vs_sth->execute($task_start,$task_end)||die "ERROR getting the SubSNPs variation id mapping: $DBI::errstr\n";
  my $variation_ids = $vs_sth->fetchall_hashref(['subsnp_id']);
 
  
  ## write values to file for later importation:
  ## variation_id,  subsnp_id, samples_id, alleles[strand appropriate], frequency, count, frequency submitter handle
  
  my %done;  ## remove duplicates


  # Fetch the alleles with frequency data
  $ss_freq_sth->execute($task_start,$task_end )||die "ERROR getting the SubSNPs with population frequency data: $DBI::errstr\n";;

  while( my $line = $ss_freq_sth->fetchrow_arrayref()){

      ## line content: [ss.subsnp_id,afbsp.pop_id, allele_id,sssl.substrand_reversed_flag,afbsp.freq,afbsp.cnt, pop.handle ]
      unless ($variation_ids->{$line->[0]}{'variation_id'}){ warn "Variation id not found for $line->[0]\n";}

     ## look-ups from ensembl db
     unless( defined $line->[6] ){ $line->[6] = 0 ; }
     if (!exists($samples{$line->[1]})) {  $samples{$line->[1]} = get_sample($line->[1], $sample_sth);}     
     if (!exists($handle{$line->[6]}))  {  $handle{$line->[6]}  = get_handle($line->[6], $handle_ext_sth, $handle_ins_sth);}
     if (!exists($alleles{$line->[2]})) {  $alleles{$line->[2]} = get_allele($line->[2], $allele_sth);}   # returns alleles as  [for, rev]    


     unless (defined $alleles{$line->[2]}->[$line->[3]]){ die "No allele for ss $line->[0] strand $line->[3]\n";}


      ## save as done by subsnp id pop and allele (don't loose 2nd allele for non-poly)
      $done{$line->[0]}{$line->[1]}{$alleles{$line->[2]}->[$line->[3]]} = 1;
      ##$done{$line->[0]}{0} = 1;  ## don't want to import records without population info if records with held


     my $row = join("\t",($variation_ids->{$line->[0]}{'variation_id'},$line->[0],$samples{$line->[1]},$alleles{$line->[2]}->[$line->[3]],$line->[4],$line->[5],$handle{$line->[6]} )) ;

     push(@output,$row);
  }



  ########## rest of alleles


  # Fetch all alleles and exclude those already written
  $ss_no_freq_sth->execute($task_start,$task_end )||die "ERROR getting the SubSNPs with or without population frequency data: $DBI::errstr\n";

  # line content: [ss.subsnp_id,  b.pop_id,uv.allele_id ,  sssl.substrand_reversed_flag,'\\N' AS frequency, \\N' AS count,0 AS handle]
  while( my $line = $ss_no_freq_sth->fetchrow_arrayref()){

      unless ($variation_ids->{$line->[0]}{'variation_id'}){ warn "Variation id not found for $line->[0]\n"; next}

      ## default pop_id to enter as null
      unless( defined $line->[1] ){ $line->[1] = 0; }        
      if (!exists($samples{$line->[1]})) {  $samples{$line->[1]} = get_sample($line->[1], $sample_sth);}

      my $sep_alleles = get_alleles_from_pattern($line->[2]);    ## Reported allele pattern [ eg A/T ]
      foreach my $sep_allele (@$sep_alleles){
	  if ($line->[3] ==1 ){
	      $line->[2] =~/^\(\w+\)/ ?	 $sep_allele = revcomp_tandem($sep_allele):	      
		  defined $QUICK_COMP{$sep_allele} ? $sep_allele = $QUICK_COMP{$sep_allele} : reverse_comp(\$sep_allele);
	  }
          ## don't add a line without frequency information if record with frequency already entered for this ss & pop & allele
          next if defined  $done{$line->[0]}{$line->[1]}{$sep_allele} ;

	  ## ens variation_id, subsnp_id, ens sample_id, allele, frequency, count, handle
	  my $row = join("\t",($variation_ids->{$line->[0]}{'variation_id'},$line->[0],$samples{$line->[1]},$sep_allele,'\N','\N','\N' )) ;
	  
	  push(@output,$row);
      }
  }
   print    Progress::location() . "Saved no-freq data\n";
 


  open(IMP,'>>',$loadfile) ||die "ERROR opening loadfile : $!\n";
  # Lock the file
  flock(IMP,LOCK_EX);
  foreach my $row (@output) {
    print IMP $row . "\n";
  }
  close(IMP);


  # We're done, return success
  return 1;
}

## extract ensembl database id for pre-imported sample
sub get_sample{

    my $pop_id     = shift;
    my $sample_sth = shift;

    $sample_sth->execute($pop_id)|| die "Failed to find sample id for $pop_id: $DBI::errstr\n";
    
    my $sample_id = $sample_sth->fetchall_arrayref();

    if(defined $sample_id->[0]->[0]){
        return $sample_id->[0]->[0];
    }
    else{
        return '\N';
    }
}

## convert dbSNP allele id to letter representation of base
## not done as join in main query due db holding 2 rows per id
sub get_allele{

    my $allele_id      = shift;
    my $allele_sth     = shift;

    $allele_sth->execute($allele_id)|| die "Failed to find allele id for $allele_id: $DBI::errstr\n";

    ## forward and reverse returned
    my  $alleles = $allele_sth->fetchall_arrayref();
    ## curently first row is A, C second (A)1 (C)1
    return $alleles->[0];

}

## get handle id from e! database or enter if new 
sub get_handle{

    my $handle         = shift;
    my $handle_ext_sth = shift;
    my $handle_ins_sth = shift;


    $handle_ext_sth->execute($handle)|| die "Failed to find handle id for $handle: $DBI::errstr\n";
    my  $id = $handle_ext_sth->fetchall_arrayref();

    unless(defined $id->[0]->[0]){ ## handle unknown => enter
       $handle_ins_sth->execute($handle)|| die "Failed to insert handle id for $handle: $DBI::errstr\n";

        $handle_ext_sth->execute($handle)|| die "Failed to find handle id for $handle: $DBI::errstr\n";
        $id = $handle_ext_sth->fetchall_arrayref();
    }

   return $id->[0]->[0];
  

}

sub calculate_gtype {
  my $self = shift;
  my $subind_table = shift;
  my $loadfile = shift;
  my $mlt_file = shift;
  my $start = shift;
  my $end = shift;
  my $mapping_file = shift;
  my $allele_file = shift;
  my $sample_file = shift;
  
  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  print $logh Progress::location() . "\tDumping from $subind_table to $loadfile, processing genotype data for submitted_ind_ids $start - $end\n";  
  
  my $dbm = $self->{'db_manager'};
  my $shared_db = $dbm->dbSNP_shared();
  my $dbVar = $dbm->dbVar();
  my $dbSNP = $dbm->dbSNP();
  
  #ÊA prepared statement for getting the genotype data from the dbSNP mirror. Note the chr_num, dbSNP is storing the gty allele for each chromosome copy? Anyway, we'll get duplicated rows if not specifying this (or using distinct)
  #ÊOrder the results by subsnp_id so that we can do the subsnp_id -> variation_id lookup more efficiently
  my $stmt = qq{
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
      ga.rev_flag,
      LEN(ug.gty_str) AS pattern_length,
      b.handle,
      si.submitted_ind_id
    FROM   
      $subind_table si JOIN 
      $shared_db..GtyAllele ga ON (
	ga.gty_id = si.gty_id
      ) JOIN
      $shared_db..UniGty ug ON (
	ug.unigty_id = ga.unigty_id
      ) JOIN
      SubmittedIndividual sind ON (
	sind.submitted_ind_id = si.submitted_ind_id
      ) JOIN Batch b ON ( b.batch_id = si.batch_id)
       
    WHERE
      si.submitted_ind_id BETWEEN ? AND ? AND
      ga.chr_num = 1
    ORDER BY
      si.subsnp_id ASC
  };
  my $subind_sth = $dbSNP->dbc()->prepare($stmt) ||die;
  
  #ÊPrepared statement to get the corresponding variation_ids for a range of subsnps from variation_synonym
  $stmt = qq{
    SELECT
      vs.subsnp_id,
      vs.variation_id,
      vs.substrand_reversed_flag
    FROM
      variation_synonym vs
    WHERE
      vs.subsnp_id = ?
  };
  my $vs_sth = $dbVar->dbc()->prepare($stmt);
  
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
  my $sample_sth = $dbVar->dbc()->prepare($stmt);
  
  # Prepared statement to get the alleles. We get each one of these when needed.
  $stmt = qq{
    SELECT
      a.allele,
      ar.allele
    FROM
      $shared_db..Allele a JOIN
      $shared_db..Allele ar ON (
	ar.allele_id = a.rev_allele_id
      )
    WHERE
      a.allele_id = ?
  };
  my $allele_sth = $dbSNP->dbc()->prepare($stmt);
  

 ## hack to add individuals missing from main import
  my %all_individuals;
  my $individual_type_id;
  if ($dbm->dbCore()->species() =~ /homo|pan|anoph/i) {  $individual_type_id = 3; }
  elsif ($dbm->dbCore()->species() =~ /mus/i) {          $individual_type_id = 1; }
  else {                                                 $individual_type_id = 2; }

  my $miss_ind_ext_sth = $dbSNP->dbc()->prepare(qq[select si.loc_ind_id_upp,
                                                          pop.loc_pop_id_upp
                                                    from SubmittedIndividual si
                                                    LEFT OUTER JOIN Population pop on (pop.pop_id = si.pop_id)
                                                    where si.submitted_ind_id =?]);


   my $ind_ins_sth = $dbVar->dbc()->prepare(qq [INSERT INTO individual 
                                               (sample_id, gender, individual_type_id)
		                               values ( ?,'Unknown', $individual_type_id)
		                                ]);

    my $sam_ins_sth = $dbVar->dbc()->prepare(qq [ INSERT INTO sample
                                                (name, description)
             	   	                         values (?,?)
             		                        ]);
  
    my $sam_ext_sth = $dbVar->dbc()->prepare(qq [ select sample_id
                                                 from sample 
             	   	                         where name =? and description=?
             		                        ]);
  
    my $samind_ext_sth = $dbVar->dbc()->prepare(qq [ select sample.sample_id
                                                 from sample, individual 
             	   	                         where sample.name =? and description=?
          	   	                         and sample.sample_id = individual.sample_id
             		                        ]);

    my $sampop_ext_sth = $dbVar->dbc()->prepare(qq [ select sample.sample_id
                                                 from sample, population 
             	   	                         where sample.name =? and description=?
          	   	                         and sample.sample_id = population.sample_id
             		                        ]);
  

    my $pop_ins_sth = $dbVar->dbc()->prepare(qq [ INSERT INTO population
                                                (sample_id)
             	   	                         values (?)
             		                        ]);
  
   
    my $indpop_ins_sth = $dbVar->dbc()->prepare(qq [ INSERT INTO individual_population
                                                (individual_sample_id, population_sample_id )
             	   	                         values (?,?)
             		                        ]);
  




  #ÊHash to hold the alleles in memory
  my %alleles;
  #ÊIf the allele file exist, we'll read alleles from it
  %alleles = %{read_alleles($allele_file)} if (defined($allele_file) && -e $allele_file);
  print $logh Progress::location();
  # The allele_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value.
  $alleles{0} = ['\N','\N'];
  # Keep track if we did any new lookups
  my $new_alleles = 0;
  
  #ÊHash to keep sample_ids in memory
  my %samples;
  #ÊIf the sample file exist, we'll read alleles from it
  %samples = %{read_samples($sample_file)} if (defined($sample_file) && -e $sample_file);
  print $logh Progress::location();
  # The individual_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value
  $samples{0} = '\N';
  # Keep track if we did any new lookups
  my $new_samples = 0;
  
  #ÊHash to keep subsnp_id to variation_id mappings in memory
  my %variation_ids;
  #ÊIf the mapping file exist, we'll read mappings from it
  %variation_ids = %{read_subsnp_mapping($mapping_file)} if (defined($mapping_file) && -e $mapping_file);
  print $logh Progress::location();
  # The subsnp_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value
  $variation_ids{0} = ['\N','\N'];
  # Keep track if we did any new lookups
  my $new_mappings = 0;
  
  #ÊKeep the genotypes in an array and write them in the end to avoid locking the file more than necessary. Separate signle and multiple bp genotypes. Use a hash having the md5sum of the row to determine whether we've already used this one 
  my @single_bp_gty;
  my @multiple_bp_gty;
  my %row_md5s;
  
  # First, get all SubSNPs
  $subind_sth->bind_param(1,$start,SQL_INTEGER);
  $subind_sth->bind_param(2,$end,SQL_INTEGER);
  
  # Fetch the import data as an arrayref
  $subind_sth->execute();
  print $logh Progress::location();
  
  #ÊNow, loop over the import data and print it to the tempfile so we can import the data. Replace the allele_id with the corresponding allele on-the-fly
  while (my $data = $subind_sth->fetchrow_arrayref()) {
    my $subsnp_id = $data->[0];
    my $ind_id = "MISS";
    if(defined $data->[1]){ $ind_id = $data->[1] ;}
    my $allele_1 = $data->[2];
    my $allele_2 = $data->[3];
    my $sub_strand = $data->[4];
    my $rev_alleles = $data->[5];
    my $pattern_length = $data->[6];
    my $handle         = $data->[7];
    my $ind_sub_id     = $data->[8];
    # If pattern length is less than 3, skip this genotype
    next if ($pattern_length < 3);
    # If ind_id is null (there is a case where this happens because one individual in SubmittedIndividual is missing from Individual (131)), skip this genotype.
    ##next if (!defined($ind_id)); e!69 - fix this

    # Look up the variation_id if necessary. This should be slow for the first chunk on each chromosome but fast for the rest..
    if (!exists($variation_ids{$subsnp_id})) {
      $vs_sth->bind_param(1,$subsnp_id,SQL_INTEGER);
      $vs_sth->execute();
      my ($vs_subsnp_id,$variation_id,$substrand_reversed_flag);
      $vs_sth->bind_columns(\$vs_subsnp_id,\$variation_id,\$substrand_reversed_flag);
      $vs_sth->fetch();
      #ÊIf, for some reason, we don't have a variation_id for the subsnp_id, skip this genotype
      next if (!defined($vs_subsnp_id));
      $variation_ids{$subsnp_id} = [$variation_id,$substrand_reversed_flag];
      $new_mappings = 1;
    }
    
    #ÊShould the alleles be flipped?
    my $reverse = ((($rev_alleles + $sub_strand + $variation_ids{$subsnp_id}->[1])%2 != 0) ? 1 : 0);
    
    #ÊIf any of the allele_ids were null, set the id to 0
    $allele_1 ||= 0;
    $allele_2 ||= 0;
    
    # Look up the alleles if necessary
    foreach my $allele_id ($allele_1,$allele_2) {
      if (!exists($alleles{$allele_id})) {
	$allele_sth->execute($allele_id);
	my ($a,$arev);
	$allele_sth->bind_columns(\$a,\$arev);
	$allele_sth->fetch();
	next if (!defined($a) || !defined($arev));
	$alleles{$allele_id} = [$a,$arev];
	$new_alleles = 1;
      }
    }
    #ÊLook up the sample id in the database if necessary
    if(!exists $all_individuals{$ind_id}{$ind_sub_id} ){
	if (exists $samples{$ind_id}){
	    $all_individuals{$ind_id}{$ind_sub_id} = $samples{$ind_id};
	}
	elsif($ind_id eq "MISS"){ 
            ## hack for non dbSNP-curated sample - enter or look up by submitter name 
	    print "Warning - missing sample checking for $subsnp_id  (submitted_ind_id $ind_sub_id) from handle $handle \n";

	    my ($ind_subname, $pop_subname);

	    ## get submitted name & pop from dbSNP 
	    $miss_ind_ext_sth->execute( $ind_sub_id)||die;
	    my $miss_dat = $miss_ind_ext_sth->fetchall_arrayref();
	    if(defined $miss_dat->[0]->[0]){
		$ind_subname = $miss_dat->[0]->[0];
	    }
	    else{
		die "Error looking up submitted_ind_id  $ind_sub_id\n";
	    }
	    if(defined $miss_dat->[0]->[1]){
		$pop_subname = $miss_dat->[0]->[1];
	    }

            ## look up or enter individual
	    my $individual_description =  $handle  . "|" . $ind_subname ;
	    $samind_ext_sth->execute( $ind_subname, $individual_description)||die;
	    my $new_sam_id = $samind_ext_sth->fetchall_arrayref();
	    unless(defined $new_sam_id->[0]->[0]){
		print "Failed to find ind $ind_subname desc =$individual_description from pop $pop_subname - entering\n";
		$sam_ins_sth->execute( $ind_subname, $individual_description)||die "Failed to enter non-curated sample:$!\n";
		$sam_ext_sth->execute( $ind_subname, $individual_description)||die;
		$new_sam_id = $sam_ext_sth->fetchall_arrayref();
		$ind_ins_sth->execute( $new_sam_id->[0]->[0] )||die;
		#print "Entered $ind_subname as ens id $new_sam_id->[0]->[0] \n";
		## check and enter population
		if(defined $pop_subname){
		   # print "Checking pop: $pop_subname\n";
		    my $population_description =  $handle  . "|" . $pop_subname ;
		    $sampop_ext_sth->execute($pop_subname, $population_description);
		    my $pop_id = $sampop_ext_sth->fetchall_arrayref();
		    unless(defined $pop_id->[0]->[0]){
			#print "Entering new pop: $pop_subname\n";
			#### enter population if missing
			$sam_ins_sth->execute( $pop_subname, $population_description)||die "Failed to enter non-curated sample:$!\n";
			$sam_ext_sth->execute( $pop_subname, $population_description)||die;	
			$pop_id = $sam_ext_sth->fetchall_arrayref();
			$pop_ins_sth->execute( $pop_id->[0]->[0])||die;	
		    }
		    #print "Linking ind:$ind_subname ($new_sam_id->[0]->[0])  to pop:$pop_subname ($pop_id->[0]->[0])\n";
		    $indpop_ins_sth->execute($new_sam_id->[0]->[0], $pop_id->[0]->[0])||die;	
		}	    		
	    }
	    $all_individuals{$ind_id}{$ind_sub_id} = $new_sam_id->[0]->[0];

	}

	else{
	    
	    $sample_sth->execute($ind_id);
	    my $sample_id;
	    $sample_sth->bind_columns(\$sample_id);
	    $sample_sth->fetch();
	     if (!defined($sample_id)){  ### check for further drop out??
		 warn "Skipping genotypes due to lack of sample id ind_id:$ind_id & ss$subsnp_id\n";
		 next;
	     } 
	    $samples{$ind_id} = $sample_id;
	    $all_individuals{$ind_id}{$ind_sub_id} = $sample_id;
	    $new_samples = 1;

	}
    }
    $allele_1 = $alleles{$allele_1}->[$reverse];
    $allele_2 = $alleles{$allele_2}->[$reverse];
    
    # Order the alleles in alphabetical order
    ($allele_1,$allele_2) = sort {uc($a) cmp uc($b)} ($allele_1,$allele_2);
    # Then in length order
    ($allele_1,$allele_2) = sort {length($a) <=> length($b)} ($allele_1,$allele_2);
    
    #ÊSkip this genotype if the alleles are N
    next if ($allele_1 eq 'N' && $allele_2 eq 'N');
    
    #my $row = join("\t",($variation_ids{$subsnp_id}->[0],$subsnp_id,$samples{$ind_id},$allele_1,$allele_2));
    my $row = join("\t",($variation_ids{$subsnp_id}->[0],$subsnp_id,$all_individuals{$ind_id}{$ind_sub_id},$allele_1,$allele_2));
    my $md5 = md5_hex($row);
    next if (exists($row_md5s{$md5}));
    $row_md5s{$md5}++;
    
    # Determine if this should go into the single or multiple genotype table
    if (length($allele_2) == 1 && length($allele_2) == 1 && $allele_1 ne '-' && $allele_2 ne '-') {
      push(@single_bp_gty,$row);
    }
    else {
      #ÊSkip this genotype if the first allele contains the string 'indeterminate'
      next if ($allele_1 =~ m/indeterminate/i);
      push(@multiple_bp_gty,$row);
    }
  }
  print $logh Progress::location();
  
  if (scalar(@single_bp_gty)) {
    # Open a file handle to the temp file that will be used for loading
    open(SGL,'>>',$loadfile) or die ("Could not open import file $loadfile for writing");
    #ÊGet a lock on the file
    flock(SGL,LOCK_EX);
    # Write the genotypes
    foreach my $gty (@single_bp_gty) {
      print SGL qq{$gty\n};
    }
    close(SGL);
    print $logh Progress::location();
  }
  
  if (scalar(@multiple_bp_gty)) {
    # Open a file handle to the file that will be used for loading the multi-bp genotypes
    open(MLT,'>>',$mlt_file) or die ("Could not open import file $mlt_file for writing");
    #ÊGet a lock on the file
    flock(MLT,LOCK_EX);
    # Write the genotypes
    foreach my $gty (@multiple_bp_gty) {
      print MLT qq{$gty\n};
    }
    close(MLT);
    print $logh Progress::location();
  }
  
  # If we had an allele file and we need to update it, do that
  delete($alleles{0});
  write_alleles($allele_file,\%alleles) if (defined($allele_file) && $new_alleles);
  print $logh Progress::location();
  # If we had a sample file and we need to update it, do that
  delete($samples{0});
  #write_samples($sample_file,\%samples) if (defined($sample_file) && $new_samples);
  # If we had a subsnp mapping file and we need to update it, do that
  delete($variation_ids{0});
  write_subsnp_mapping($mapping_file,\%variation_ids) if (defined($mapping_file) && $new_mappings);
  
  print $logh Progress::location() . "\tFinished dumping from $subind_table to $loadfile\n";
  return 1;
}

sub load_data_infile {
  my $self = shift;
  my $loadfile = shift;
  my $dst_table = shift;
  my @args = @_;
  unless($dst_table ){ die "Error from load_data_infile : destination table not set\n";}
  my $logh = $self->{'log'};  
  
  my $dbm = $self->{'db_manager'};
  my $dbVar = $dbm->dbVar();

  #ÊDisable the keys on the destination table before loading
  my $stmt = qq {
    ALTER TABLE
      $dst_table
    DISABLE KEYS
  };
  $dbVar->dbc()->do($stmt);
  print $logh Progress::location();

  # Load the data from the infile
  print $logh Progress::location() . "\tLoads the infile $loadfile to $dst_table\n";
  loadfile($loadfile,$dbVar->dbc(),$dst_table,@args);
  print $logh Progress::location();

  # Enable the keys after the inserts
  $stmt = qq {
    ALTER TABLE
      $dst_table
    ENABLE KEYS
  };
  $dbVar->dbc()->do($stmt);
  print $logh Progress::location();

  return 1;
}

sub read_alleles {
  my $allelefile = shift;
  
  my %alleles;
  
  # Get a filehandle
  open(FH,'<',$allelefile);
  
  # Lock the file for shared access
  flock(FH,LOCK_SH);
  # Parse the file
  while (<FH>) {
    chomp;
    my ($id,$a,$arev) = split;
    $alleles{$id} = [$a,$arev];
  }
  
  # Close the file
  close(FH);
  
  return \%alleles;
}
sub write_alleles {
  my $allele_file = shift;
  my $alleles = shift;

  # Get a filehandle  
  open(FH,'>',$allele_file);
  
  # Lock the file for exclusive access
  flock(FH,LOCK_EX);
  
  #ÊWrite each allele_id and the forward and reverse alleles
  while (my ($id,$a) = each(%{$alleles})) {
    print FH join("\t",($id,@{$a})) . "\n";
  }
  
  # Close the filehandle and release the lock
  close(FH);
}

sub read_subsnp_mapping {
  my $mappingfile = shift;
  
  # We just use the same method as the alleles since the data structures are identical
  return read_alleles($mappingfile);
}
sub write_subsnp_mapping {
  my $mappingfile = shift;
  my $mappings = shift;
  
  # We just use the same method as the alleles since the data structures are identical
  write_alleles($mappingfile,$mappings);
}

sub read_samples {
  my $samplefile = shift;
  
  my @samples;
  
  # Get a filehandle
  open(FH,'<',$samplefile);
  
  # Lock the file for shared access
  flock(FH,LOCK_SH);
  # Parse the file
  while (<FH>) {
    chomp;
    push(@samples,split);
  }
  
  # Close the file
  close(FH);
  
  my %s = (@samples);
  return \%s;
}
sub write_samples {
  my $sample_file = shift;
  my $samples = shift;

  # Get a filehandle  
  open(FH,'>',$sample_file);
  
  # Lock the file for exclusive access
  flock(FH,LOCK_EX);
  
  #ÊWrite each population_id/individual_id and the corresponding sample_id
  while (my @row = each(%{$samples})) {
    print FH join("\t",@row) . "\n";
  }
  
  # Close the filehandle and release the lock
  close(FH);
}

1;


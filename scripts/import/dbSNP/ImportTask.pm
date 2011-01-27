use strict;
use warnings;

#generic object for the dbSNP data. Contains the general methods to dump the data into the new Variation database. Any change in the methods
# will need to overload the correspondent method in the subclass for the specie

package dbSNP::ImportTask;

use ImportUtils qw(dumpSQL debug create_and_load load loadfile get_create_statement);
use Progress;
use DBI qw(:sql_types);
use Fcntl qw( LOCK_SH LOCK_EX );
use Digest::MD5 qw ( md5_hex );

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


sub allele_table {
  my $self = shift;
  my $loadfile = shift;
  my $task_start = shift;
  my $task_end = shift;
  my $allelefile = shift;
  my $samplefile = shift;
  
  my $stmt;
  my $dbm = $self->{'db_manager'};
  my $shared_db = $dbm->dbSNP_shared();
  
  my $logh = $self->{'log'};
  print $logh Progress::location() . "\tImporting for SubSNP ids from $task_start to $task_end\n";
  
  #�Prepared statement for getting the SubSNPs with or without population frequency data
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
      $shared_db..ObsVariation ov ON (
	ov.var_id = ss.variation_id
      ) JOIN
      $shared_db..UniVariAllele uv ON (
	uv.univar_id = ov.univar_id
      )
    WHERE
      ss.subsnp_id BETWEEN ? AND ?
    ORDER BY
      ss.subsnp_id ASC,
      uv.allele_id ASC,
      b.pop_id ASC
  };
  my $ss_sth = $dbm->dbSNP()->dbc->prepare($stmt);
  
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
  my $ss_freq_sth = $dbm->dbSNP()->dbc->prepare($stmt);
  
  #�Prepared statement to get the corresponding variation_ids for a range of subsnps from variation_synonym
  $stmt = qq{
    SELECT
      vs.subsnp_id,
      vs.variation_id
    FROM
      variation_synonym vs
    WHERE
      vs.subsnp_id BETWEEN ? AND ?
  };
  my $vs_sth = $dbm->dbVar()->dbc->prepare($stmt);

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
  my $allele_sth = $dbm->dbSNP()->dbc->prepare($stmt);
  
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
  my $sample_sth = $dbm->dbVar()->dbc->prepare($stmt);
  
  #�Hash to hold the alleles in memory
  my %alleles;
  #�If the allele file exist, we'll read alleles from it
  %alleles = %{read_alleles($allelefile)} if (defined($allelefile) && -e $allelefile);
  print $logh Progress::location();
  # The allele_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value
  $alleles{0} = ['\N','\N'];
  # Keep track if we did any new lookups
  my $new_alleles = 0;
  
  #�Hash to keep sample_ids in memory
  my %samples;
  #�If the sample file exist, we'll read alleles from it
  %samples = %{read_samples($samplefile)} if (defined($samplefile) && -e $samplefile);
  print $logh Progress::location();
  # The pop_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value
  $samples{0} = '\N';
  # Keep track if we did any new lookups
  my $new_samples = 0;
  
  # Hash to hold the output until it's time to write to disk
  my @output;
  
  # We will start at the beginning of the range this task was given
  my $min_ss = $task_start;
  my $max_ss = $task_end;
  
  # First, get all SubSNPs
  $ss_sth->bind_param(1,$min_ss,SQL_INTEGER);
  $ss_sth->bind_param(2,$max_ss,SQL_INTEGER);
  $ss_freq_sth->bind_param(1,$min_ss,SQL_INTEGER);
  $ss_freq_sth->bind_param(2,$max_ss,SQL_INTEGER);
  $vs_sth->bind_param(1,$min_ss,SQL_INTEGER);
  $vs_sth->bind_param(2,$max_ss,SQL_INTEGER);
  
  #�Fetch the query result as a hashref
  $ss_sth->execute();
  print $logh Progress::location();
  my $subsnp_alleles = $ss_sth->fetchall_hashref(['subsnp_id','pop_id','allele_id']);
  
  # Fetch the population frequency alleles as an arrayref
  $ss_freq_sth->execute();
  print $logh Progress::location();
  my $subsnp_freq_alleles = $ss_freq_sth->fetchall_arrayref();
  
  # Add new alleles with frequency data and update the ones with data missing in the first hash (where subsnp_id, pop_id and allele_id match).
  map {
    $subsnp_alleles->{$_->[0]}{$_->[1]}{$_->[2]} = {'substrand_reversed_flag' => $_->[3], 'frequency' => $_->[4], 'count' => $_->[5]};
  } @{$subsnp_freq_alleles};
  
  # Fetch the subsnp_id to variation_id mapping and store as hashref
  $vs_sth->execute();
  print $logh Progress::location();
  my $variation_ids = $vs_sth->fetchall_hashref(['subsnp_id']);
  
  #�Now, loop over the subsnp hash and print it to the tempfile so we can import the data. Replace the allele_id with the corresponding allele on-the-fly
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
	    $new_alleles = 1;
	  }
	  #�Look up the sample id in the database if necessary
	  if (!exists($samples{$pop_id})) {
	    $sample_sth->execute($pop_id);
	    my $sample_id;
	    $sample_sth->bind_columns(\$sample_id);
	    $sample_sth->fetch();
	    # If the pop_id was not found (dbSNP actually has broken relationships), set it to null
	    $sample_id ||= '\N';
	    $samples{$pop_id} = $sample_id;
	    $new_samples = 1;
	  }
	  my $row = join("\t",($variation_ids->{$subsnp_id}{'variation_id'},$subsnp_id,$samples{$pop_id},$alleles{$allele_id}->[$allele_data->{'substrand_reversed_flag'}],$allele_data->{'frequency'},$allele_data->{'count'}));
	  push(@output,$row);
	}
    }
  }

  print $logh Progress::location();
  open(IMP,'>>',$loadfile);
  # Lock the file
  flock(IMP,LOCK_EX);
  foreach my $row (@output) {
    print IMP $row . "\n";
  }
  close(IMP);
  print $logh Progress::location();
  
  write_alleles($allelefile,\%alleles) if (defined($allelefile) && $new_alleles);
  print $logh Progress::location();
  delete($samples{0});
  write_samples($samplefile,\%samples) if (defined($samplefile) && $new_samples);
  print $logh Progress::location();
  
  # We're done, return success
  return 1;
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
  
  #�Put the log filehandle in a local variable
  my $logh = $self->{'log'};
  print $logh Progress::location() . "\tDumping from $subind_table to $loadfile, processing genotype data for submitted_ind_ids $start - $end\n";  
  
  my $dbm = $self->{'db_manager'};
  my $shared_db = $dbm->dbSNP_shared();
  my $dbVar = $dbm->dbVar();
  my $dbSNP = $dbm->dbSNP();
  
  #�A prepared statement for getting the genotype data from the dbSNP mirror. Note the chr_num, dbSNP is storing the gty allele for each chromosome copy? Anyway, we'll get duplicated rows if not specifying this (or using distinct)
  #�Order the results by subsnp_id so that we can do the subsnp_id -> variation_id lookup more efficiently
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
      LEN(ug.gty_str) AS pattern_length
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
      )
    WHERE
      si.submitted_ind_id BETWEEN ? AND ? AND
      ga.chr_num = 1
    ORDER BY
      si.subsnp_id ASC
  };
  my $subind_sth = $dbSNP->dbc()->prepare($stmt);
  
  #�Prepared statement to get the corresponding variation_ids for a range of subsnps from variation_synonym
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
  
  #�Hash to hold the alleles in memory
  my %alleles;
  #�If the allele file exist, we'll read alleles from it
  %alleles = %{read_alleles($allele_file)} if (defined($allele_file) && -e $allele_file);
  print $logh Progress::location();
  # The allele_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value.
  $alleles{0} = ['\N','\N'];
  # Keep track if we did any new lookups
  my $new_alleles = 0;
  
  #�Hash to keep sample_ids in memory
  my %samples;
  #�If the sample file exist, we'll read alleles from it
  %samples = %{read_samples($sample_file)} if (defined($sample_file) && -e $sample_file);
  print $logh Progress::location();
  # The individual_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value
  $samples{0} = '\N';
  # Keep track if we did any new lookups
  my $new_samples = 0;
  
  #�Hash to keep subsnp_id to variation_id mappings in memory
  my %variation_ids;
  #�If the mapping file exist, we'll read mappings from it
  %variation_ids = %{read_subsnp_mapping($mapping_file)} if (defined($mapping_file) && -e $mapping_file);
  print $logh Progress::location();
  # The subsnp_id = 0 is a replacement for NULL but since it's used for a key in the hash below, we need it to have an actual numerical value
  $variation_ids{0} = ['\N','\N'];
  # Keep track if we did any new lookups
  my $new_mappings = 0;
  
  #�Keep the genotypes in an array and write them in the end to avoid locking the file more than necessary. Separate signle and multiple bp genotypes. Use a hash having the md5sum of the row to determine whether we've already used this one 
  my @single_bp_gty;
  my @multiple_bp_gty;
  my %row_md5s;
  
  # First, get all SubSNPs
  $subind_sth->bind_param(1,$start,SQL_INTEGER);
  $subind_sth->bind_param(2,$end,SQL_INTEGER);
  
  # Fetch the import data as an arrayref
  $subind_sth->execute();
  print $logh Progress::location();
  
  #�Now, loop over the import data and print it to the tempfile so we can import the data. Replace the allele_id with the corresponding allele on-the-fly
  while (my $data = $subind_sth->fetchrow_arrayref()) {
    my $subsnp_id = $data->[0];
    my $ind_id = $data->[1];
    my $allele_1 = $data->[2];
    my $allele_2 = $data->[3];
    my $sub_strand = $data->[4];
    my $rev_alleles = $data->[5];
    my $pattern_length = $data->[6];
    
    # If pattern length is less than 3, skip this genotype
    next if ($pattern_length < 3);
    # If ind_id is null (there is a case where this happens because one individual in SubmittedIndividual is missing from Individual (131)), skip this genotype.
    next if (!defined($ind_id));

    # Look up the variation_id if necessary. This should be slow for the first chunk on each chromosome but fast for the rest..
    if (!exists($variation_ids{$subsnp_id})) {
      $vs_sth->bind_param(1,$subsnp_id,SQL_INTEGER);
      $vs_sth->execute();
      my ($vs_subsnp_id,$variation_id,$substrand_reversed_flag);
      $vs_sth->bind_columns(\$vs_subsnp_id,\$variation_id,\$substrand_reversed_flag);
      $vs_sth->fetch();
      #�If, for some reason, we don't have a variation_id for the subsnp_id, skip this genotype
      next if (!defined($vs_subsnp_id));
      $variation_ids{$subsnp_id} = [$variation_id,$substrand_reversed_flag];
      $new_mappings = 1;
    }
    
    #�Should the alleles be flipped?
    my $reverse = ((($rev_alleles + $sub_strand + $variation_ids{$subsnp_id}->[1])%2 != 0) ? 1 : 0);
    
    #�If any of the allele_ids were null, set the id to 0
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
    #�Look up the sample id in the database if necessary
    if (!exists($samples{$ind_id})) {
      $sample_sth->execute($ind_id);
      my $sample_id;
      $sample_sth->bind_columns(\$sample_id);
      $sample_sth->fetch();
      next if (!defined($sample_id));
      $samples{$ind_id} = $sample_id;
      $new_samples = 1;
    }
    
    $allele_1 = $alleles{$allele_1}->[$reverse];
    $allele_2 = $alleles{$allele_2}->[$reverse];
    
    # Order the alleles in alphabetical order
    ($allele_1,$allele_2) = sort {uc($a) cmp uc($b)} ($allele_1,$allele_2);
    # Then in length order
    ($allele_1,$allele_2) = sort {length($a) <=> length($b)} ($allele_1,$allele_2);
    
    #�Skip this genotype if the alleles are N
    next if ($allele_1 eq 'N' && $allele_2 eq 'N');
    
    my $row = join("\t",($variation_ids{$subsnp_id}->[0],$subsnp_id,$samples{$ind_id},$allele_1,$allele_2));
    my $md5 = md5_hex($row);
    next if (exists($row_md5s{$md5}));
    $row_md5s{$md5}++;
    
    # Determine if this should go into the single or multiple genotype table
    if (length($allele_2) == 1 && length($allele_2) == 1 && $allele_1 ne '-' && $allele_2 ne '-') {
      push(@single_bp_gty,$row);
    }
    else {
      #�Skip this genotype if the first allele contains the string 'indeterminate'
      next if ($allele_1 =~ m/indeterminate/i);
      push(@multiple_bp_gty,$row);
    }
  }
  print $logh Progress::location();
  
  if (scalar(@single_bp_gty)) {
    # Open a file handle to the temp file that will be used for loading
    open(SGL,'>>',$loadfile) or die ("Could not open import file $loadfile for writing");
    #�Get a lock on the file
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
    #�Get a lock on the file
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
  write_samples($sample_file,\%samples) if (defined($sample_file) && $new_samples);
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

  my $logh = $self->{'log'};  
  
  my $dbm = $self->{'db_manager'};
  my $dbVar = $dbm->dbVar();
  
  #�Disable the keys on the destination table before loading
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
  
  #�Write each allele_id and the forward and reverse alleles
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
  
  #�Write each population_id/individual_id and the corresponding sample_id
  while (my @row = each(%{$samples})) {
    print FH join("\t",@row) . "\n";
  }
  
  # Close the filehandle and release the lock
  close(FH);
}

1;


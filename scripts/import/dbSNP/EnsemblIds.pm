use strict;
use warnings;
#object that contains the specific methods to dump data when the specie is a HUMAN (adds HGVbase and TSC information). 
package dbSNP::EnsemblIds;

#use dbSNP::GenericContig;
use dbSNP::GenericChromosome;

use vars qw(@ISA);
use ImportUtils qw(debug load dumpSQL create_and_load);
use Progress;

#@ISA = ('dbSNP::GenericContig');
@ISA = ('dbSNP::GenericChromosome');

sub dump_dbSNP{
    my $self = shift;
    
    my $start;
    my $end;
    my $duration;
    
    #first, dump all dbSNP data as usual
    $start = time();
    $self->SUPER::dump_dbSNP();
    print Progress::location();
    $end = time();
    $duration = Progress::time_format($end-$start);
    print $duration->{'weeks'} . " weeks, " . $duration->{'days'} . " days, " . $duration->{'hours'} . " hours, " . $duration->{'minutes'} . " minutes and " . $duration->{'seconds'} . " seconds spent in SUPER::dump_dbSNP()\n";
    
    #then, get ENS IDs from dbSNP
    $start = time();
    $self->dump_ENSIDs() if $self->{'dbCore'}->species =~ /hum|homo|rat|mouse|platypus|tetraodon/i;
    print Progress::location();
    $end = time();
    $duration = Progress::time_format($end-$start);
    print $duration->{'weeks'} . " weeks, " . $duration->{'days'} . " days, " . $duration->{'hours'} . " hours, " . $duration->{'minutes'} . " minutes and " . $duration->{'seconds'} . " seconds spent in dump_ENSIDs()\n";
    
    $start = time();
    $self->dump_AFFYIDs() if $self->{'dbCore'}->species =~ /hum|homo/i;
    print Progress::location();
    $end = time();
    $duration = Progress::time_format($end-$start);
    print $duration->{'weeks'} . " weeks, " . $duration->{'days'} . " days, " . $duration->{'hours'} . " hours, " . $duration->{'minutes'} . " minutes and " . $duration->{'seconds'} . " seconds spent in dump_AFFYIDs()\n";
  
}

#get ENSIDs from dbSNP
sub dump_ENSIDs{
  my $self = shift;

  my %ens_batch_name = ('rat' => ['RAT_COMPUTATIONAL_CELERA','RAT_COMPUTATIONAL_STAR'],
			'platypus' => ['2008FEB_platypus_assembly','2008FEB_platypus_reads'],
		       );

  #add the ENS source to the table
  foreach my $species (keys %ens_batch_name) {
    if ($self->{'dbCore'}->species =~ /$species/i) {
      foreach my $batch_name (@{$ens_batch_name{$species}}) {
	debug("Processing $batch_name...");
	$self->{'dbVar'}->do(qq{INSERT INTO source (name,version,description) values ("ENSEMBL_$batch_name",1,"Variation called by ENSEMBL")});
	my $source_id = $self->{'dbVar'}->db_handle->{'mysql_insertid'}; #get the last autoinc id in the database (the one from the ENS source)
	#and finally add the ENS ids to the synonyms table
	debug("Dumping ENS information from dbSNP");
	
	my $stmt = qq{
	  SELECT
	    'rs'+LTRIM(STR(s.snp_id)),
	    s.subsnp_id,
	    $source_id,
	    s.loc_snp_id
	  FROM
	    SubSNP s,
	    Batch b
	  WHERE
	    s.batch_id = b.batch_id AND
	    b.handle = 'ENSEMBL' AND
	    b.loc_batch_id = '$batch_name'
	};
	dumpSQL($self->{'dbSNP'},$stmt);
	debug("Loading ENS ids into temporary table");
	create_and_load($self->{'dbVar'},"tmp_rs_ENS_$batch_name","rsID *","ss_id i","source_id i","ENSid");
	debug("Loading ENS ids into variation_synonym table");
	$self->{'dbVar'}->do(qq{ INSERT IGNORE INTO variation_synonym (variation_id, subsnp_id, source_id, name)
				     SELECT v.variation_id, trt.ss_id, trt.source_id, trt.ENSid 
				     FROM variation v, tmp_rs_ENS_$batch_name trt
				     WHERE v.name = trt.rsID 
				}
			    );
	#and finally, remove the temporary table
	#$self->{'dbVar'}->do(qq{DROP TABLE tmp_rs_ENS_$batch_name});
      }
    }
  }
}


# This method uses hardcoded databases on ens-genomics1
sub dump_AFFYIDs{

  my $self = shift;
  my ($source_name,$set_name);
    
  debug(localtime() . "\tCreate separate db connection to ens-genomics\n");
  my $dbh = DBI->connect("DBI:mysql:timeout=2678200;host=ens-genomics1","ensro","") or die("Could not connect to database $DBI::errstr");
  my $stmt;
  
  debug("Dumping AFFY information from dbSNP");
  $stmt = qq{
    SELECT
	'rs'+LTRIM(STR(s.snp_id)),
	s.loc_snp_id,
	loc_batch_id
    FROM
	SubSNP s,
	Batch b
    WHERE
	s.batch_id = b.batch_id AND
	b.handle = 'AFFY'
  };
  dumpSQL($self->{'dbSNP'},$stmt);
  
  debug("Loading  ids into temporary table");
  create_and_load($self->{'dbVar'},"tmp_rs_AFFY","rsID *","AFFYid", "affy_name");
  print Progress::location();
  
  foreach my $table ("yuan_aff_100k_var_46","yuan_aff_500k_var_46","yuan_aff_genome6_var_47") {
    
    $stmt = qq{
	USE $table
    };
    $dbh->do($stmt);
    
    if ($table =~ /100k/i) {
      $source_name = "Affy GeneChip 100K Array";
      $set_name = "Mapping50K";
    }
    elsif ($table =~ /500k/i) {
      $source_name = "Affy GeneChip 500K Array";
      $set_name = "Mapping250K";
    }
    elsif ($table =~ /genome6/i) {
      $source_name = "Affy GenomeWideSNP_6.0";
      $set_name = "6.0";
    }

    debug(localtime() . "\tGet name_pair create statements\n");
    $stmt = qq{
	SHOW CREATE TABLE
	    name_pair
    };
    my $create_stmt = $dbh->selectall_arrayref($stmt)->[0][1];

    debug(localtime() . "\tCreating name_pair table with source_name $source_name...");
    $self->{'dbVar'}->do($create_stmt);
    print Progress::location();
    $stmt = qq{
	RENAME TABLE
	    name_pair
	TO
	    $table\_name_pair
    };
    $self->{'dbVar'}->do($stmt);
    print Progress::location();
    
    debug(localtime() . "\tDumping AFFY id name pairs from ens-genomics\n");
    $stmt = qq{
	SELECT
	    *
	FROM
	    name_pair
    };
    dumpSQL($dbh,$stmt);
    
    debug(localtime() . "\tLoading AFFY id name pairs into $table\_name_pair\n");
    load($self->{'dbVar'},"$table\_name_pair","affy_name","rs_name");
    print Progress::location();
    
    debug(localtime() . "\tCreating a temporary table for AFFY ids\n");
    $stmt = qq{
	CREATE TABLE
	    tmp
	SELECT
	    a.AFFYid AS affy_name,
	    a.rsID AS rs_name
	FROM
	    tmp_rs_AFFY a,
	    $table\_name_pair c
	WHERE
	    a.AFFYid = c.affy_name
    };
    $self->{'dbVar'}->do($stmt);
    print Progress::location();
    
    debug(localtime() . "\tAdding AFFY ids from temporary table\n");
    $stmt = qq{
	INSERT IGNORE INTO
	    $table\_name_pair (
		affy_name,
		rs_name
	    )
	SELECT
	    affy_name,
	    rs_name
	FROM
	    tmp
    };
    $self->{'dbVar'}->do($stmt);
    print Progress::location();
    
    debug(localtime() . "\tDropping temporary table\n");
    $stmt = qq{
	DROP TABLE
	    tmp
    };
    $self->{'dbVar'}->do($stmt);
    print Progress::location();
    
    
    # The code below has been replaced by the dump-and-load code above
    #$self->{'dbVar'}->do(qq{CREATE TABLE $table\_name_pair like $table.name_pair});
    #$self->{'dbVar'}->do(qq{insert into $table\_name_pair select * from $table.name_pair});
    #$self->{'dbVar'}->do(qq{insert ignore into $table\_name_pair
    #                        select a.AFFYid as affy_name,a.rsID as rs_name 
    #                        from tmp_rs_AFFY a, $table.name_pair c 
    #                        where a.AFFYid=c.affy_name});

    if ($table =~ /genome6/i) {
	debug(localtime() . "\tUpdating $table\_name_pair\n");
	$stmt = qq{
	    UPDATE
		$table\_name_pair
	    SET
		affy_name = REPLACE(affy_name,'AFFY_6_1M_','')
	};
	$self->{'dbVar'}->do($stmt);
	print Progress::location();
    }
    
    my $source_id_ref = $self->{'dbVar'}->db_handle->selectall_arrayref(qq{
                     SELECT source_id from source where name = "$source_name"});
    my $source_id = $source_id_ref->[0][0];

    if (!$source_id) {
      $self->{'dbVar'}->do(qq{insert into source (name) values("$source_name")});
	print Progress::location();
      $source_id = $self->{'dbVar'}->db_handle->{'mysql_insertid'};
    }

    debug("Inserting in variation_synonym table from $table\_name_pair...");
    $self->{'dbVar'}->do(qq{ INSERT IGNORE INTO variation_synonym (variation_id, source_id, name)
                                     SELECT v.variation_id, $source_id as source_id, t.affy_name as name 
                                     FROM variation v, $table\_name_pair t
                                     WHERE v.name = t.rs_name
                                }
			);
    print Progress::location();

    #update rs_ID to rsCurrent from rsHigh
    #$self->{'dbVar'}->do(qq{update tmp_rs_AFFY_test t, rsHist h set t.rsID=h.rsCurrent where t.rsID=h.rsHigh});
    debug("Inserting in variation_synonym table from table  tmp_rs_AFFY with $source_name...");
    $self->{'dbVar'}->do(qq{ INSERT IGNORE INTO variation_synonym (variation_id, source_id, name)
				     SELECT v.variation_id, $source_id as source_id, t.AFFYid as name 
				     FROM variation v, tmp_rs_AFFY t
				     WHERE v.name = t.rsID
                                     AND t.affy_name like "%$set_name%"
				}
			);
    print Progress::location();

   #and finally, remove the temporary table
   #$self->{'dbVar'}->do(qq{DROP TABLE $table\_name_pair});

  }
}

1;

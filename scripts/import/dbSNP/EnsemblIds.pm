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
    #$self->dump_ENSIDs() if $self->{'dbCore'}->species =~ /hum|homo|rat|mouse|platypus|tetraodon/i;
    print Progress::location();
    $end = time();
    $duration = Progress::time_format($end-$start);
    print $duration->{'weeks'} . " weeks, " . $duration->{'days'} . " days, " . $duration->{'hours'} . " hours, " . $duration->{'minutes'} . " minutes and " . $duration->{'seconds'} . " seconds spent in dump_ENSIDs()\n";
    
    # Note that the dump_AFFIDs routine has some database/table names hard-coded within!!
		# /!\ Now has to be launched in a separate script for post-processing: import_affy.pl
    $start = time();
#    $self->dump_AFFYIDs() if $self->{'dbm'}->dbCore()->species =~ /hum|homo/i;
    print Progress::location();
    $end = time();
    $duration = Progress::time_format($end-$start);
    print $duration->{'weeks'} . " weeks, " . $duration->{'days'} . " days, " . $duration->{'hours'} . " hours, " . $duration->{'minutes'} . " minutes and " . $duration->{'seconds'} . " seconds spent in dump_AFFYIDs()\n";
	
	$start = time();
    $self->dump_LSDBIDs() if $self->{'dbm'}->dbCore()->species =~ /hum|homo/i;
    print Progress::location();
    $end = time();
    $duration = Progress::time_format($end-$start);
    print $duration->{'weeks'} . " weeks, " . $duration->{'days'} . " days, " . $duration->{'hours'} . " hours, " . $duration->{'minutes'} . " minutes and " . $duration->{'seconds'} . " seconds spent in dump_LSDBIDs()\n";
	
  
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


# This method uses the pontus_affy_array_mapping database on ens-variation
sub dump_AFFYIDs{

  my $self = shift;
  my ($source_name,$source_description,$source_url,$set_name);
    
  my $stmt;
  $source_url = "http://www.affymetrix.com/";
  
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
  
  my $db = "pontus_dbsnp132_human_external_data";
  foreach my $table ("affy_array_name_pair_100k","affy_array_name_pair_500k","affy_array_name_pair_g6") {
    
    if ($table =~ /100k/i) {
      $source_name = "Affy GeneChip 100K Array";
      $source_description = "Variants from the Affymetrix GeneChip Human Mapping 100K Array Set";
      $set_name = "Mapping50K";
    }
    elsif ($table =~ /500k/i) {
      $source_name = "Affy GeneChip 500K Array";
      $source_description = "Variants from the Affymetrix GeneChip Human Mapping 500K Array Set";
      $set_name = "Mapping250K";
    }
    elsif ($table =~ /g6/i) {
      $source_name = "Affy GenomeWideSNP_6.0";
      $source_description = "Variants from the Affymetrix Genome-Wide Human SNP Array 6.0";
      $set_name = "6.0";
    }

    debug(localtime() . "\tCreating name_pair table with source_name $source_name...");
    $stmt = qq{
	CREATE TABLE
	    $table
	LIKE
	    $db\.$table
    };
    $self->{'dbVar'}->do($stmt);
    print Progress::location();
    
    debug(localtime() . "\tDumping AFFY id name pairs from $db\.$table\n");
    $stmt = qq{
	INSERT INTO
	    $table (
		affy_name,
		rs_name
	    )
	SELECT
	    affy_name,
	    rs_name
	FROM
	    $db\.$table
    };
    $self->{'dbVar'}->do($stmt);
    
    debug(localtime() . "\tCreating a temporary table for AFFY ids\n");
    $stmt = qq{
	CREATE TABLE
	    tmp
	SELECT
	    a.AFFYid AS affy_name,
	    a.rsID AS rs_name
	FROM
	    tmp_rs_AFFY a,
	    $table c
	WHERE
	    a.AFFYid = c.affy_name
    };
    $self->{'dbVar'}->do($stmt);
    print Progress::location();
    
    debug(localtime() . "\tAdding AFFY ids from temporary table\n");
    $stmt = qq{
	INSERT IGNORE INTO
	    $table (
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

    if ($table =~ /g6/i) {
	debug(localtime() . "\tUpdating $table\n");
	$stmt = qq{
	    UPDATE
		$table
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
      $self->{'dbVar'}->do(qq{insert into source (name,description,url,type) values("$source_name","$source_description","$source_url","chip")});
	print Progress::location();
      $source_id = $self->{'dbVar'}->db_handle->{'mysql_insertid'};
    }

    debug("Inserting in variation_synonym table from $table\...");
    $self->{'dbVar'}->do(qq{ INSERT IGNORE INTO variation_synonym (variation_id, source_id, name)
                                     SELECT v.variation_id, $source_id as source_id, t.affy_name as name 
                                     FROM variation v, $table t
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
    $self->{'dbVar'}->do(qq{DROP TABLE $table});
    print Progress::location();

  }
}


# gets LSDB local IDs
sub dump_LSDBIDs {
	my $self = shift;
	
    print Progress::location();
    debug(localtime() . "\tAdding/getting source ID for LSDB\n");
	
	my $source_name = 'LSDB';
	my $source_id_ref = $self->{'dbVar'}->db_handle->selectall_arrayref(qq{
		SELECT source_id from source where name = "$source_name"
	});
    my $source_id = $source_id_ref->[0][0];
	
	# add source row if not there
	if (!$source_id) {

	    my ($species,$tax_id,$version) = $self->{'snp_dbname'} =~ m/^(.+)?\_([0-9]+)\_([0-9]+)$/;

		$self->{'dbVar'}->do(qq{insert into source (name, version, description,url) 
                                        values("$source_name",
                                                $version,
                                               "Variants dbSNP annotates as being from LSDBs",
                                               "http://www.ncbi.nlm.nih.gov/projects/SNP/")});

		$source_id = $self->{'dbVar'}->db_handle->{'mysql_insertid'};
	}
	
    print Progress::location();
    debug(localtime() . "\tDumping LSDB local IDs from dbSNP\n");
	
	# dump IDs from dbSNP DB
	my $stmt = qq{
		SELECT
			DISTINCT ss.loc_snp_id,
			ss.subsnp_id,
			sl.snp_id
			
		FROM SubSNP ss
			JOIN SNPSubSNPLink sl ON sl.subsnp_id = ss.subsnp_id
			JOIN Batch b ON b.batch_id = ss.batch_id
			JOIN $self->{'dbSNP_share_db'}..Method m ON b.method_id = m.method_id AND m.method_class = 109;
	};
	dumpSQL($self->{'dbSNP'}, $stmt);
	
    print Progress::location();
    debug(localtime() . "\tCopying IDs to synonym table\n");
	
	create_and_load($self->{'dbVar'}, "tmp_lsdb_ids", 'name', 'ssid i*', 'rsid i*');
	
	$self->{'dbVar'}->db_handle->do(qq{
		INSERT IGNORE INTO
			variation_synonym(variation_id, subsnp_id, source_id, name)
		SELECT
			v.variation_id,
			t.ssid,
			$source_id,
			t.name
		FROM
			variation v,
			tmp_lsdb_ids t
		WHERE
			v.name = concat('rs', t.rsid);
	});
	
    print Progress::location();
    debug(localtime() . "\tDropping temporary table\n");
	
	$self->{'dbVar'}->db_handle->do(qq{DROP TABLE tmp_lsdb_ids});
}

1;

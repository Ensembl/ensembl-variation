use strict;
use warnings;

#generic object for the dbSNP data. Contains the general methods to dump the data into the new Variation database. Any change in the methods
# will need to overload the correspondent method in the subclass for the specie

package dbSNP::GenericContig;

use POSIX;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use ImportUtils qw(dumpSQL debug create_and_load load get_create_statement);
use Progress;
use DBI qw(:sql_types);

#creates the object and assign the attributes to it (connections, basically)
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbSNP, $dbCore, $dbVar, $snp_dbname, $species, $tmp_dir, $tmp_file, $limit, $mapping_file_dir, $dbSNP_BUILD_VERSION, $shared_db, $ASSEMBLY_VERSION, $skip_routines, $logh) =
        rearrange([qw(DBSNP DBCORE DBVAR SNP_DBNAME SPECIES TMPDIR TMPFILE LIMIT MAPPING_FILE_DIR DBSNP_VERSION SHARED_DB ASSEMBLY_VERSION SKIP_ROUTINES LOG)],@_);

#  my $dbSNP_share_db = "dbSNP_$dbSNP_BUILD_VERSION\_shared";
#  $dbSNP_share_db =~ s/dbSNP_b/dbSNP_/;
#  my $dbSNP_share_db = 'new_shared';
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
		'log' => $logh}, $class;
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
    'individual_genotypes',
    'population_genotypes',
    'allele_table',
    'flanking_sequence_table',
    'variation_feature',
    'cleanup'
  );
  
  #ÊThe GenericContig object has an array where routines that should be skipped can be specified. For now, add create_coredb and cleanup by default
  push(@{$self->{'skip_routines'}},('create_coredb','cleanup'));
  
  # This is just for experimenting, should be removed in any "real" import
  # push(@{$self->{'skip_routines'}},('allele_table','individual_genotypes','population_genotypes','flanking_sequence_table','variation_feature'));
  
  #ÊThis is for resuming if the script crashed in allele_table
  # push(@{$self->{'skip_routines'}},('source_table','population_table','individual_table','variation_table'));
  
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
    $self->$subroutine();
    print $logh $clock->duration();
  }

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
    
    $self->{'dbVar'}->do(qq{INSERT INTO source (source_id,name,version,description) VALUES (1,"$dbname",$version,"Variants (including SNPs and indels) imported from dbSNP [http://www.ncbi.nlm.nih.gov/projects/SNP/]")});
    print $logh Progress::location();

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
  my $count = $self->{'dbSNP'}->selectall_arrayref($stmt);
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
    
# create table rsHist to store rs history
# Table vwRsMergeArch only exists for human. Filter this as a temporary solution and talk to Yuan about it
    if ($self->{'dbCore'}->species =~ m/human|homo/i) {
      dumpSQL($self->{'dbSNP'},(qq{SELECT * FROM vwRsMergeArch})) unless ($resume_at_subsnp_id > 0);
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
   dumpSQL($self->{'dbSNP'},$stmt);
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
  if ($self->{'dbCore'}->species =~ /homo|pan|anoph/i) {
    $individual_type_id = 3;
  }
  elsif ($self->{'dbCore'}->species =~ /mus/i) {
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
  my ($min_ss,$max_ss) = @{$self->{'dbSNP'}->selectall_arrayref($stmt)->[0]};
  
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
  if($self->{'dbCore'}->species =~ /human|homo/i) {
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
  $assembly_version = 1 if $self->{'dbCore'}->species =~ /ornith/i;
  $assembly_version = 3 if $self->{'dbCore'}->species =~ /rerio/i;

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

#
# loads variation_group and variation_group_variation tables from the
# contents of the HapSet and HapSetSnpList tables
#
sub variation_group {
    my $self = shift;

  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  debug(localtime() . "\tDumping HapSet data");

  dumpSQL($self->{'dbSNP'}, qq{SELECT  CONCAT(hs.handle, ':', hs.hapset_name),
                    hs.hapset_id, hssl.subsnp_id
             FROM HapSet hs, HapSetSnpList hssl   ###, SubSNP ss 
             WHERE hs.hapset_id = hssl.hapset_id});

  create_and_load($self->{'dbVar'}, 'tmp_var_grp', 'name', 'hapset_id i*', 'subsnp_id i*');
  print $logh Progress::location();

  $self->{'dbVar'}->do("ALTER TABLE variation_group add column hapset_id int");
  print $logh Progress::location();

  debug(localtime() . "\tLoading variation_group");

  $self->{'dbVar'}->do(qq{INSERT INTO variation_group (name, source_id, type, hapset_id)
                SELECT name, 1, 'haplotype', hapset_id
                FROM tmp_var_grp
                GROUP BY hapset_id});
  print $logh Progress::location();

  $self->{'dbVar'}->do("ALTER TABLE variation_group ADD INDEX hapset_id(hapset_id)");
  print $logh Progress::location();

  debug(localtime() . "\tLoading variation_group_variation");

  # there are a few hapsets in dbSNP which have two subsnps which have been
  # merged into the same refsnp.  Thus the group by clause.

  $self->{'dbVar'}->do(qq{INSERT INTO variation_group_variation
                     (variation_group_id, variation_id)
                SELECT vg.variation_group_id, vs.variation_id
                FROM   variation_group vg, variation_synonym vs,
                       tmp_var_grp tvg
                WHERE  tvg.hapset_id = vg.hapset_id
                AND    tvg.subsnp_id = vs.subsnp_id
                GROUP BY variation_group_id, variation_id});
  print $logh Progress::location();

  $self->{'dbVar'}->do("DROP TABLE tmp_var_grp");
  print $logh Progress::location();
}

#
# loads allele_group table
#
sub allele_group {
    my $self = shift;

  #ÊPut the log filehandle in a local variable
  my $logh = $self->{'log'};
  
  debug(localtime() . "\tDumping Hap data");

  dumpSQL($self->{'dbSNP'}, qq{SELECT  h.hap_id, h.hapset_id, h.loc_hap_id,
                    hsa.snp_allele, hsa.subsnp_id
             FROM   Hap h, HapSnpAllele hsa   ###, SubSNP ss
             WHERE  hsa.hap_id = h.hap_id});

  create_and_load($self->{'dbVar'}, 'tmp_allele_group_allele','hap_id i*','hapset_id i*',
                  'name','snp_allele', 'subsnp_id i*');
  print $logh Progress::location();

  $self->{'dbVar'}->do(qq{ALTER TABLE allele_group ADD COLUMN hap_id int});
  print $logh Progress::location();

  debug(localtime() . "\tLoading allele_group");

  $self->{'dbVar'}->do(qq{INSERT INTO allele_group (variation_group_id, name, source_id, hap_id)
                SELECT vg.variation_group_id, tag.name, 1, tag.hap_id
                FROM   variation_group vg, tmp_allele_group_allele tag
                WHERE  vg.hapset_id = tag.hapset_id
                GROUP BY hap_id});
  print $logh Progress::location();

  $self->{'dbVar'}->do(qq{ALTER TABLE allele_group ADD INDEX hap_id(hap_id)});
  print $logh Progress::location();

  debug(localtime() . "\tLoading allele_group_allele");

  # there are a few haps in dbSNP which have two subsnps which have been
  # merged into the same refsnp.  Thus the group by clause.

  $self->{'dbVar'}->do(qq{INSERT INTO allele_group_allele (allele_group_id,variation_id, allele)
                SELECT ag.allele_group_id, vs.variation_id, taga.snp_allele
                FROM   allele_group ag, tmp_allele_group_allele taga,
                       variation_synonym vs
                WHERE  ag.hap_id = taga.hap_id
                AND    vs.subsnp_id = taga.subsnp_id
                GROUP BY ag.allele_group_id, vs.variation_id});
  print $logh Progress::location();

  $self->{'dbVar'}->do("DROP TABLE tmp_allele_group_allele");
  print $logh Progress::location();

  return;
}



#
# loads individual genotypes into the individual_genotype table,need variation_synonym and sample tables
#
sub individual_genotypes {
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
   @subind_tables = map{$_->[0]} @{$self->{'dbSNP'}->selectall_arrayref($stmt)};

   if ($self->{'dbCore'}->species !~ /homo|hum/i) {
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
  my $ind_gty_stmt = get_create_statement($self->{'dbVar'},'tmp_individual_genotype_single_bp',$self->{'master_schema_db'});
  print $logh Progress::location();

  #ÊDrop the tmp_individual.. table if it exists
  $stmt = qq{
    DROP TABLE IF EXISTS
      tmp_individual_genotype_single_bp
  };
  $self->{'dbVar'}->do($stmt);
  print $logh Progress::location();
  
   my $insert = 'INSERT';
   if ($self->{'dbCore'}->species !~ /homo|hum/i) {
    
    #ÊCreate the tmp_individual_genotype_single_bp table
    $self->{'dbVar'}->do($ind_gty_stmt);
    print $logh Progress::location();
    
     ###only non human one needs unique index
     $insert = "INSERT IGNORE";

     $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX ind_genotype_idx ON tmp_individual_genotype_single_bp(variation_id,subsnp_id,sample_id,allele_1,allele_2)});

     calculate_gtype($self,$self->{'dbVar'},"SubInd_tmp_gty","tmp_individual_genotype_single_bp",$insert);
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
	 calculate_gtype($self,$self->{'dbVar'},$table1,$table2,$insert);
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

  $dbVariation->do(qq{alter TABLE $table2 engine = MyISAM}) if ($self->{'dbCore'}->species =~ /homo|hum/i);
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

     create_and_load($self->{'dbVar'}, "tmp_pop_gty", 'subsnp_id i*', 'pop_id i*', 'freq',
 			'allele_1', 'allele_2');
  print $logh Progress::location();

   $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX pop_genotype_idx ON population_genotype(variation_id,subsnp_id,frequency,sample_id,allele_1,allele_2)});


    $self->{'dbVar'}->do(qq{INSERT IGNORE INTO population_genotype (variation_id,subsnp_id,allele_1, allele_2, frequency, sample_id)
			                  SELECT vs.variation_id,vs.subsnp_id,tra1.rev_allele as allele_1,tra2.rev_allele as allele_2,tg.freq,s.sample_id
					  FROM   variation_synonym vs, tmp_pop_gty tg,tmp_rev_allele tra1,tmp_rev_allele tra2, sample s
					  WHERE  vs.subsnp_id = tg.subsnp_id
                                          AND    tg.allele_1 = tra1.allele
                                          AND    tg.allele_2 = tra2.allele
                                          AND    vs.substrand_reversed_flag = 1
					  AND    s.pop_id = tg.pop_id});
  print $logh Progress::location();

    $self->{'dbVar'}->do(qq{INSERT IGNORE INTO population_genotype (variation_id,subsnp_id,allele_1, allele_2, frequency, sample_id)
					  SELECT vs.variation_id,vs.subsnp_id,tg.allele_1,tg.allele_2,tg.freq,s.sample_id
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

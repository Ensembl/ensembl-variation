use strict;
use warnings;

#generic object for the dbSNP data. Contains the general methods to dump the data into the new Variation database. Any change in the methods
# will need to overload the correspondent method in the subclass for the specie

package dbSNP::GenericContig;

use POSIX;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use ImportUtils qw(dumpSQL debug create_and_load load);

#creates the object and assign the attributes to it (connections, basically)
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbSNP, $dbCore, $dbVar, $snp_dbname, $tmp_dir, $tmp_file, $limit, $mapping_file_dir, $dbSNP_BUILD_VERSION, $ASSEMBLY_VERSION) =
        rearrange([qw(DBSNP DBCORE DBVAR SNP_DBNAME TMPDIR TMPFILE LIMIT MAPPING_FILE_DIR DBSNP_VERSION ASSEMBLY_VERSION)],@_);

  return bless {'dbSNP' => $dbSNP,
		'dbCore' => $dbCore,
		'dbVar' => $dbVar, ##this is a dbconnection
		'snp_dbname' => $snp_dbname,
		'tmpdir' => $tmp_dir,
		'tmpfile' => $tmp_file,
		'limit' => $limit,
                'mapping_file_dir' => $mapping_file_dir,
		'dbSNP_version' => $dbSNP_BUILD_VERSION,
		'assembly_version' => $ASSEMBLY_VERSION}, $class;
}

#main and only function in the object that dumps all dbSNP data 
sub dump_dbSNP{ 

  my $self = shift;

  $self->source_table();
  $self->population_table();
  $self->individual_table();
  $self->variation_table();
  $self->individual_genotypes();
  $self->population_genotypes();
  $self->allele_table();
  $self->flanking_sequence_table();
  $self->variation_feature();
  if ($self->{'dbCore'}->species =~ /homo/i) {
   $self->variation_group();
   $self->allele_group();
  }

  $self->cleanup();

}

sub source_table {
    my $self = shift;
    #my ($dbname,$version) = split /\_/,$self->{'dbSNP'}->dbname(); #get the version of the dbSNP release
    my ($dbname,$version) = split /\_/,$self->{'snp_dbname'};
    
    $self->{'dbVar'}->do(qq{INSERT INTO source (source_id,name,version) VALUES (1,"$dbname",$version)});

}


# filling of the variation table from SubSNP and SNP
# creating of a link table variation_id --> subsnp_id
sub variation_table {
    my $self = shift;

   $self->{'dbVar'}->do( "ALTER TABLE variation add column snp_id int" );

   debug("Dumping RefSNPs");
    
   my $count = $self->{'dbSNP'}->selectall_arrayref(qq{SELECT COUNT(*) FROM SNPAncestralAllele});
   ###if we have data in table SNPAncestralAllele, like human, will use it, otherwise goto else
   if ($count->[0][0]) {
     dumpSQL($self->{'dbSNP'},
	      qq{SELECT 1, concat( "rs", snp.snp_id), if(snp.validation_status = 0,NULL,
							 snp.validation_status), a.allele, snp.snp_id
		 FROM SNP snp
		 LEFT JOIN SNPAncestralAllele snpa ON snp.snp_id = snpa.snp_id
		 LEFT JOIN Allele a on snpa.ancestral_allele_id = a.allele_id
		 $self->{'limit'}
		}
	     );
   }
   else {
     dumpSQL($self->{'dbSNP'},  
	      qq{SELECT 1, concat( "rs", snp.snp_id), if(snp.validation_status = 0,NULL,
							 snp.validation_status), NULL, snp.snp_id
		 FROM SNP snp
		 $self->{'limit'}
		}
	     );
   }
    
    debug("Loading RefSNPs into variation table");
    
    load( $self->{'dbVar'}, "variation", "source_id", "name", "validation_status", "ancestral_allele", "snp_id" );
    
    $self->{'dbVar'}->do( "ALTER TABLE variation ADD INDEX snpidx( snp_id )" );
    
    # create a temp table of subSNP info
    
    debug("Dumping SubSNPs");
    
    $self->dump_subSNPs;

    create_and_load( $self->{'dbVar'}, "tmp_var_allele", "subsnp_id i*", "refsnp_id i*",
		     "pop_id i", "allele", "substrand_reversed_flag i", "moltype");
    
    # load the synonym table with the subsnp identifiers
    
   debug("loading variation_synonym table with subsnps");
   
   $self->{'dbVar'}->do(qq{ALTER TABLE variation_synonym add column subsnp_id int});
   $self->{'dbVar'}->do(qq{ALTER TABLE variation_synonym add column substrand_reversed_flag tinyint});
    
   $self->{'dbVar'}->do( qq{INSERT INTO variation_synonym (variation_id, source_id, name,
						  moltype, subsnp_id, substrand_reversed_flag )
		       SELECT v.variation_id, 1,
		       CONCAT('ss',tv.subsnp_id), tv.moltype, tv.subsnp_id,
		       tv.substrand_reversed_flag
		       FROM tmp_var_allele tv, variation v
		       WHERE tv.refsnp_id = v.snp_id
		       GROUP BY tv.subsnp_id
		   });
    
    $self->{'dbVar'}->do("ALTER TABLE variation_synonym ADD INDEX subsnp_id(subsnp_id)");
    
    ### FIX: Not sure if all RefSNPs have subsnps, and if ones which do not
    ### should possibly be eliminated

    return;
}


#
# dumps subSNPs and associated allele information
#
sub dump_subSNPs {
    my $self = shift;

    my $sth = $self->{'dbSNP'}->prepare
	(qq{SELECT subsnp.subsnp_id, subsnplink.snp_id, b.pop_id, ov.pattern,
	    subsnplink.substrand_reversed_flag, b.moltype
		FROM SubSNP subsnp, SNPSubSNPLink subsnplink, ObsVariation ov, Batch b
		WHERE subsnp.batch_id = b.batch_id
		AND   subsnp.subsnp_id = subsnplink.subsnp_id
		AND   ov.var_id = subsnp.variation_id
	    $self->{'limit'}}, {mysql_use_result => 1});

    
    $sth->execute();

    open ( FH, ">" . $self->{'tmpdir'} . "/" . $self->{'tmpfile'} );
    
  my ($row, $prefix);
  while($row = $sth->fetchrow_arrayref()) {
    my @alleles = split('/', $row->[3]);
    if ($row->[3] =~ /^(\(.*\))\d+\/\d+/) {
      $prefix = $1;
    }
    my @row = map {(defined($_)) ? $_ : '\N'} @$row;

    # split alleles into multiple rows
    foreach my $a (@alleles) {
      if ($prefix and $a !~ /\(.*\)/) {#incase (CA)12/13/14 CHANGE TO (CA)12/(CA)13/(CA)14
	$a = $prefix.$a;
	$prefix = "";
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
sub population_table {
    my $self = shift;
    
  $self->{'dbVar'}->do("ALTER TABLE sample ADD column pop_id int");   
  $self->{'dbVar'}->do("ALTER TABLE sample ADD column pop_class_id int"); 

  # load PopClassCode data as populations

  debug("Dumping population class data");

  dumpSQL($self->{'dbSNP'}, qq{SELECT pop_class, pop_class_id, pop_class_text
				   FROM PopClassCode 
			   });

  load($self->{'dbVar'}, 'sample', 'name', 'pop_class_id', 'description');

  $self->{'dbVar'}->do(qq{ALTER TABLE sample ADD INDEX pop_class_id (pop_class_id)});

  debug("Dumping population data");

  # load Population data as populations
    $self->{'dbSNP'}->do("SET SESSION group_concat_max_len = 10000");

    dumpSQL($self->{'dbSNP'}, qq{SELECT DISTINCT concat(p.handle, ':', p.loc_pop_id),
                    p.pop_id, pc.pop_class_id, GROUP_CONCAT(pl.line)
             FROM   Population p
             LEFT JOIN PopClass pc ON p.pop_id = pc.pop_id
	     LEFT JOIN PopLine pl ON p.pop_id = pl.pop_id
	     GROUP BY p.pop_id
	 });

    debug("Loading sample data");

    create_and_load( $self->{'dbVar'}, "tmp_pop", "name", "pop_id i*", "pop_class_id i*", "description l" );

    #populate the Sample table with the populations

    $self->{'dbVar'}->do(qq{INSERT INTO sample (name, pop_id,description)
                 SELECT tp.name, tp.pop_id, description
                 FROM   tmp_pop tp
                 GROUP BY tp.pop_id
                 });

    $self->{'dbVar'}->do(qq{ALTER TABLE sample ADD INDEX pop_id (pop_id)});

    #and copy the data from the sample to the Population table
    debug("Loading population table with data from Sample");

    $self->{'dbVar'}->do(qq{INSERT INTO population (sample_id)
				  SELECT sample_id
			          FROM sample});

     debug("Loading population_synonym table");

     # build super/sub population relationships
     $self->{'dbVar'}->do(qq{INSERT INTO population_structure (super_population_sample_id,sub_population_sample_id)
 				    SELECT p1.sample_id, p2.sample_id
 				    FROM tmp_pop tp, sample p1, sample p2
 				    WHERE tp.pop_class_id = p1.pop_class_id
 				    AND   tp.pop_id = p2.pop_id});
    

     #load population_synonym table with dbSNP population id
     $self->{'dbVar'}->do(qq{INSERT INTO sample_synonym (sample_id,source_id,name)
 				      SELECT sample_id, 1, pop_id
 				      FROM sample
 				      WHERE pop_id is NOT NULL
 				  });
    
     $self->{'dbVar'}->do("DROP TABLE tmp_pop");
}



# loads the individual table
#
sub individual_table {
    my $self = shift;

  # load individuals into the population table

  debug("Dumping Individual data");

  # a few submitted  individuals have the same individual or no individual
  # we ignore this problem with a group by
  #there were less individuals in the individual than in the individual_genotypes table, the reason is that some individuals do not have
  #assigned a specie for some reason in the individual table, but they do have in the SubmittedIndividual table
  #to solve the problem, get the specie information from the SubmittedIndividual table
  dumpSQL($self->{'dbSNP'}, qq{ SELECT IF(si.loc_ind_alias = '' ,si.loc_ind_id, si.loc_ind_alias), i.descrip, i.ind_id
				   FROM   SubmittedIndividual si, Individual i
				   WHERE  si.ind_id = i.ind_id
				   GROUP BY i.ind_id
			    });

  create_and_load($self->{'dbVar'}, 'tmp_ind', 'loc_ind_id', 'description', 'ind_id i*');

  # load pedigree into seperate tmp table because there are no
  # indexes on it in dbsnp and it makes the left join b/w tables v. slow
  # one individual has 2 (!) pedigree rows, thus the group by

  dumpSQL($self->{'dbSNP'}, qq{ SELECT ind_id, pa_ind_id, ma_ind_id, sex
              FROM PedigreeIndividual GROUP BY ind_id});

  create_and_load($self->{'dbVar'}, 'tmp_ped', 'ind_id i*', 'pa_ind_id i', 'ma_ind_id i', 'sex');

  debug("Loading individuals into individual table");

  # to make things easier keep dbSNPs individual.ind_id as our individual_id

  #add the individual_id column in the sample table
  $self->{'dbVar'}->do("ALTER TABLE sample ADD column individual_id int, add index ind_idx (individual_id)");   

  #and the individual data in the sample table
  $self->{'dbVar'}->do(qq{INSERT INTO sample (individual_id, name, description)
				  SELECT ti.ind_id, ti.loc_ind_id, ti.description
			          FROM tmp_ind ti
			      });

  #decide which individual_type should this species be
  my $individual_type_id;
  if ($self->{'dbCore'}->species =~ /homo|pan/i) {
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

    $self->{'dbVar'}->do("DROP table tmp_ind");

    #need to convert the father_individual_id and mother_individual_id in father_sample_individual_id and mother_sample_individual_id
    $self->{'dbVar'}->do("UPDATE individual i, sample s, individual i2 set i.father_individual_sample_id = s.sample_id where i.father_individual_sample_id = s.individual_id and i2.sample_id = s.sample_id");
   $self->{'dbVar'}->do("UPDATE individual i, sample s, individual i2 set i.mother_individual_sample_id = s.sample_id where i.mother_individual_sample_id = s.individual_id and i2.sample_id = s.sample_id");

    $self->{'dbVar'}->do("DROP table tmp_ped");
    
    #necessary to fill in the individual_population table with the relation between individual and populations
    dumpSQL($self->{'dbSNP'}, qq{ SELECT si.pop_id, i.ind_id
				      FROM   SubmittedIndividual si, Individual i
				      WHERE  si.ind_id = i.ind_id
			      });
    
    create_and_load($self->{'dbVar'}, 'tmp_ind_pop', 'pop_id i*', 'ind_id i*');

    debug("Loading individuals_population table");

    $self->{'dbVar'}->do(qq(INSERT INTO individual_population (individual_sample_id, population_sample_id)
				      SELECT s1.sample_id, s2.sample_id
				      FROM tmp_ind_pop tip, sample s1, sample s2
				      WHERE tip.pop_id = s2.pop_id
				      AND s1.individual_id = tip.ind_id
				  ));
    
    $self->{'dbVar'}->do("DROP table tmp_ind_pop");

    #necessary to fill in the sample_synonym table with the relation between individual_id and sample_id
    $self->{'dbVar'}->do(qq{INSERT INTO sample_synonym (sample_id,source_id,name)
 				      SELECT sample_id, 1, individual_id
 				      FROM sample
 				      WHERE individual_id is NOT NULL
 				  });

    return;
}


#
# loads the allele table
#
sub allele_table {
    my $self = shift;

    my $allele_table_ref = $self->{'dbVar'}->db_handle->selectall_arrayref(qq{show tables like "tmp_rev_allele"});
    my $allele_table = $allele_table_ref->[0][0];
    if (! $allele_table) {
      debug("Dumping allele data");

      dumpSQL($self->{'dbSNP'}, qq(SELECT a1.allele_id, a1.allele, a2.allele
                                   FROM Allele a1, Allele a2
                                   WHERE a1.rev_allele_id = a2.allele_id));
    
      create_and_load($self->{'dbVar'}, "tmp_rev_allele", "allele_id i*","allele *", "rev_allele");
    }
    #first load the allele data for alleles that we have population and
    #frequency data for

    dumpSQL($self->{'dbSNP'}, qq(SELECT afsp.subsnp_id, afsp.pop_id, afsp.allele_id, afsp.freq
			       FROM   AlleleFreqBySsPop afsp, SubSNP ss
			       WHERE    afsp.subsnp_id = ss.subsnp_id
			       $self->{'limit'}));

   debug("Loading allele frequency data");

   create_and_load($self->{'dbVar'}, "tmp_allele", "subsnp_id i*", "pop_id i*",
		    "allele_id i*", "freq");

    debug("Creating allele table");
  
  #necessary to create a unique index to simulate the GROUP BY clause
  $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX unique_allele_idx ON allele (variation_id,allele(2),frequency,sample_id)});

   $self->{'dbVar'}->do(qq(INSERT IGNORE INTO allele (variation_id, allele,frequency, sample_id)
               SELECT vs.variation_id,
                      IF(vs.substrand_reversed_flag,
                         tra.rev_allele,tra.allele) as allele,
                      ta.freq, s.sample_id
               FROM   tmp_allele ta, tmp_rev_allele tra, variation_synonym vs,
                      sample s
               WHERE  ta.subsnp_id = vs.subsnp_id
	       AND    ta.allele_id = tra.allele_id
	       AND    ta.pop_id = s.pop_id),{mysql_use_result => 1} );

   $self->{'dbVar'}->do("ALTER TABLE allele ENABLE KEYS"); #after ignoring in the insertion, we must enable keys again
   
   $self->{'dbVar'}->do("DROP TABLE tmp_allele");    
   #going to add the other allele for the variations with 1 allele (have frequency 1 but no frequency for the other allele)
   debug("Loading allele data without frequency");
   
   $self->{'dbVar'}->do("CREATE TABLE tmp_allele (variation_id int, allele text, primary key (variation_id,allele(10)))");
   $self->{'dbVar'}->do("INSERT IGNORE INTO tmp_allele SELECT variation_id, allele FROM allele");
   $self->{'dbVar'}->do(qq{CREATE TABLE tmp_unique_allele 
				      SELECT ta.variation_id,  ta.allele, vs.snp_id
				         FROM variation vs,
				                    (SELECT variation_id, allele 
						     FROM tmp_allele
						     GROUP BY variation_id 
						     HAVING COUNT(*) = 1) as ta
					 WHERE ta.variation_id = vs.variation_id});

   $self->{'dbVar'}->do("CREATE INDEX tmp_unique_allele_idx on tmp_unique_allele (variation_id)");
   $self->{'dbVar'}->do("DROP TABLE tmp_allele");
   #create table with unique alleles from dbSNP very slow with unique index??????
   $self->{'dbVar'}->do("CREATE TABLE tmp_allele (refsnp_id int, allele text, primary key (refsnp_id,allele(10)))");

   $self->{'dbVar'}->do(qq{INSERT IGNORE INTO tmp_allele
				      SELECT tva.refsnp_id, IF (tva.substrand_reversed_flag, tra.rev_allele,tva.allele) as allele
				      FROM tmp_var_allele tva, tmp_rev_allele tra, tmp_unique_allele tua
				      WHERE tva.allele = tra.allele
				      AND tua.snp_id = tva.refsnp_id
				  });

   $self->{'dbVar'}->do(qq{INSERT IGNORE INTO allele (variation_id, allele, frequency)
						  SELECT tua.variation_id, ta.allele,0
						  FROM tmp_unique_allele tua, tmp_allele ta
						  WHERE tua.snp_id = ta.refsnp_id
						  AND tua.allele <> ta.allele
						  });
   #remove the index
   $self->{'dbVar'}->do("DROP INDEX unique_allele_idx ON allele");
   $self->{'dbVar'}->do("DROP TABLE tmp_unique_allele");
   $self->{'dbVar'}->do("DROP TABLE tmp_allele");

   # load remaining allele data which we do not have frequency data for
   # this will not import alleles without frequency for a variation which
   # already has frequency
   
   debug("Loading other allele data");

   $self->{'dbVar'}->do(qq{CREATE TABLE tmp_allele
		   SELECT vs.variation_id as variation_id, tva.pop_id,
                         IF(vs.substrand_reversed_flag,
                            tra.rev_allele, tra.allele) as allele
		   FROM   variation_synonym vs, tmp_var_allele tva,
                         tmp_rev_allele tra
		   WHERE  tva.subsnp_id = vs.subsnp_id
                   AND    tva.allele = tra.allele 
                   AND    NOT EXISTS ## excluding alleles that already in allele table
                       (SELECT * FROM allele a where a.variation_id = vs.variation_id
                        AND a.allele = IF(tva.substrand_reversed_flag,tra.rev_allele, tra.allele))});


   $self->{'dbVar'}->do("ALTER TABLE tmp_allele ADD INDEX pop_id(pop_id)");

   $self->{'dbVar'}->do(qq{INSERT INTO allele (variation_id, allele,
                                      frequency, sample_id)
                  SELECT ta.variation_id, ta.allele, null, s.sample_id
                  FROM   tmp_allele ta
                  LEFT JOIN sample s ON s.pop_id = ta.pop_id
                  GROUP BY ta.variation_id, s.sample_id, ta.allele });

    ###generate allele_string table for variation_feature allele_string column

    dumpSQL($self->{'dbSNP'}, qq(select snp.snp_id, uv.var_str from SNP snp, UniVariation uv 
                                 where snp.univar_id = uv.univar_id));

    open ( IN, $self->{'tmpdir'} . "/" . $self->{'tmpfile'} );
    open OUT, ">" . $self->{'tmpdir'} . "/allele_string.txt";
    my $prefix; ##in case like (AT)/10/11/12->(AT)10/(AT)11/..because this is how it appear in tmp_rev_allele table
    while (<IN>) {
      my ($snp_id,$var_str) = split /\t/,$_;
      my @alleles = split /\//, $var_str;
      if ($var_str =~ /^(\(.*\))\d+\/\d+/) {
	$prefix = $1;
      }
      foreach my $allele (@alleles) {
	if ($prefix and $allele !~ /\(.*\)/) {
	  $allele = $prefix.$allele ;
	}
	print OUT "rs$snp_id\t$allele\n";
      }
    }
    my $call = "mv ". $self->{'tmpdir'}."/allele_string.txt ". $self->{'tmpdir'} . "/" . $self->{'tmpfile'};
    print "$call\n";
    system($call);
    create_and_load($self->{'dbVar'},"tmp_allele_string","snp_name *","allele");
    
    $self->{'dbVar'}->do(qq{CREATE TABLE allele_string select v.variation_id as variation_id, tas.allele as allele                                  FROM variation v, tmp_allele_string tas
                                  WHERE v.name = tas.snp_name});

    $self->{'dbVar'}->do(qq{ALTER TABLE allele_string ADD INDEX variation_idx(variation_id)});

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

  $self->{'dbVar'}->do(qq{CREATE TABLE tmp_seq (variation_id int,
						      subsnp_id int,
						      line_num int,
						      type enum ('5','3'),
						      line varchar(255),
						      revcom tinyint)
				MAX_ROWS = 100000000});

  # import both the 5prime and 3prime flanking sequence tables

  
  foreach my $type ('3','5') {
    debug("Dumping $type' flanking sequence");
    
    dumpSQL($self->{'dbSNP'}, qq{SELECT seq.subsnp_id, seq.line_num, seq.line
				 FROM SubSNPSeq$type seq, SNP snp
				 WHERE snp.exemplar_subsnp_id = seq.subsnp_id
				 $self->{'limit'}});
    

    $self->{'dbVar'}->do(qq{CREATE TABLE tmp_seq_$type (
							      subsnp_id int,
							      line_num int,
							      line varchar(255),
							      KEY subsnp_id_idx(subsnp_id))
				  MAX_ROWS = 100000000 });
    
    load($self->{'dbVar'}, "tmp_seq_$type", "subsnp_id", "line_num", "line");

    # merge the tables into a single tmp table
    $self->{'dbVar'}->do(qq{INSERT INTO tmp_seq (variation_id, subsnp_id,
						       line_num, type, line, revcom)
				  SELECT vs.variation_id, ts.subsnp_id, ts.line_num, '$type',
				  ts.line, vs.substrand_reversed_flag
				  FROM   tmp_seq_$type ts, variation_synonym vs
				  WHERE  vs.subsnp_id = ts.subsnp_id});
    #drop tmp table to free space
    $self->{'dbVar'}->do(qq{DROP TABLE tmp_seq_$type});
  }
      
  $self->{'dbVar'}->do("ALTER TABLE tmp_seq ADD INDEX idx (subsnp_id, type, line_num)");

  my $sth = $self->{'dbVar'}->prepare(qq{SELECT ts.variation_id, ts.subsnp_id, ts.type,
					       ts.line, ts.revcom
					       FROM   tmp_seq ts FORCE INDEX (idx)
					       ORDER BY ts.subsnp_id, ts.type, ts.line_num},{mysql_use_result => 1});
  
  $sth->execute();

  my ($vid, $ssid, $type, $line, $revcom);

  $sth->bind_columns(\$vid, \$ssid, \$type, \$line, \$revcom);

  open(FH, ">" . $self->{'tmpdir'} . "/" . $self->{'tmpfile'});
  my $upstream = '';
  my $dnstream = '';
  my $cur_vid;
  my $cur_revcom;

  debug("Rearranging flanking sequence data");


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

  debug("Loading flanking sequence data");

  # import the generated data
  load($self->{'dbVar'},"flanking_sequence","variation_id","up_seq","down_seq");

  unlink($self->{'tmpdir'} . "/" . $self->{'tmpfile'});

  return;
}



sub variation_feature {
  my $self = shift;

  debug("Dumping seq_region data");

  dumpSQL($self->{'dbCore'}->dbc()->db_handle, qq{SELECT sr.seq_region_id, sr.name
 				      FROM   seq_region sr, coord_system cs
                                      WHERE  sr.coord_system_id = cs.coord_system_id
                                      AND    cs.attrib like "%default_version%"});
    
  debug("Loading seq_region data");
  create_and_load($self->{'dbVar'}, "tmp_seq_region", "seq_region_id", "name *");

  debug("Dumping SNPLoc data");

  my ($tablename1,$tablename2,$row);

  my ($assembly_version) =  $self->{'assembly_version'} =~ /^[a-zA-Z]+(\d+)\.*.*$/;
  print "assembly_version again is $assembly_version\n";

  my $sth = $self->{'dbSNP'}->prepare(qq{SHOW TABLES LIKE 
					 '$self->{'dbSNP_version'}\_SNPContigLoc\_$assembly_version\__'});
  $sth->execute();

  while($row = $sth->fetchrow_arrayref()) {
    $tablename1 = $row->[0];
  }

  my $sth1 = $self->{'dbSNP'}->prepare(qq{SHOW TABLES LIKE 
					   '$self->{'dbSNP_version'}\_ContigInfo\_$assembly_version\__'});
  $sth1->execute();

  while($row = $sth1->fetchrow_arrayref()) {
    $tablename2 = $row->[0];
  }
  debug("table_name1 is $tablename1 table_name2 is $tablename2");
  #note the contig based cordinate is 0 based, ie. start at 0, lc_ngbr+2, t1.rc_ngbr
      dumpSQL($self->{'dbSNP'}, qq{SELECT t1.snp_id, t2.contig_acc,
 				  t1.lc_ngbr+2,t1.rc_ngbr,
   				  IF(t1.orientation, -1, 1)
 				  FROM   $tablename1 t1, $tablename2 t2
 				  WHERE t1.ctg_id=t2.ctg_id
 				  $self->{'limit'}});
    
    
  debug("Loading SNPLoc data");

  create_and_load($self->{'dbVar'}, "tmp_contig_loc", "snp_id i*", "contig *", "start i", 
		  "end i", "strand i");
    
  debug("Creating genotyped variations");
    
  #creating the temporary table with the genotyped variations
  #$self->{'dbVar'}->do(qq{CREATE TABLE tmp_genotyped_var SELECT DISTINCT variation_id FROM tmp_individual_genotype_single_bp});
  #$self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX variation_idx ON tmp_genotyped_var (variation_id)});
  #$self->{'dbVar'}->do(qq{INSERT IGNORE INTO  tmp_genotyped_var SELECT DISTINCT variation_id FROM individual_genotype_multiple_bp});

    
  debug("Creating tmp_variation_feature data");
    
  dumpSQL($self->{'dbVar'},qq{SELECT v.variation_id, ts.seq_region_id, tcl.start, tcl.end, tcl.strand, v.name, v.source_id, v.validation_status
					  FROM variation v, tmp_contig_loc tcl, tmp_seq_region ts
					  WHERE v.snp_id = tcl.snp_id
					  AND ts.name = tcl.contig});
    
  create_and_load($self->{'dbVar'},'tmp_variation_feature',"variation_id *","seq_region_id i", "seq_region_start i", "seq_region_end i", "seq_region_strand i", "variation_name", "source_id i", "validation_status i");
    
  debug("Dumping data into variation_feature table");
  $self->{'dbVar'}->do(qq{INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,variation_name, flags, source_id, validation_status)
				SELECT tvf.variation_id, tvf.seq_region_id, tvf.seq_region_start, tvf.seq_region_end, tvf.seq_region_strand,tvf.variation_name,IF(tgv.variation_id,'genotyped',NULL), tvf.source_id, tvf.validation_status
				FROM tmp_variation_feature tvf LEFT JOIN tmp_genotyped_var tgv ON tvf.variation_id = tgv.variation_id
				 });
    
    #$self->{'dbVar'}->do("DROP TABLE tmp_contig_loc");
    #$self->{'dbVar'}->do("DROP TABLE tmp_seq_region");
    #$self->{'dbVar'}->do("DROP TABLE tmp_genotyped_var");
    #$self->{'dbVar'}->do("DROP TABLE tmp_variation_feature");
}

#
# loads variation_group and variation_group_variation tables from the
# contents of the HapSet and HapSetSnpList tables
#
sub variation_group {
    my $self = shift;

  debug("Dumping HapSet data");

  dumpSQL($self->{'dbSNP'}, qq{SELECT  CONCAT(hs.handle, ':', hs.hapset_name),
                    hs.hapset_id, hssl.subsnp_id
             FROM HapSet hs, HapSetSnpList hssl   ###, SubSNP ss 
             WHERE hs.hapset_id = hssl.hapset_id});

  create_and_load($self->{'dbVar'}, 'tmp_var_grp', 'name', 'hapset_id i*', 'subsnp_id i*');

  $self->{'dbVar'}->do("ALTER TABLE variation_group add column hapset_id int");

  debug("Loading variation_group");

  $self->{'dbVar'}->do(qq{INSERT INTO variation_group (name, source_id, type, hapset_id)
                SELECT name, 1, 'haplotype', hapset_id
                FROM tmp_var_grp
                GROUP BY hapset_id});

  $self->{'dbVar'}->do("ALTER TABLE variation_group ADD INDEX hapset_id(hapset_id)");

  debug("Loading variation_group_variation");

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

  $self->{'dbVar'}->do("DROP TABLE tmp_var_grp");
}

#
# loads allele_group table
#
sub allele_group {
    my $self = shift;

  debug("Dumping Hap data");

  dumpSQL($self->{'dbSNP'}, qq{SELECT  h.hap_id, h.hapset_id, h.loc_hap_id,
                    hsa.snp_allele, hsa.subsnp_id
             FROM   Hap h, HapSnpAllele hsa   ###, SubSNP ss
             WHERE  hsa.hap_id = h.hap_id});

  create_and_load($self->{'dbVar'}, 'tmp_allele_group_allele','hap_id i*','hapset_id i*',
                  'name','snp_allele', 'subsnp_id i*');

  $self->{'dbVar'}->do(qq{ALTER TABLE allele_group ADD COLUMN hap_id int});

  debug("Loading allele_group");

  $self->{'dbVar'}->do(qq{INSERT INTO allele_group (variation_group_id, name, source_id, hap_id)
                SELECT vg.variation_group_id, tag.name, 1, tag.hap_id
                FROM   variation_group vg, tmp_allele_group_allele tag
                WHERE  vg.hapset_id = tag.hapset_id
                GROUP BY hap_id});

  $self->{'dbVar'}->do(qq{ALTER TABLE allele_group ADD INDEX hap_id(hap_id)});

  debug("Loading allele_group_allele");

  # there are a few haps in dbSNP which have two subsnps which have been
  # merged into the same refsnp.  Thus the group by clause.

  $self->{'dbVar'}->do(qq{INSERT INTO allele_group_allele (allele_group_id,variation_id, allele)
                SELECT ag.allele_group_id, vs.variation_id, taga.snp_allele
                FROM   allele_group ag, tmp_allele_group_allele taga,
                       variation_synonym vs
                WHERE  ag.hap_id = taga.hap_id
                AND    vs.subsnp_id = taga.subsnp_id
                GROUP BY ag.allele_group_id, vs.variation_id});

  $self->{'dbVar'}->do("DROP TABLE tmp_allele_group_allele");

  return;
}



#
# loads individual genotypes into the individual_genotype table,need variation_synonym and sample tables
#
sub individual_genotypes {
   my $self = shift;

   #
   # load SubInd individual genotypes into genotype table
   #

   my $TMP_DIR = $self->{'tmpdir'};
   my $TMP_FILE = $self->{'tmpfile'};

    debug("Dumping SubInd and ObsGenotype data");

    debug("Time starting to dump tmp2_gty table: ",scalar(localtime(time)));
    dumpSQL($self->{'dbSNP'}, qq{select si.subsnp_id, sind.ind_id, length(og.obs) as length_pat,
				 IF (length(og.obs)=3, SUBSTRING(og.obs,1,1),SUBSTRING_INDEX(og.obs,'/',1)) as allele_1,
				 IF (length(og.obs)=3, SUBSTRING(og.obs,3,1),SUBSTRING_INDEX(og.obs,'/',-1)) as allele_2,
				 si.submitted_strand_code
				 FROM   SubInd si, ObsGenotype og, SubmittedIndividual sind
				 WHERE  og.gty_id = si.gty_id 
				 AND    sind.submitted_ind_id = si.submitted_ind_id
                                 AND    og.obs != 'N/N'
				 $self->{'limit'}});

    create_and_load($self->{'dbVar'}, "tmp2_gty", 'subsnp_id i*', 'ind_id i', 'length_pat i','allele_1', 'allele_2','submitted_strand i');

    debug("Time finishing to dump tmp2_gty table: ",scalar(localtime(time)));
    debug("Loading individual_genotype table");

    # load a temp table that can be used to reverse compliment alleles
    # we place subsnps in the same orientation as the refSNP

   my $allele_table_ref = $self->{'dbVar'}->db_handle->selectall_arrayref(qq{show tables like "tmp_rev_allele"});
   my $allele_table = $allele_table_ref->[0][0];
   if (! $allele_table) {
     debug("Dumping allele data");
     dumpSQL($self->{'dbSNP'}, qq(SELECT a1.allele_id, a1.allele, a2.allele
                                 FROM Allele a1, Allele a2
                                 WHERE a1.rev_allele_id = a2.allele_id));

     create_and_load($self->{'dbVar'}, "tmp_rev_allele", "allele_id i*","allele *", "rev_allele");
   }

    #we have truncated the individual_genotype table, one contains the genotypes single bp, and the other, the rest
    #necessary to create a unique index to remove duplicated genotypes

    $self->{'dbVar'}->do(qq{CREATE INDEX individual_id on sample(individual_id)});

    $self->{'dbVar'}->do(qq{CREATE TABLE tmp_individual_genotype_single_bp (
                           variation_id int not null,allele_1 varchar(255),allele_2 varchar(255),sample_id int,
                           key variation_idx(variation_id),
                           key sample_idx(sample_id)
                           ) MAX_ROWS = 100000000
    });

    if ($self->{'dbCore'}->species !~ /homo|hum|mouse|mus/i) {
      $self->{'dbVar'}->do(qq{CREATE UNIQUE INDEX ind_genotype_idx ON tmp_individual_genotype_single_bp(variation_id,sample_id,allele_1,allele_2)});
      my $insert = 'INSERT IGNORE';
      calculate_gtype($self,$self->{'dbVar'},"tmp2_gty","tmp_individual_genotype_single_bp",$insert);
    }
   else {
     my $total_row_ref = $self->{'dbVar'}->selectall_arrayref(qq{select count(*) from tmp2_gty});
     my $total_row = $total_row_ref->[0][0];
     my $limit_start = 0;
     my $num_tables = 4;
     my $insert = 'INSERT';
     my %rec_pid;
     my $table_row = int($total_row/$num_tables);
     #my $table_row = int(10000/$num_tables); for debugging
     debug("splitting tmp2_gty table into $num_tables tables");
     foreach my $num (1..$num_tables+1) {
       my $table1 = "tmp2_gty\_$num";
       my $table2 = "tmp_individual_genotype_single_bp\_$num";
       debug("table1 is $table1 and table2 is $table2");
       $self->{'dbVar'}->do(qq{CREATE TABLE $table1 like tmp2_gty});
       $self->{'dbVar'}->do(qq{CREATE TABLE $table2 like tmp_individual_genotype_single_bp});
       $self->{'dbVar'}->do(qq{ALTER TABLE $table2 engine = innodb});
       $self->{'dbVar'}->do(qq{DROP INDEX ind_genotype_idx ON $table2});
       $self->{'dbVar'}->do(qq{INSERT INTO $table1 select * from tmp2_gty limit $limit_start, $table_row});
       
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
       $limit_start += $table_row;
       debug("next limit_start = $limit_start");
       $rec_pid{$pid}++;
     }
     #the parent waits for the child to finish;
     foreach my $p (keys %rec_pid) {
       waitpid($p,0);
     }
     #system("cat $TMP_DIR/tmp_individual_genotype_single_bp_*\.su >$TMP_DIR/$TMP_FILE");
     #load ($self->{'dbVar'},"tmp_individual_genotype_single_bp","variation_id","allele_1", "allele_2","sample_id");
     #$self->{'dbVar'}->do(qq{DROP INDEX ind_genotype_idx ON individual_genotype_single_bp});
     debug("Time finishing to insert to tmp_individual_genotype_single_bp table: ",scalar(localtime(time)));

     debug("Time start to insert to tmp_individual_genotype_multiple_bp table");
     $self->{'dbVar'}->do(qq{INSERT INTO individual_genotype_multiple_bp (variation_id, sample_id, allele_1, allele_2) 
				  SELECT vs.variation_id, s.sample_id,
				  SUBSTRING_INDEX(tg.allele_1,'/',1) as allele_1,
				  SUBSTRING_INDEX(tg.allele_2,'/',-1) as allele_2 
				  FROM   tmp2_gty tg, variation_synonym vs, sample s
				  WHERE  tg.subsnp_id = vs.subsnp_id
				  AND    tg.length_pat > 3
				  AND    tg.ind_id = s.individual_id});

     #$self->{'dbVar'}->do("DROP TABLE tmp2_gty");
     debug("Time finishing to insert to tmp_individual_genotype_multiple_bp table: ",scalar(localtime(time)));
     return;
   }
}

sub calculate_gtype {
  
  my ($self,$dbVariation,$table1,$table2,$insert) = @_;


  debug("Time starting to insert into $table2 : ",scalar(localtime(time)));
  $dbVariation->do(qq{$insert INTO $table2 (variation_id, sample_id, allele_1, allele_2) 
				 SELECT vs.variation_id, s.sample_id, tra1.allele, tra2.allele
				 FROM $table1 tg, variation_synonym vs, tmp_rev_allele tra1, tmp_rev_allele tra2, sample s
				 WHERE tg.submitted_strand IN (1,3,4)
				 AND vs.substrand_reversed_flag =1
				 AND tg.subsnp_id = vs.subsnp_id
				 AND tra1.allele = tg.allele_1
				 AND tra2.allele = tg.allele_2
				 AND tg.length_pat = 3
				 AND tg.ind_id = s.individual_id});
  debug("Time starting to insert into $table2 table 2: ",scalar(localtime(time)));

  $dbVariation->do(qq{$insert INTO $table2 (variation_id, sample_id, allele_1, allele_2)
				 SELECT vs.variation_id, s.sample_id, tra1.rev_allele, tra2.rev_allele
				 FROM $table1 tg, variation_synonym vs, tmp_rev_allele tra1, tmp_rev_allele tra2, sample s
				 WHERE tg.submitted_strand IN (1,3,4)
				 AND vs.substrand_reversed_flag =0
				 AND tg.subsnp_id = vs.subsnp_id
				 AND tra1.allele = tg.allele_1
				 AND tra2.allele = tg.allele_2
				 AND tg.length_pat = 3
				 AND tg.ind_id = s.individual_id});

  debug("Time starting to insert into $table2 table 3: ",scalar(localtime(time)));

  $dbVariation->do(qq{$insert INTO $table2 (variation_id, sample_id, allele_1, allele_2)

				 SELECT vs.variation_id, s.sample_id, tra1.rev_allele, tra2.rev_allele
				 FROM $table1 tg, variation_synonym vs, tmp_rev_allele tra1, tmp_rev_allele tra2, sample s
				 WHERE tg.submitted_strand IN (0,2,5)
				 AND vs.substrand_reversed_flag =1
				 AND tg.subsnp_id = vs.subsnp_id
				 AND tra1.allele = tg.allele_1
				 AND tra2.allele = tg.allele_2
				 AND tg.length_pat = 3
				 AND tg.ind_id = s.individual_id});

  debug("Time starting to insert into $table2 table 4: ",scalar(localtime(time)));

  $dbVariation->do(qq{$insert INTO $table2 (variation_id, sample_id, allele_1, allele_2)

				 SELECT vs.variation_id, s.sample_id, tra1.allele, tra2.allele
				 FROM $table1 tg, variation_synonym vs, tmp_rev_allele tra1, tmp_rev_allele tra2, sample s
				 WHERE tg.submitted_strand IN (0,2,5)
				 AND vs.substrand_reversed_flag =0
				 AND tg.subsnp_id = vs.subsnp_id
				 AND tra1.allele = tg.allele_1
				 AND tra2.allele = tg.allele_2
				 AND tg.length_pat = 3
				 AND tg.ind_id = s.individual_id});

  if ($table2 ne "tmp_individual_genotype_single_bp") {
    
    my $tmp_dir = $self->{'tmpdir'};
    my $user = $dbVariation->username;
    my $pass = $dbVariation->password;
    my $host = $dbVariation->host;
    my $dbname = $dbVariation->dbname;
    system("chmod 777 $tmp_dir");
    my $call = "mysqldump -T $tmp_dir -u $user -p$pass -h $host $dbname $table2";
    system($call);
    system("sort $tmp_dir/$table2\.txt |uniq >$tmp_dir/$table2\.su");
  }
}

#     if ($self->{'dbCore'}->species =~ /mouse|mus/i) {
#       $self->{'dbVar'}->do(qq{CREATE TABLE tmp_strain_gtype_final (
# 								      variation_id int not null, 
# 								      allele_string char(3) not null,
# 								      sample_id1 int(11) not null,
# 								      sample_id2 int(11) not null,
#                                                                       sample_name1 varchar(50) not null,
#                                                                       sample_name2 varchar(50) not null,
# 								      unique key variation_idx(variation_id,sample_id1,sample_id2),
# 								      key sample_idx(sample_id1,sample_id2)) MAX_ROWS = 100000000
# 				});
#       $self->{'dbVar'}->do(qq{INSERT IGNORE into tmp_strain_gtype_final
# 				  SELECT ig1.variation_id,concat(ig1.allele_1,"/",ig2.allele_1) as allele_string,ip1.population_sample_id as sample_id1,ip2.population_sample_id as sample_id2, s1.name as sample_name1, s2.name as sample_name2
# 				  FROM population p1, population p2,tmp_individual_genotype_single_bp ig1, tmp_individual_genotype_single_bp ig2,individual_population ip1,individual_population ip2, sample s1, sample s2
# 				  WHERE p1.is_strain=1 
# 				  AND p2.is_strain=1 
# 				  AND ig1.allele_1=ig1.allele_2 
# 				  AND ig2.allele_1 = ig2.allele_2 
# 				  AND ig1.variation_id=ig2.variation_id 
# 				  AND ig1.sample_id=ip1.individual_sample_id 
# 				  AND ig2.sample_id=ip2.individual_sample_id 
# 				  AND p1.sample_id = ip1.population_sample_id 
# 				  AND p2.sample_id = ip2.population_sample_id 
# 				  AND ig1.allele_1 != ig2.allele_1 
# 				  AND ig1.allele_1 != 'N' 
# 				  AND ig2.allele_1 != 'N' 
# 				  AND p2.sample_id != p1.sample_id
#                                   AND s1.sample_id = p1.sample_id
#                                   AND s2.sample_id = p2.sample_id
# 				 });
#Query OK, 141999022 rows affected (13 hours 54 min 5.95 sec)
#Records: 142921258  Duplicates: 922236  Warnings: 0
#       $self->{'dbVar'}->do(qq{CREATE TABLE tmp_strain_gtype_count_final 
# 				  SELECT sg1.sample_id1 as sample_id1,sg2.sample_id2 as sample_id2,count(*) as count 
# 				  FROM tmp_strain_gtype_final sg1 
#                                        LEFT JOIN tmp_ids1 t1 on sg1.variation_id = t1.variation_id1,  
#                                      tmp_strain_gtype_final sg2
# 				  WHERE sg1.variation_id=sg2.variation_id 
# 				  AND sg1.sample_id1 != sg2.sample_id2 
# 				  AND sg1.sample_id1 = sg2.sample_id1 
# 				  AND sg1.sample_id2 = sg2.sample_id2 
#                                 AND t1.variation_id1 is NULL #to avoid duplicate variation from both dbsnp and sanger
# 				  GROUP BY sg1.sample_id1,sg2.sample_id2});
#     $self->{'dbVar'}->do("ALTER TABLE tmp_strain_gtype_count_final ADD INDEX sample_idx(sample_id1,sample_id2)");
#Query OK, 6578 rows affected (29 min 33.50 sec)
#Records: 6578  Duplicates: 0  Warnings: 0
#     #$self->{'dbVariation'}->do(qq{CREATE TABLE strain_gtype_final ( variation_id int not null, allele_string char(3) not null, sample_id1 int(11) not null, sample_id2 int(11) not null, sample_name1 varchar(50) not null, sample_name2 varchar(50) not null, count int, key sample_idx(sample_id1,sample_id2))});
#     $self->{'dbVariation'}->do(qq{CREATE TABLE strain_gtype_final select s.*,c.count from tmp_strain_gtype_final s, tmp_strain_gtype_count_final c where s.sample_id1 = c.sample_id1 and s.sample_id2 = c.sample_id2 and c.count < 110000});
#     $self->{'dbVariation'}->do("ALTER TABLE strain_gtype_final ADD INDEX sample_id1(sample_id1), add INDEX sample_id2(sample_id2), add index variation_id(variation_id)");
#     create tmp_sample table contail only 46 selected sample
#     create table tmp_sample select s.sample_id,s.name from sample s, meta m, population p where s.name=m.meta_value and s.sample_id=p.sample_id and p.is_strain=1;
# only use tmp_sample instead of sample and population tables
#INSERT IGNORE into tmp_strain_gtype_final SELECT ig1.variation_id,concat(ig1.allele_1,"/",ig2.allele_1) as allele_string,ip1.population_sample_id as sample_id1,ip2.population_sample_id as sample_id2, s1.name as sample_name1, s2.name as sample_name2 FROM tmp_individual_genotype_single_bp ig1, tmp_individual_genotype_single_bp ig2,individual_population ip1,individual_population ip2, tmp_sample s1, tmp_sample s2 WHERE ig1.allele_1=ig1.allele_2 AND ig2.allele_1 = ig2.allele_2 AND ig1.variation_id=ig2.variation_id AND ig1.sample_id=ip1.individual_sample_id AND ig2.sample_id=ip2.individual_sample_id AND s1.sample_id = ip1.population_sample_id AND s2.sample_id = ip2.population_sample_id AND ig1.allele_1 != ig2.allele_1 AND s2.sample_id != s1.sample_id;
#updated see strain_table.sql in import directory
#     } ###end of if ($self->{'dbCore'}->species =~ /mouse|mus/i)



#
# loads population genotypes into the population_genotype table
#
sub population_genotypes {
    my $self = shift;

    my $allele_table_ref = $self->{'dbVar'}->db_handle->selectall_arrayref(qq{show tables like "tmp_rev_allele"});
    my $allele_table = $allele_table_ref->[0][0];
    if (! $allele_table) {
      debug("Dumping allele data");
      
      dumpSQL($self->{'dbSNP'}, qq(SELECT a1.allele_id, a1.allele, a2.allele
                                   FROM Allele a1, Allele a2
                                   WHERE a1.rev_allele_id = a2.allele_id));
    
      create_and_load($self->{'dbVar'}, "tmp_rev_allele", "allele_id i*","allele *", "rev_allele");
    }
    debug("Dumping GtyFreqBySsPop and UniGty data");
    
    dumpSQL($self->{'dbSNP'},qq{SELECT gtfsp.subsnp_id, gtfsp.pop_id, gtfsp.freq,a1.allele, a2.allele
					FROM   GtyFreqBySsPop gtfsp, UniGty ug, Allele a1, Allele a2
					WHERE  gtfsp.unigty_id = ug.unigty_id
					AND    ug.allele_id_1 = a1.allele_id
					AND    ug.allele_id_2 = a2.allele_id
					$self->{'limit'}
				    });
	
    debug("loading genotype data");
	
    create_and_load($self->{'dbVar'}, "tmp_gty", 'subsnp_id i*', 'pop_id i*', 'freq',
			'allele_1', 'allele_2');

    $self->{'dbVar'}->do(qq{INSERT INTO population_genotype (variation_id,allele_1, allele_2, frequency, sample_id)
					  SELECT vs.variation_id,IF(vs.substrand_reversed_flag,tra.rev_allele,tra.allele) as allele_1,IF(vs.substrand_reversed_flag,tra.rev_allele,tra.allele) as allele_2,tg.freq,s.sample_id
					  FROM   variation_synonym vs, tmp_gty tg, tmp_rev_allele tra, sample s
					  WHERE  vs.subsnp_id = tg.subsnp_id
                                          AND    tg.allele_1 = tra.allele
                                          AND    tg.allele_2 = tra.allele
					  AND    s.pop_id = tg.pop_id});
	
    $self->{'dbVar'}->do("DROP TABLE tmp_gty");
    $self->{'dbVar'}->do("DROP TABLE tmp_rev_allele");
}

# cleans up some of the necessary temporary data structures after the
# import is complete
sub cleanup {
    my $self = shift;

    debug("In cleanup...");
  #remove populations that are not present in the Individual or Allele table for the specie
     $self->{'dbVar'}->do('CREATE TABLE tmp_pop (sample_id int PRIMARY KEY)'); #create a temporary table with unique populations
     $self->{'dbVar'}->do('INSERT IGNORE INTO tmp_pop SELECT distinct(sample_id) FROM allele'); #add the populations from the alleles
    $self->{'dbVar'}->do('INSERT IGNORE INTO tmp_pop SELECT distinct(sample_id) FROM population_genotype'); #add the populations from the population_genotype
     $self->{'dbVar'}->do('INSERT IGNORE INTO tmp_pop SELECT population_sample_id FROM individual_population'); #add the populations from the individuals
     $self->{'dbVar'}->do(qq{INSERT IGNORE INTO tmp_pop SELECT super_population_sample_id 
 				      FROM population_structure ps, tmp_pop tp 
 				      WHERE tp.sample_id = ps.sub_population_sample_id}); #add the populations from the super-populations
    
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
    # populations not present
    $self->{'dbVar'}->do($sql_2); #delete from sample_synonym
    $self->{'dbVar'}->do($sql_3); #and delete from the population_structure table
    $self->{'dbVar'}->do($sql_4); #and delete from Sample table the ones that are not in population

    $self->{'dbVar'}->do('DROP TABLE tmp_pop'); #and finally remove the temporary table
    
    $self->{'dbVar'}->do('ALTER TABLE variation  DROP COLUMN snp_id');
    $self->{'dbVar'}->do('RENAME TABLE variation_synonym TO variation_synonym_old');
    $self->{'dbVar'}->do('CREATE TABLE variation_synonym LIKE variation_synonym_old');
    $self->{'dbVar'}->do('ALTER TABLE variation_synonym DROP COLUMN subsnp_id, DROP COLUMN substrand_reversed_flag');
    $self->{'dbVar'}->do('ALTER TABLE sample DROP COLUMN pop_class_id, DROP COLUMN pop_id, DROP COLUMN individual_id');
    $self->{'dbVar'}->do('ALTER TABLE variation_group DROP COLUMN hapset_id') if ($self->{'dbCore'}->species =~ /homo|hum/i);
    $self->{'dbVar'}->do('ALTER TABLE allele_group DROP COLUMN hap_id') if ($self->{'dbCore'}->species =~ /homo|hum/i);

}

1;

use strict;
use warnings;

#generic object for the dbSNP data. Contains the general methods to dump the data into the new Variation database. Any change in the methods
# will need to overload the correspondent method in the subclass for the specie

package dbSNP::GenericContig;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use ImportUtils qw(dumpSQL debug create_and_load load);

#creates the object and assign the attributes to it (connections, basically)
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbSNP, $dbCore, $dbVariation, $tmp_dir, $tmp_file, $limit, $alldiff_file, $taxID, $species_prefix) =
        rearrange([qw(DBSNP DBCORE DBVARIATION TMPDIR TMPFILE LIMIT ALLDIFF TAXID SPECIES_PREFIX)],@_);


  return bless {'dbSNP' => $dbSNP,
		'dbCore' => $dbCore,
		'dbVariation' => $dbVariation,
		'tmpdir' => $tmp_dir,
		'tmpfile' => $tmp_file,
		'limit' => $limit,
                'alldiff' => $alldiff_file,
		'taxID' => $taxID,
		'species_prefix' => $species_prefix}, $class;
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
    $self->variation_group();
    $self->allele_group();
    
    $self->cleanup();

}

sub source_table {
    my $self = shift;
    my ($dbname,$version) = split /\_/,$self->{'dbSNP'}->dbname(); #get the version of the dbSNP release
    $self->{'dbVariation'}->do(qq{INSERT INTO source (source_id,name,version) VALUES (1,"$dbname",$version)});

}


# filling of the variation table from SubSNP and SNP
# creating of a link table variation_id --> subsnp_id
sub variation_table {
    my $self = shift;

    $self->{'dbVariation'}->do( "ALTER TABLE variation add column snp_id int" );

    # load refSNPs into the variation table
    
    debug("Dumping RefSNPs");
    
    dumpSQL($self->{'dbSNP'},  qq{
	SELECT 1, concat( "rs", snp_id), if(validation_status = 0,NULL,validation_status), snp_id
	    FROM SNP
	    WHERE tax_id = $self->{'taxID'}
	    $self->{'limit'}
          }
	    );
    
    debug("Loading RefSNPs into variation table");
    
    load( $self->{'dbVariation'}, "variation", "source_id", "name", "validation_status", "snp_id" );
    
    $self->{'dbVariation'}->do( "ALTER TABLE variation ADD INDEX snpidx( snp_id )" );
    
    # create a temp table of subSNP info
    
    debug("Dumping SubSNPs");
    
    $self->dump_subSNPs($self->{'dbSNP'},$self->{'taxID'},$self->{'limit'},$self->{'tmpdir'},$self->{'tmpfile'});
    
    create_and_load( $self->{'dbVariation'}, "tmp_var_allele", "subsnp_id i*", "refsnp_id i*",
		     "pop_id i", "allele", "substrand_reversed_flag i");
    
    # load the synonym table with the subsnp identifiers
    
    debug("loading variation_synonym table with subsnps");
    
    $self->{'dbVariation'}->do(qq{ALTER TABLE variation_synonym add column subsnp_id int});
    $self->{'dbVariation'}->do(qq{ALTER TABLE variation_synonym add column substrand_reversed_flag tinyint});
    
    $self->{'dbVariation'}->do( qq{INSERT INTO variation_synonym (variation_id, source_id, name,
						  subsnp_id, substrand_reversed_flag )
		       SELECT v.variation_id, 1,
		       CONCAT('ss',tv.subsnp_id), tv.subsnp_id,
		       tv.substrand_reversed_flag
		       FROM tmp_var_allele tv, variation v
		       WHERE tv.refsnp_id = v.snp_id
		       GROUP BY tv.subsnp_id
		   });
    
    $self->{'dbVariation'}->do("ALTER TABLE variation_synonym ADD INDEX subsnp_id(subsnp_id)");
    
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
	(qq{SELECT subsnp.subsnp_id, subsnplink.snp_id, b.pop_id, ov.pattern,subsnplink.substrand_reversed_flag
		FROM SubSNP subsnp, SNPSubSNPLink subsnplink, ObsVariation ov, Batch b
		WHERE subsnp.batch_id = b.batch_id
		AND   subsnp.subsnp_id = subsnplink.subsnp_id
		AND   ov.var_id = subsnp.variation_id
		AND   b.tax_id = $self->{'taxID'}
	    $self->{'limit'}
	} );
    
    $sth->execute();
    
    open ( FH, ">" . $self->{'tmpdir'} . "/" . $self->{'tmpfile'} );
    
  my $row;
  while($row = $sth->fetchrow_arrayref()) {
    my @alleles = split('/', $row->[3]);

    my @row = map {(defined($_)) ? $_ : '\N'} @$row;

    # split alleles into multiple rows
    foreach my $a (@alleles) {
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
    
  $self->{'dbVariation'}->do("ALTER TABLE population ADD column pop_id int");   
  $self->{'dbVariation'}->do("ALTER TABLE population ADD column pop_class_id int"); 

  # load PopClassCode data as populations

  debug("Dumping population class data");

  dumpSQL($self->{'dbSNP'}, qq{SELECT pop_class, pop_class_id, pop_class_text
				   FROM PopClassCode 
			   });

  load($self->{'dbVariation'}, 'population', 'name', 'pop_class_id', 'description');

  $self->{'dbVariation'}->do(qq{ALTER TABLE population ADD INDEX pop_class_id (pop_class_id)});

  debug("Dumping population data");

  # load Population data as populations

  dumpSQL($self->{'dbSNP'}, qq{SELECT DISTINCT concat(p.handle, ':', p.loc_pop_id),
                    p.pop_id, pc.pop_class_id
             FROM   Population p
             LEFT JOIN PopClass pc ON p.pop_id = pc.pop_id
	 });

  debug("Loading population data");

  create_and_load( $self->{'dbVariation'}, "tmp_pop", "name", "pop_id i*", "pop_class_id i*" );

  $self->{'dbVariation'}->do(qq{INSERT INTO population (name, pop_id)
                SELECT tp.name, tp.pop_id
                FROM   tmp_pop tp
                GROUP BY tp.pop_id});

    $self->{'dbVariation'}->do(qq{ALTER TABLE population ADD INDEX pop_id (pop_id)});

    
    debug("Loading population_synonym table");

    # build super/sub population relationships
    $self->{'dbVariation'}->do(qq{INSERT INTO population_structure (super_population_id,sub_population_id)
				    SELECT p1.population_id, p2.population_id
				    FROM tmp_pop tp, population p1, population p2
				    WHERE tp.pop_class_id = p1.pop_class_id
				    AND   tp.pop_id = p2.pop_id});
    

    #load population_synonym table with dbSNP population id
    $self->{'dbVariation'}->do(qq{INSERT INTO population_synonym (population_id,source_id,name)
				      SELECT population_id, 1, pop_id
				      FROM population
				      WHERE pop_id is NOT NULL
				  });
    
    $self->{'dbVariation'}->do("DROP TABLE tmp_pop");
}


#
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
  dumpSQL($self->{'dbSNP'}, qq{ SELECT si.loc_ind_id, i.descrip, i.ind_id
				   FROM   SubmittedIndividual si, Individual i
				   WHERE  si.ind_id = i.ind_id
				   AND    si.tax_id = $self->{'taxID'}
				   GROUP BY i.ind_id
			    });

  create_and_load($self->{'dbVariation'}, 'tmp_ind', 'loc_ind_id', 'description', 'ind_id i*');

  # load pedigree into seperate tmp table because there are no
  # indexes on it in dbsnp and it makes the left join b/w tables v. slow
  # one individual has 2 (!) pedigree rows, thus the group by

  dumpSQL($self->{'dbSNP'}, qq{ SELECT ind_id, pa_ind_id, ma_ind_id, sex
              FROM PedigreeIndividual GROUP BY ind_id});

  create_and_load($self->{'dbVariation'}, 'tmp_ped', 'ind_id i*', 'pa_ind_id i', 'ma_ind_id i', 'sex');

  debug("Loading individuals into individual table");

  # to make things easier keep dbSNPs individual.ind_id as our individual_id

  $self->{'dbVariation'}->do(qq{INSERT INTO individual (individual_id, name, description,gender,father_individual_id, mother_individual_id)
				    SELECT ti.ind_id, ti.loc_ind_id, ti.description,
				    IF(tp.sex = 'M', 'Male',
				       IF(tp.sex = 'F', 'Female', 'Unknown')),
				    IF(tp.pa_ind_id > 0, tp.pa_ind_id, null),
				    IF(tp.ma_ind_id > 0, tp.ma_ind_id, null)
				    FROM   tmp_ind ti
				    LEFT JOIN tmp_ped tp ON ti.ind_id = tp.ind_id});
    
    $self->{'dbVariation'}->do("DROP table tmp_ind");
    $self->{'dbVariation'}->do("DROP table tmp_ped");
    
    #necessary to fill in the individual_population table with the relation between individual and populations
    dumpSQL($self->{'dbSNP'}, qq{ SELECT si.pop_id, i.ind_id
				      FROM   SubmittedIndividual si, Individual i
				      WHERE  si.ind_id = i.ind_id
				      AND    si.tax_id = $self->{'taxID'}
			      });
    
    create_and_load($self->{'dbVariation'}, 'tmp_ind_pop', 'pop_id i*', 'ind_id i*');

    debug("Loading individuals_population table");

    $self->{'dbVariation'}->do(qq(INSERT INTO individual_population (individual_id, population_id)
				      SELECT tip.ind_id, p.population_id
				      FROM tmp_ind_pop tip, population p
				      WHERE tip.pop_id = p.pop_id));
    
    $self->{'dbVariation'}->do("DROP table tmp_ind_pop");
    return;
}


#
# loads the allele table
#
sub allele_table {
    my $self = shift;

  debug("Dumping allele data");

    # load a temp table that can be used to reverse compliment alleles
    # we place subsnps in the same orientation as the refSNP
    dumpSQL($self->{'dbSNP'}, qq(SELECT a1.allele, a2.allele
				 FROM Allele a1, Allele a2
				 WHERE a1.rev_allele_id = a2.allele_id));

  create_and_load($self->{'dbVariation'}, "tmp_rev_allele", "allele *", "rev_allele");

  # first load the allele data for alleles that we have population and
  # frequency data for

  dumpSQL($self->{'dbSNP'}, qq(SELECT afsp.subsnp_id, afsp.pop_id, a.allele_id, a.allele, afsp.freq
			       FROM   AlleleFreqBySsPop afsp, Allele a, SubSNP ss
			       WHERE  afsp.allele_id = a.allele_id
			       AND    afsp.subsnp_id = ss.subsnp_id
			       AND    ss.tax_id = $self->{'taxID'}
			       $self->{'limit'}));

  debug("Loading allele frequency data");

  create_and_load($self->{'dbVariation'}, "tmp_allele", "subsnp_id i*", "pop_id i*",
                  "allele_id i*", "allele", "freq");

  debug("Creating allele table");

  #necessary to create a unique index to simulate the GROUP BY clause
  $self->{'dbVariation'}->do(qq{CREATE UNIQUE INDEX unique_allele_idx ON allele (variation_id,allele(2),frequency,population_id)});

  $self->{'dbVariation'}->do(qq(INSERT IGNORE INTO allele (variation_id, allele,frequency, population_id)
                SELECT vs.variation_id,
                       IF(vs.substrand_reversed_flag,
                          tra.rev_allele,tra.allele) as allele,
                       ta.freq, p.population_id
                FROM   tmp_allele ta, tmp_rev_allele tra, variation_synonym vs,
                       population p
                WHERE  ta.subsnp_id = vs.subsnp_id
                AND    ta.allele = tra.allele
		AND    ta.pop_id = p.pop_id));
    #remove the index
    $self->{'dbVariation'}->do("DROP INDEX unique_allele_idx ON allele");

  # load remaining allele data which we do not have frequency data for
  # this will not import alleles without frequency for a variation which
  # already has frequency



   $self->{'dbVariation'}->do("ALTER TABLE allele ENABLE KEYS"); #after ignoring in the insertion, we must enable keys again

   #get SNPs with 1 allele, because the frequency is 1 and there is no entry in the AFSP table for the other allele
   debug("Getting alleles with no frequency");

   $self->{'dbVariation'}->do(qq{CREATE TABLE tmp_allele_no_freq 
 				    SELECT a.variation_id, vs.subsnp_id, a.population_id, 
				        IF(vs.substrand_reversed_flag, tra.rev_allele, a.allele) as allele, vs.substrand_reversed_flag
 				    FROM allele a, variation_synonym vs, tmp_rev_allele tra
 				    WHERE a.frequency = 1
 				    AND   a.variation_id = vs.variation_id
				    AND   tra.allele = a.allele
 				    GROUP BY a.variation_id, a.population_id, vs.subsnp_id
 				    HAVING count(*) = 1
 				});

     $self->{'dbVariation'}->do(qq{INSERT INTO allele (variation_id, allele, frequency, population_id)
 				      SELECT taf.variation_id, IF(tva.substrand_reversed_flag,tra.rev_allele,tra.allele) as allele, 
				       0, taf.population_id				       
 				      FROM tmp_allele_no_freq taf, tmp_var_allele tva, tmp_rev_allele tra
 				      WHERE taf.subsnp_id = tva.subsnp_id
				      AND   tva.allele = tra.allele
 				      AND   taf.allele <> tva.allele
				      AND   taf.substrand_reversed_flag = tva.substrand_reversed_flag
				      GROUP BY taf.variation_id, allele, taf.population_id
 				  });

     debug("Loading other allele data");

   $self->{'dbVariation'}->do("DROP TABLE tmp_allele");    

   $self->{'dbVariation'}->do(qq{CREATE TABLE tmp_allele
                 SELECT vs.variation_id as variation_id, tva.pop_id,
                        IF(vs.substrand_reversed_flag,
                           tra.rev_allele, tra.allele) as allele
                 FROM   variation_synonym vs, tmp_var_allele tva,
                        tmp_rev_allele tra
                 LEFT JOIN allele a ON a.variation_id = vs.variation_id
                 WHERE  tva.subsnp_id = vs.subsnp_id
                 AND    tva.allele = tra.allele
                 AND    a.allele_id is NULL});

   $self->{'dbVariation'}->do("ALTER TABLE tmp_allele ADD INDEX pop_id(pop_id)");

   $self->{'dbVariation'}->do(qq{INSERT INTO allele (variation_id, allele,
                                     frequency, population_id)
                 SELECT ta.variation_id, ta.allele, null, p.population_id
                 FROM   tmp_allele ta
                 LEFT JOIN population p ON p.pop_id = ta.pop_id
                 GROUP BY ta.variation_id, p.population_id, ta.allele });

  $self->{'dbVariation'}->do("DROP TABLE tmp_allele_no_freq");
  $self->{'dbVariation'}->do("DROP TABLE tmp_rev_allele");
  $self->{'dbVariation'}->do("DROP TABLE tmp_var_allele");
  $self->{'dbVariation'}->do("DROP TABLE tmp_allele");
}


#
# loads the flanking sequence table
#
sub flanking_sequence_table {
  my $self = shift;

  $self->{'dbVariation'}->do(qq{CREATE TABLE tmp_seq (variation_id int,
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
				 FROM SubSNPSeq$type seq, SNPFlankStatus sfs
				 WHERE sfs.subsnp_id = seq.subsnp_id
				 $self->{'limit'}});
    

    $self->{'dbVariation'}->do(qq{CREATE TABLE tmp_seq_$type (
							      subsnp_id int,
							      line_num int,
							      line varchar(255),
							      KEY subsnp_id_idx(subsnp_id))
				  MAX_ROWS = 100000000 });
    
    load($self->{'dbVariation'}, "tmp_seq_$type", "subsnp_id", "line_num", "line");

    # merge the tables into a single tmp table
    $self->{'dbVariation'}->do(qq{INSERT INTO tmp_seq (variation_id, subsnp_id,
						       line_num, type, line, revcom)
				  SELECT vs.variation_id, ts.subsnp_id, ts.line_num, '$type',
				  ts.line, vs.substrand_reversed_flag
				  FROM   tmp_seq_$type ts, variation_synonym vs
				  WHERE  vs.subsnp_id = ts.subsnp_id});
    #drop tmp table to free space
    $self->{'dbVariation'}->do(qq{DROP TABLE tmp_seq_$type});
  }
      
  $self->{'dbVariation'}->do("ALTER TABLE tmp_seq ADD INDEX idx (subsnp_id, type, line_num)");

  my $sth = $self->{'dbVariation'}->prepare(qq{SELECT ts.variation_id, ts.subsnp_id, ts.type,
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
  $self->{'dbVariation'}->do("DROP TABLE tmp_seq");

  debug("Loading flanking sequence data");

  # import the generated data
  $self->{'dbVariation'}->do(qq{LOAD DATA LOCAL INFILE '$self->{'tmpdir'}/$self->{'tmpfile'}' INTO TABLE flanking_sequence });

  unlink($self->{'tmpdir'} . "/" . $self->{'tmpfile'});

  return;
}



sub variation_feature {
    my $self = shift;

  ### TBD not sure if variations with map_weight > 1 or 2 should be
  ### imported.

     debug("Dumping seq_region data");

     dumpSQL($self->{'dbCore'}, qq{SELECT sr.seq_region_id, sr.name
 				      FROM   seq_region sr});
    
     debug("Loading seq_region data");
     create_and_load($self->{'dbVariation'}, "tmp_seq_region", "seq_region_id", "name *");
    
     debug("Dumping SNPLoc data");
    
     my $tablename = $self->{'species_prefix'} . 'SNPContigLoc';
    
     dumpSQL($self->{'dbSNP'}, qq{SELECT snp_id, contig_acc,
 				 asn_from, 
				 IF(loc_type = 3,  asn_from - 1, asn_to), # 3 = between
 				 IF(orientation, -1, 1)
 				     FROM   $tablename
 				     $self->{'limit'}});
    
    
     debug("Loading SNPLoc data");

     create_and_load($self->{'dbVariation'}, "tmp_contig_loc", "snp_id i*", "contig *", "start i", 
 		    "end i", "strand i");
    
     debug("Creating genotyped variations");
    
     #creating the temporary table with the genotyped variations
     $self->{'dbVariation'}->do(qq{CREATE TABLE tmp_genotyped_var SELECT DISTINCT variation_id FROM individual_genotype_single_bp});
     $self->{'dbVariation'}->do(qq{INSERT IGNORE INTO  tmp_genotyped_var SELECT DISTINCT variation_id FROM individual_genotype_multiple_bp});
     $self->{'dbVariation'}->do(qq{CREATE INDEX variation_idx ON tmp_genotyped_var (variation_id)});
    
    debug("Creating tmp_variation_feature data");
    
    dumpSQL($self->{'dbVariation'},qq{SELECT v.variation_id, ts.seq_region_id, tcl.start, tcl.end, tcl.strand, v.name, v.source_id, v.validation_status
					  FROM variation v, tmp_contig_loc tcl, tmp_seq_region ts
					  WHERE v.snp_id = tcl.snp_id
					  AND ts.name = tcl.contig});
    
    create_and_load($self->{'dbVariation'},'tmp_variation_feature',"variation_id *","seq_region_id", "seq_region_start", "seq_region_end", "seq_region_strand", "variation_name", "source_id", "validation_status");
    
    debug("Dumping data into variation_feature table");
    $self->{'dbVariation'}->do(qq{INSERT INTO variation_feature (variation_id, seq_region_id,seq_region_start, seq_region_end, seq_region_strand,
								 variation_name, flags, source_id, validation_status)
				  SELECT tvf.variation_id, tvf.seq_region_id, tvf.seq_region_start, tvf.seq_region_end, tvf.seq_region_strand,
				  tvf.variation_name,IF(tgv.variation_id,'genotyped',NULL), tvf.source_id, tvf.validation_status
				  FROM tmp_variation_feature tvf LEFT JOIN tmp_genotyped_var tgv ON tvf.variation_id = tgv.variation_id
				 });
    
    $self->{'dbVariation'}->do("DROP TABLE tmp_contig_loc");
    $self->{'dbVariation'}->do("DROP TABLE tmp_seq_region");
    $self->{'dbVariation'}->do("DROP TABLE tmp_genotyped_var");
    $self->{'dbVariation'}->do("DROP TABLE tmp_variation_feature");
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
             FROM HapSet hs, HapSetSnpList hssl, SubSNP ss #, Batch b
             WHERE hs.hapset_id = hssl.hapset_id
             AND   hssl.subsnp_id = ss.subsnp_id
#             AND   ss.batch_id = b.batch_id
             AND   ss.tax_id = $self->{'taxID'}});

  create_and_load($self->{'dbVariation'}, 'tmp_var_grp', 'name', 'hapset_id i*', 'subsnp_id i*');

  $self->{'dbVariation'}->do("ALTER TABLE variation_group add column hapset_id int");

  debug("Loading variation_group");

  $self->{'dbVariation'}->do(qq{INSERT INTO variation_group (name, source_id, type, hapset_id)
                SELECT name, 1, 'haplotype', hapset_id
                FROM tmp_var_grp
                GROUP BY hapset_id});

  $self->{'dbVariation'}->do("ALTER TABLE variation_group ADD INDEX hapset_id(hapset_id)");

  debug("Loading variation_group_variation");

  # there are a few hapsets in dbSNP which have two subsnps which have been
  # merged into the same refsnp.  Thus the group by clause.

  $self->{'dbVariation'}->do(qq{INSERT INTO variation_group_variation
                     (variation_group_id, variation_id)
                SELECT vg.variation_group_id, vs.variation_id
                FROM   variation_group vg, variation_synonym vs,
                       tmp_var_grp tvg
                WHERE  tvg.hapset_id = vg.hapset_id
                AND    tvg.subsnp_id = vs.subsnp_id
                GROUP BY variation_group_id, variation_id});

  $self->{'dbVariation'}->do("DROP TABLE tmp_var_grp");
}

#
# loads allele_group table
#
sub allele_group {
    my $self = shift;

  debug("Dumping Hap data");

  dumpSQL($self->{'dbSNP'}, qq{SELECT  h.hap_id, h.hapset_id, h.loc_hap_id,
                    hsa.snp_allele, hsa.subsnp_id
             FROM   Hap h, HapSnpAllele hsa, SubSNP ss #, Batch b
             WHERE  hsa.hap_id = h.hap_id
             AND    hsa.subsnp_id = ss.subsnp_id
#             AND    ss.batch_id = b.batch_id
             AND    ss.tax_id = $self->{'taxID'}});

  create_and_load($self->{'dbVariation'}, 'tmp_allele_group_allele','hap_id i*','hapset_id i*',
                  'name','snp_allele', 'subsnp_id i*');

  $self->{'dbVariation'}->do(qq{ALTER TABLE allele_group ADD COLUMN hap_id int});

  debug("Loading allele_group");

  $self->{'dbVariation'}->do(qq{INSERT INTO allele_group (variation_group_id, name, source_id, hap_id)
                SELECT vg.variation_group_id, tag.name, 1, tag.hap_id
                FROM   variation_group vg, tmp_allele_group_allele tag
                WHERE  vg.hapset_id = tag.hapset_id
                GROUP BY hap_id});

  $self->{'dbVariation'}->do(qq{ALTER TABLE allele_group ADD INDEX hap_id(hap_id)});

  debug("Loading allele_group_allele");

  # there are a few haps in dbSNP which have two subsnps which have been
  # merged into the same refsnp.  Thus the group by clause.

  $self->{'dbVariation'}->do(qq{INSERT INTO allele_group_allele (allele_group_id,variation_id, allele)
                SELECT ag.allele_group_id, vs.variation_id, taga.snp_allele
                FROM   allele_group ag, tmp_allele_group_allele taga,
                       variation_synonym vs
                WHERE  ag.hap_id = taga.hap_id
                AND    vs.subsnp_id = taga.subsnp_id
                GROUP BY ag.allele_group_id, vs.variation_id});

  $self->{'dbVariation'}->do("DROP TABLE tmp_allele_group_allele");

  return;
}



#
# loads individual genotypes into the individual_genotype table
#
sub individual_genotypes {
    my $self = shift;

  #
  # load SubInd individual genotypes into genotype table
  #
  debug("Dumping SubInd and ObsGenotype data");
  dumpSQL($self->{'dbSNP'}, qq{SELECT si.subsnp_id, sind.ind_id, og.obs
             FROM   SubInd si, ObsGenotype og, SubmittedIndividual sind
             WHERE  og.gty_id = si.gty_id 
             AND    sind.submitted_ind_id = si.submitted_ind_id
             AND    sind.tax_id = $self->{'taxID'}
             $self->{'limit'}});

  create_and_load($self->{'dbVariation'}, "tmp_gty", 'subsnp_id i*', 'ind_id i', 'genotype');

  debug("Loading individual_genotype table");

    #we have truncated the individual_genotype table, one contains the genotypes single bp, and the other, the rest
    $self->{'dbVariation'}->do(qq{INSERT INTO individual_genotype_single_bp (variation_id, individual_id, allele_1, allele_2) 
				      SELECT vs.variation_id, tg.ind_id, SUBSTRING_INDEX(tg.genotype,'/',1) as allele_1, 
				               SUBSTRING_INDEX(tg.genotype,'/',-1) as allele_2
				      FROM   tmp_gty tg, variation_synonym vs
				      WHERE  tg.subsnp_id = vs.subsnp_id
				      AND    length(tg.genotype) = 3});

    $self->{'dbVariation'}->do(qq{INSERT INTO individual_genotype_multiple_bp (variation_id, individual_id, allele_1, allele_2) 
				      SELECT vs.variation_id, tg.ind_id, SUBSTRING_INDEX(tg.genotype,'/',1) as allele_1,
				             SUBSTRING_INDEX(tg.genotype,'/',-1) as allele_2
				      FROM   tmp_gty tg, variation_synonym vs
				      WHERE  tg.subsnp_id = vs.subsnp_id
				      AND    length(tg.genotype) > 3});
    
    $self->{'dbVariation'}->do("DROP TABLE tmp_gty");
    
    return;
}


#
# loads population genotypes into the population_genotype table
#
sub population_genotypes {
    my $self = shift;

  debug("Dumping GtyFreqBySsPop and UniGty data");
    
  #create an IN statement to select the genotypes only in the species determined by the tax
    my $sth = $self->{'dbSNP'}->prepare(qq{SELECT DISTINCT pop_id
				FROM SubmittedIndividual
				WHERE tax_id = $self->{'taxID'}     
    });
    $sth->execute();
    my $pop_ids = $sth->fetchall_arrayref([0]);
    my @ids;
    map {push @ids,$_->[0]} @{$pop_ids};
    my $in_str = " IN (" . join(',', @ids). ")";
    #only dump data if there is population genotype information
    if (defined $ids[0]){
	dumpSQL($self->{'dbSNP'},qq{SELECT gtfsp.subsnp_id, gtfsp.pop_id, gtfsp.freq,a1.allele, a2.allele
					FROM   GtyFreqBySsPop gtfsp, UniGty ug, Allele a1, Allele a2 #, SubmittedIndividual si
					WHERE  gtfsp.unigty_id = ug.unigty_id
					AND    ug.allele_id_1 = a1.allele_id
					AND    ug.allele_id_2 = a2.allele_id
					AND    gtfsp.pop_id $in_str
					$self->{'limit'}
				    });
	
	debug("loading genotype data");
	
	create_and_load($self->{'dbVariation'}, "tmp_gty", 'subsnp_id i*', 'pop_id i*', 'freq',
			'allele_1', 'allele_2');
	
	$self->{'dbVariation'}->do(qq{INSERT INTO population_genotype (variation_id,allele_1, allele_2, frequency, population_id)
					  SELECT vs.variation_id, tg.allele_1, tg.allele_2, tg.freq,
					  p.population_id
					  FROM   variation_synonym vs, tmp_gty tg, population p
					  WHERE  vs.subsnp_id = tg.subsnp_id
					  AND    p.pop_id = tg.pop_id});
	
	$self->{'dbVariation'}->do("DROP TABLE tmp_gty");
    }
}



# cleans up some of the necessary temporary data structures after the
# import is complete
sub cleanup {
    my $self = shift;

  #remove populations that are not present in the Individual or Allele table for the specie
    $self->{'dbVariation'}->do('CREATE TABLE tmp_pop (population_id int PRIMARY KEY)'); #create a temporary table with unique populations
    $self->{'dbVariation'}->do('INSERT IGNORE INTO tmp_pop SELECT population_id FROM allele'); #add the populations from the alleles
    $self->{'dbVariation'}->do('INSERT IGNORE INTO tmp_pop SELECT population_id FROM individual_population'); #add the populations from the individuals
    $self->{'dbVariation'}->do(qq{INSERT IGNORE INTO tmp_pop SELECT super_population_id 
				      FROM population_structure ps, tmp_pop tp 
				      WHERE tp.population_id = ps.sub_population_id}); #add the populations from the super-populations
    
    #necessary to difference between MySQL 4.0 and MySQL 4.1
    my $sql;
    my $sql_2;

    my $sth = $self->{'dbVariation'}->prepare(qq{SHOW VARIABLES LIKE 'version'});
    $sth->execute();
    my $row_ref = $sth->fetch();
    $sth->finish();
    #check if the value in the version contains the 4.1
    if ($row_ref->[1] =~ /4\.1/){
    
	$sql = qq{DELETE FROM p, ps USING population p 
		      LEFT JOIN tmp_pop tp ON p.population_id = tp.population_id, 
		      population_synonym ps LEFT JOIN tmp_pop tp1 on ps.population_id = tp1.population_id 
		      WHERE tp.population_id is null AND tp1.population_id is null};
	$sql_2 = qq{DELETE FROM ps USING population_structure ps 
			LEFT JOIN tmp_pop tp ON ps.sub_population_id = tp.population_id 
			WHERE tp.population_id is null};
    }
    else{

	$sql = qq{DELETE population, population_synonym FROM population p 
		      LEFT JOIN tmp_pop tp ON p.population_id = tp.population_id, 
		      population_synonym ps LEFT JOIN tmp_pop tp1 on ps.population_id = tp1.population_id 
		      WHERE tp.population_id is null AND tp1.population_id is null};
	$sql_2 = qq{DELETE population_structure FROM population_structure ps 
			LEFT JOIN tmp_pop tp ON ps.sub_population_id = tp.population_id 
			WHERE tp.population_id is null};	
    }

    $self->{'dbVariation'}->do($sql); #delete from population and population_synonym
    # populations not present
    $self->{'dbVariation'}->do($sql_2); #and delete from the population_structure table

    $self->{'dbVariation'}->do('DROP TABLE tmp_pop'); #and finally remove the temporary table
    
    $self->{'dbVariation'}->do('ALTER TABLE variation  DROP COLUMN snp_id');
    $self->{'dbVariation'}->do('ALTER TABLE variation_synonym DROP COLUMN subsnp_id');
    $self->{'dbVariation'}->do('ALTER TABLE variation_synonym DROP COLUMN substrand_reversed_flag');
    $self->{'dbVariation'}->do('ALTER TABLE population DROP COLUMN pop_class_id');
    $self->{'dbVariation'}->do('ALTER TABLE population DROP COLUMN pop_id');
    $self->{'dbVariation'}->do('ALTER TABLE variation_group DROP COLUMN hapset_id');
    $self->{'dbVariation'}->do('ALTER TABLE allele_group DROP COLUMN hap_id');

}

1;

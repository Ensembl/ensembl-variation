# this is the experimental script to fill the new variation 
# schema with data from dbSNP
# we use the local mysql copy of dbSNP at the sanger center

use strict;
use warnings;
use DBI;

my $tmp_dir = "/ecs2/scratch5/ensembl/mcvicker/dbSNP";

my $dbSNP = DBI->connect( "DBI:mysql:host=cbi2.internal.sanger.ac.uk;dbname=dbSNP_120", "dbsnpro" ) or die("Could not connect to dbSNP db: $!");

my $dbVar = DBI->connect( "DBI:mysql:host=ecs4.internal.sanger.ac.uk;dbname=mcvicker_variation;port=3352","ensadmin", "ensembl" ) or die("Could not connect to variation database: $!");

my $dbCore = DBI->connect( "DBI:mysql:host=ecs2.internal.sanger.ac.uk;dbname=mus_musculus_core_22_32b;port=3364","ensro" ) or die("Could not connect to core database: $!");


#my $TAX_ID = 9606; # human
my $TAX_ID = 10090; # mouse

my $LIMIT = '';
#my $LIMIT = ' LIMIT 100000 ';

source_table();
population_table();
individual_table();
variation_table();
individual_genotypes();
population_genotypes();
allele_table();
flanking_sequence_table();
variation_feature();
variation_group();
allele_group();

cleanup();


sub source_table {
  $dbVar->do(qq(INSERT INTO source SET source_id = 1, name = "dbSNP"));
}


# filling of the variation table from SubSNP and SNP
# creating of a link table variation_id --> subsnp_id
sub variation_table {
  $dbVar->do( "ALTER TABLE variation add column snp_id int" );

  # load refSNPs into the variation table

  debug("Dumping RefSNPs");

  dumpSQL( qq{
           SELECT 1, concat( "rs", snp_id), snp_id
           FROM SNP
           WHERE tax_id = $TAX_ID
           $LIMIT
          }
      );

  debug("Loading RefSNPs into variation table");

  load( "variation", "source_id", "name", "snp_id" );

  $dbVar->do( "ALTER TABLE variation ADD INDEX snpidx( snp_id )" );

  # create a temp table of subSNP info
  # containing RefSNP id, SubSNP id and validation status

  debug("Dumping SubSNPs");

  dump_subSNPs();

  create_and_load( "tmp_var_allele", "subsnp_id i*", "refsnp_id i*",
                   "pop_id i", "allele","valid", "substrand_reversed_flag i");

  # load the synonym table with the subsnp identifiers

  debug("loading variation_synonym table with subsnps");

  $dbVar->do("ALTER TABLE variation_synonym add column subsnp_id int");
  $dbVar->do(qq{ALTER TABLE variation_synonym 
                add column substrand_reversed_flag tinyint});
  $dbVar->do(qq{ALTER TABLE variation_synonym 
                add column validated enum ('VALIDATED', 'NOT_VALIDATED')});

  $dbVar->do( qq{INSERT INTO variation_synonym (variation_id, source_id, name,
                               subsnp_id, validated, substrand_reversed_flag )
                 SELECT v.variation_id, 1,
                        CONCAT('ss',tv.subsnp_id), tv.subsnp_id, tv.valid,
                        tv.substrand_reversed_flag
                 FROM tmp_var_allele tv, variation v
                 WHERE tv.refsnp_id = v.snp_id
                 GROUP BY tv.subsnp_id
                });

  $dbVar->do("ALTER TABLE variation_synonym ADD INDEX subsnp_id(subsnp_id)");

  # set the validation status of the RefSNPs.  A refSNP is validated if
  # it has a valid subsnp

  debug("Reloading RefSNPs with validation status set");

  my $sth = $dbVar->prepare
    (qq{SELECT v.variation_id, v.source_id, v.name, v.snp_id,
               vs.validated
        FROM   variation v, variation_synonym vs
        WHERE  vs.variation_id = v.variation_id
        ORDER BY v.variation_id, vs.validated});

  $sth->execute();

  my $cur_variation_id = undef;
  my $validated = 0;
  my $arr;

  # dump RefSNPs to tmp file with validation status set

  open ( FH, ">$tmp_dir/tabledump.txt" );

  while($arr = $sth->fetchrow_arrayref()) {
    if(!defined($cur_variation_id) || $arr->[0] != $cur_variation_id) {
      my @arr = map {(defined($_)) ? $_ : '\N' } @$arr;
      print FH join("\t", @arr), "\n";
    }
    $cur_variation_id = $arr->[0];
  }

  close(FH);

  # remove RefSNPs from db and reload them (faster than individual updates)

  $dbVar->do("DELETE FROM variation");

  load("variation", "variation_id", "source_id",
       "name", "snp_id", "validation_status");

  return;
}


#
# dumps subSNPs and associated allele information
#
sub dump_subSNPs {

  my $sth = $dbSNP->prepare
    (qq{SELECT subsnp.subsnp_id, subsnplink.snp_id, b.pop_id, ov.pattern,
             if ( subsnp.validation_status > 0, "VALIDATED", "NOT_VALIDATED" ),
             subsnplink.substrand_reversed_flag
        FROM SubSNP subsnp, SNPSubSNPLink subsnplink, ObsVariation ov,
             Batch b
        WHERE subsnp.batch_id = b.batch_id
        AND   subsnp.subsnp_id = subsnplink.subsnp_id
        AND   ov.var_id = subsnp.variation_id
        AND   b.tax_id = $TAX_ID
        $LIMIT
       } );

  $sth->execute();

  open ( FH, ">$tmp_dir/tabledump.txt" );

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

  $dbVar->do("ALTER TABLE population ADD column pop_id int");
  $dbVar->do("ALTER TABLE population ADD column pop_class_id int");

  # load PopClassCode data as populations

  debug("Dumping population class data");

  dumpSQL("SELECT pop_class, pop_class_id, pop_class_text FROM PopClassCode");

  load('population', 'name', 'pop_class_id', 'description');

  $dbVar->do(qq{ALTER TABLE population ADD INDEX pop_class_id (pop_class_id)});

  debug("Dumping population data");

  # load Population data as populations

  dumpSQL(qq{SELECT concat(p.handle, ':', p.loc_pop_id),
                    p.pop_id, pc.pop_class_id
             FROM   Population p
             LEFT JOIN PopClass pc ON p.pop_id = pc.pop_id});

  debug("Loading population data");

  create_and_load( "tmp_pop", "name", "pop_id i*", "pop_class_id i*" );

  $dbVar->do(qq{INSERT INTO population (name, pop_id)
                SELECT tp.name, tp.pop_id
                FROM   tmp_pop tp
                GROUP BY tp.pop_id});

  $dbVar->do(qq{ALTER TABLE population ADD INDEX pop_id (pop_id)});


  # build super/sub population relationships
  $dbVar->do(qq{INSERT INTO population_structure (super_population_id,
                                                  sub_population_id)
                SELECT p1.population_id, p2.population_id
                FROM tmp_pop tp, population p1, population p2
                WHERE tp.pop_class_id = p1.pop_class_id
                AND   tp.pop_id = p2.pop_id});


  $dbVar->do("DROP TABLE tmp_pop");
}


#
# loads the individual table
#
sub individual_table {
  # load individuals into the population table

  debug("Dumping Individual data");

  # a few submitted  individuals have the same individual or no individual
  # we ignore this problem with a group by
  dumpSQL(qq{ SELECT si.pop_id, si.loc_ind_id, i.descrip, i.ind_id
             FROM   SubmittedIndividual si, Individual i
             WHERE  si.ind_id = i.ind_id
             AND    i.tax_id = $TAX_ID
             GROUP BY i.ind_id});

  create_and_load('tmp_ind', 'pop_id i*', 'loc_ind_id', 'description',
                  'ind_id i*');

  # load pedigree into seperate tmp table because there are no
  # indexes on it in dbsnp and it makes the left join b/w tables v. slow
  # one individual has 2 (!) pedigree rows, thus the group by

  dumpSQL(qq{ SELECT ind_id, pa_ind_id, ma_ind_id, sex
              FROM PedigreeIndividual GROUP BY ind_id});

  create_and_load('tmp_ped', 'ind_id i*', 'pa_ind_id i', 'ma_ind_id i', 'sex');

  debug("Loading individuals into individual table");

  # to make things easier keep dbSNPs individual.ind_id as our individual_id

  $dbVar->do(qq{INSERT INTO individual (individual_id, name, description,
                                        population_id, gender,
                                    father_individual_id, mother_individual_id)
                SELECT ti.ind_id, ti.loc_ind_id, ti.description,
                       p.population_id,
                       IF(tp.sex = 'M', 'Male',
                         IF(tp.sex = 'F', 'Female', 'Unknown')),
                       IF(tp.pa_ind_id > 0, tp.pa_ind_id, null),
                       IF(tp.ma_ind_id > 0, tp.ma_ind_id, null)
                FROM   tmp_ind ti, population p
                LEFT JOIN tmp_ped tp ON ti.ind_id = tp.ind_id
                WHERE  ti.pop_id = p.pop_id});

  $dbVar->do("DROP table tmp_ind");
  $dbVar->do("DROP table tmp_ped");

  return;
}


#
# loads the allele table
#
sub allele_table {
  debug("Dumping allele data");

  # load a temp table that can be used to reverse compliment alleles
  # we place subsnps in the same orientation as the refSNP
  dumpSQL(qq(SELECT a1.allele, a2.allele
             FROM Allele a1, Allele a2
             WHERE a1.rev_allele_id = a2.allele_id));

  create_and_load("tmp_rev_allele", "allele *", "rev_allele");

  # first load the allele data for alleles that we have population and
  # frequency data for

  dumpSQL(qq(SELECT afsp.subsnp_id, afsp.pop_id, a.allele_id, a.allele,
                    afsp.freq
             FROM   AlleleFreqBySsPop afsp, Allele a, SubSNP ss, Batch b
             WHERE  afsp.allele_id = a.allele_id
             AND    afsp.subsnp_id = ss.subsnp_id
             AND    ss.batch_id = b.batch_id
             AND    b.tax_id = $TAX_ID
             $LIMIT));

  debug("Loading allele frequency data");

  create_and_load("tmp_allele", "subsnp_id i*", "pop_id i*",
                  "allele_id i*", "allele", "freq");

  debug("Creating allele table");

  $dbVar->do(qq(INSERT INTO allele (variation_id, allele,
                                    frequency, population_id)
                SELECT vs.variation_id,
                       IF(vs.substrand_reversed_flag,
                          tra.rev_allele,tra.allele) as allele,
                       ta.freq, p.population_id
                FROM   tmp_allele ta, tmp_rev_allele tra, variation_synonym vs,
                       population p
                WHERE  ta.subsnp_id = vs.subsnp_id
                AND    ta.allele = tra.allele
                AND    ta.pop_id = p.pop_id
                GROUP BY vs.variation_id, p.population_id, allele, ta.freq));

  # load remaining allele data which we do not have frequence data for
  # this will not import alleles without frequency for a variation which
  # already has frequency

  $dbVar->do("DROP TABLE tmp_allele");

  debug("Loading other allele data");


  $dbVar->do(qq{CREATE TABLE tmp_allele
                SELECT vs.variation_id as variation_id, tva.pop_id,
                       IF(vs.substrand_reversed_flag,
                          tra.rev_allele, tra.allele) as allele
                FROM   variation_synonym vs, tmp_var_allele tva,
                       tmp_rev_allele tra
                LEFT JOIN allele a ON a.variation_id = vs.variation_id
                WHERE  tva.subsnp_id = vs.subsnp_id
                AND    tva.allele = tra.allele
                AND    a.allele_id is NULL});

  $dbVar->do("ALTER TABLE tmp_allele ADD INDEX pop_id(pop_id)");

  $dbVar->do(qq{INSERT INTO allele (variation_id, allele,
                                    frequency, population_id)
                SELECT ta.variation_id, ta.allele, null, p.population_id
                FROM   tmp_allele ta
                LEFT JOIN population p ON p.pop_id = ta.pop_id
                GROUP BY ta.variation_id, p.population_id, ta.allele });

  $dbVar->do("DROP TABLE tmp_rev_allele");
  $dbVar->do("DROP TABLE tmp_var_allele");
  $dbVar->do("DROP TABLE tmp_allele");
}




#
# loads the flanking sequence table
#
sub flanking_sequence_table {
  ### TBD - need to reverse compliment flanking sequence if subsnp has
  ### reverse orientation to refsnp

  $dbVar->do(qq{CREATE TABLE tmp_seq (variation_id int,
                                      subsnp_id int,
                                      line_num int,
                                      type enum ('5','3'),
                                      line varchar(255),
                                      revcom tinyint)});

  # import both the 5prime and 3prime flanking sequence tables

  foreach my $type ('3','5') {

    debug("Dumping $type' flanking sequence");

    dumpSQL(qq{SELECT seq.subsnp_id, seq.line_num, seq.line
               FROM SubSNPSeq$type seq, Batch b, SubSNP ss
               WHERE ss.subsnp_id = seq.subsnp_id
               AND   ss.batch_id = b.batch_id
               AND   b.tax_id = $TAX_ID
               $LIMIT});
    create_and_load("tmp_seq_$type", "subsnp_id i*", "line_num i", "line");

    # merge the tables into a single tmp table
    $dbVar->do(qq{INSERT INTO tmp_seq (variation_id, subsnp_id,
                                       line_num, type, line, revcom)
                  SELECT vs.variation_id, ts.subsnp_id, ts.line_num, '$type',
                         ts.line, vs.substrand_reversed_flag
                  FROM   tmp_seq_$type ts, variation_synonym vs
                  WHERE  vs.subsnp_id = ts.subsnp_id});
  }

  $dbVar->do("ALTER TABLE tmp_seq ADD INDEX idx (variation_id, type, line_num)");

  my $sth = $dbVar->prepare(qq{SELECT ts.variation_id, ts.subsnp_id, ts.type,
                                      ts.line, ts.revcom
                               FROM   tmp_seq ts
                               ORDER BY ts.subsnp_id, ts.type, ts.line_num},
                            { mysql_use_result => 1 });

  $sth->execute();

  my ($vid, $ssid, $type, $line, $revcom);

  $sth->bind_columns(\$vid, \$ssid, \$type, \$line, \$revcom);

  open(FH, ">$tmp_dir/flankingdump.txt");

  my $upstream = '';
  my $dnstream = '';
  my $cur_vid;
  my $cur_ssid;

  debug("Rearranging flanking sequence data");

  my $longest_up = '';
  my $longest_dn = '';

  # dump sequences to file that can be imported all at once
  while($sth->fetch()) {
    if(defined($cur_ssid) && $cur_ssid != $ssid) {
      # if subsnp in reverse orientation to refsnp,
      # reverse compliment flanking sequence
      if($revcom) {
        ($upstream, $dnstream) = ($dnstream, $upstream);
        reverse_comp(\$upstream);
        reverse_comp(\$dnstream);
      }

      # take only the longest total flanking from all of the subsnps
      if(length($upstream) + length($dnstream) >
         length($longest_up) + length($longest_dn)) {
        $longest_up = $upstream;
        $longest_dn = $dnstream;
      }

      if($cur_vid != $vid) {
        $upstream = '\N' if(!$upstream); # null
        $dnstream = '\N' if(!$dnstream);
        print FH join("\t", $cur_vid, $longest_up, $longest_dn), "\n";
        $longest_up = '';
        $longest_dn = '';
      }

      $upstream = '';
      $dnstream = '';
    }
    $cur_ssid  = $ssid;
    $cur_vid   = $vid;

    if($type == 5) {
      $upstream .= $line;
    } else {
      $dnstream .= $line;
    }
  }

  # do not forget last row...
  print FH join("\t", $cur_vid, $longest_up, $longest_dn), "\n";

  $sth->finish();

  close FH;

  debug("Loading flanking sequence data");

  # import the generated data
  $dbVar->do(qq{LOAD DATA LOCAL INFILE '$tmp_dir/flankingdump.txt'
              INTO TABLE flanking_sequence});

  unlink(">$tmp_dir/flankingdump.txt");
  $dbVar->do("DROP TABLE tmp_seq_3");
  $dbVar->do("DROP TABLE tmp_seq_5");
  $dbVar->do("DROP TABLE tmp_seq");

  return;
}



sub variation_feature {

  ### TBD not sure if variations with map_weight > 1 or 2 should be
  ### imported. If they are then the map_weight needs to be set.

  debug("Dumping seq_region data");

  dumpSQL( qq{SELECT sr.seq_region_id, sr.name
              FROM   seq_region sr},
           $dbCore);

  debug("Loading seq_region data");
  create_and_load("tmp_seq_region", "seq_region_id", "name *");

  debug("Dumping SNPLoc data");
  dumpSQL( qq{SELECT snp_id, CONCAT(contig_acc, '.', contig_ver),
                     asn_from, asn_to, IF(orientation, -1, 1)
              FROM   SNPContigLoc
              $LIMIT});


  debug("Loading SNPLoc data");

  create_and_load("tmp_contig_loc", "snp_id i*", "contig *", "start i", 
                  "end i", "strand i");

  debug("Creating variation_feature data");

  $dbVar->do(qq{INSERT INTO variation_feature 
                       (variation_id, seq_region_id,
                        seq_region_start, seq_region_end, seq_region_strand,
                        variation_name)
                SELECT v.variation_id, ts.seq_region_id, tcl.start, tcl.end,
                       tcl.strand, v.name
                FROM   variation v, tmp_contig_loc tcl, tmp_seq_region ts
                WHERE  v.snp_id = tcl.snp_id
                AND    tcl.contig = ts.name});

  $dbVar->do("DROP TABLE tmp_contig_loc");
}

#
# loads variation_group and variation_group_variation tables from the
# contents of the HapSet and HapSetSnpList tables
#
sub variation_group {
  debug("Dumping HapSet data");

  dumpSQL(qq{SELECT  CONCAT(hs.handle, ':', hs.hapset_name),
                    hs.hapset_id, hssl.subsnp_id
             FROM HapSet hs, HapSetSnpList hssl, SubSNP ss, Batch b
             WHERE hs.hapset_id = hssl.hapset_id
             AND   hssl.subsnp_id = ss.subsnp_id
             AND   ss.batch_id = b.batch_id
             AND   b.tax_id = $TAX_ID});

  create_and_load('tmp_var_grp', 'name', 'hapset_id i*', 'subsnp_id i*');

  $dbVar->do("ALTER TABLE variation_group add column hapset_id int");

  debug("Loading variation_group");

  $dbVar->do(qq{INSERT INTO variation_group (name, source_id, type, hapset_id)
                SELECT name, 1, 'haplotype', hapset_id
                FROM tmp_var_grp
                GROUP BY hapset_id});

  $dbVar->do("ALTER TABLE variation_group ADD INDEX hapset_id(hapset_id)");

  debug("Loading variation_group_variation");

  $dbVar->do(qq{INSERT INTO variation_group_variation
                     (variation_group_id, variation_id)
                SELECT vg.variation_group_id, vs.variation_id
                FROM   variation_group vg, variation_synonym vs,
                       tmp_var_grp tvg
                WHERE  tvg.hapset_id = vg.hapset_id
                AND    tvg.subsnp_id = vs.subsnp_id});

  $dbVar->do("DROP TABLE tmp_var_grp");
}

#
# loads allele_group table
#
sub allele_group {
  debug("Dumping Hap data");

  dumpSQL(qq{SELECT  h.hap_id, h.hapset_id, h.loc_hap_id,
                    hsa.snp_allele, hsa.subsnp_id
             FROM   Hap h, HapSnpAllele hsa, SubSNP ss, Batch b
             WHERE  hsa.hap_id = h.hap_id
             AND    hsa.subsnp_id = ss.subsnp_id
             AND    ss.batch_id = b.batch_id
             AND    b.tax_id = $TAX_ID});

  create_and_load('tmp_allele_group_allele','hap_id i*','hapset_id i*',
                  'loc_hap_id','snp_allele', 'subsnp_id i*');

  $dbVar->do(qq{ALTER TABLE allele_group ADD COLUMN hap_id int});

  debug("Loading allele_group");

  $dbVar->do(qq{INSERT INTO allele_group (variation_group_id, name, source_id,
                                          hap_id)
                SELECT vg.variation_group_id, tag.name, 1, tag.hap_id
                FROM   variation_group vg, tmp_allele_group_allele tag
                WHERE  vg.hapset_id = tag.hapset_id
                GROUP BY hap_id});

  $dbVar->do(qq{ALTER TABLE allele_group ADD INDEX hap_id(hap_id)});

  debug("Loading allele_group_allele");

  $dbVar->do(qq{INSERT INTO allele_group_allele (allele_group_id,
                                                 variation_id, allele)
                SELECT ag.allele_group_id, vs.variation_id, taga.snp_allele
                FROM   allele_group ag, tmp_allele_group_allele taga,
                       variation_synonym vs
                WHERE  ag.hap_id = taga.hap_id
                AND    vs.subsnp_id = taga.subsnp_id});

  $dbVar->do("DROP TABLE tmp_allele_group_allele");

  return;
}



#
# loads individual genotypes into the individual_genotype table
#
sub individual_genotypes {

  #
  # load SubInd individual genotypes into genotype table
  #
  debug("Dumping SubInd and ObsGenotype data");
  dumpSQL(qq{SELECT si.subsnp_id, sind.ind_id, og.obs
             FROM   SubInd si, ObsGenotype og, SubmittedIndividual sind
             WHERE  og.gty_id = si.gty_id
             AND    sind.submitted_ind_id = si.submitted_ind_id
             $LIMIT});

  create_and_load("tmp_gty", 'subsnp_id i*', 'ind_id i', 'genotype');

  # dump to file and split apart the genotype strings
  my $sth = $dbVar->prepare(qq{SELECT vs.variation_id, tg.ind_id, tg.genotype
                               FROM   tmp_gty tg, variation_synonym vs
                               WHERE  tg.subsnp_id = vs.subsnp_id},
                            {mysql_use_result => 1});

  $sth->execute();

  open ( FH, ">$tmp_dir/tabledump.txt" );

  my $row;
  while($row = $sth->fetchrow_arrayref()) {
    my @row = @$row;
    ($row[2], $row[3]) = split('/', $row[2]);
    @row = map {(defined($_)) ? $_ : '\N'} @row;  # convert undefined to NULL;
    print FH join("\t", @row), "\n";
  }

  $sth->finish();
  close(FH);

  debug("Loading individual_genotype table");

  load("individual_genotype", 'variation_id', 'individual_id',
       'allele_1', 'allele_2');


  $dbVar->do("DROP TABLE tmp_gty");

  return;
}


#
# loads population genotypes into the population_genotype table
#
sub population_genotypes {
  debug("Dumping GtyFreqBySsPop and UniGty data");

  dumpSQL(qq{SELECT gtfsp.subsnp_id, gtfsp.pop_id, gtfsp.freq,
                    a1.allele, a2.allele
             FROM   GtyFreqBySsPop gtfsp, UniGty ug, Allele a1, Allele a2
             WHERE  gtfsp.unigty_id = ug.unigty_id
             AND    ug.allele_id_1 = a1.allele_id
             AND    ug.allele_id_2 = a2.allele_id
             $LIMIT});

  debug("loading genotype data");

  create_and_load("tmp_gty", 'subsnp_id i*', 'pop_id i*', 'freq',
                  'allele_1', 'allele_2');

  $dbVar->do(qq{INSERT INTO population_genotype (variation_id,
                                                 allele_1, allele_2,
                                                 frequency, population_id)
                SELECT vs.variation_id, tg.allele_1, tg.allele_2, tg.freq,
                       p.population_id
                FROM   variation_synonym vs, tmp_gty tg, population p
                WHERE  vs.subsnp_id = tg.subsnp_id
                AND    p.pop_id = tg.pop_id});

  $dbVar->do("DROP TABLE tmp_gty");
}



# cleans up some of the necessary temporary data structures after the
# import is complete
sub cleanup {
  $dbVar->do('ALTER TABLE variation  DROP COLUMN snp_id');
  $dbVar->do('ALTER TABLE variation_synonym DROP COLUMN subsnp_id');
  $dbVar->do('ALTER TABLE variation_synonym DROP COLUMN validated');
  $dbVar->do(qq{ALTER TABLE variation_synonym
                DROP COLUMN substrand_reversed_flag});
  $dbVar->do('ALTER TABLE population DROP COLUMN pop_class_id');
  $dbVar->do('ALTER TABLE population DROP COLUMN pop_id');
  $dbVar->do('ALTER TABLE variation_group DROP COLUMN hapset_id');
  $dbVar->do('ALTER TABLE allele_group DROP COLUMN hap_id');
  $dbVar->do("DROP TABLE tmp_seq_region");
}



# successive dumping and loading of tables is typical for this process
# dump does effectively a select into outfile without server file system access
sub dumpSQL {
  my $sql = shift;
  my $db  = shift;

  $db ||= $dbSNP;

  local *FH;

  open FH, ">$tmp_dir/tabledump.txt";

  my $sth = $db->prepare( $sql, { mysql_use_result => 1 });
  $sth->execute();
  my $first;
  while ( my $aref = $sth->fetchrow_arrayref() ) {
    my @a = map {defined($_) ? $_ : '\N'} @$aref;
    print FH join("\t", @a), "\n";
  }

  close FH;

  $sth->finish();
}


# load imports a table, optionally not all columns
# if table doesnt exist, create a varchar(255) for each column
sub load {
  my $tablename = shift;
  my @colnames = @_;

  my $cols = join( ",", @colnames );

  local *FH;
  open FH, "<$tmp_dir/tabledump.txt";
  my $sql;

  if ( @colnames ) {

    $sql = qq{
              LOAD DATA LOCAL INFILE '$tmp_dir/tabledump.txt' 
              INTO TABLE $tablename( $cols )
             };
  } else {
    $sql = qq{
              LOAD DATA LOCAL INFILE '$tmp_dir/tabledump.txt' 
              INTO TABLE $tablename
             };
  }

  $dbVar->do( $sql );

  unlink( "$tmp_dir/tabledump.txt" );
}


#
# creates a table with specified columns and loads data that was dumped
# to a tmp file into the table.
#
# by default all columns are VARCHAR(255), but an 'i' may be added after the
# column name to make it an INT.  Additionally a '*' means add an index to
# the column.
#
# e.g.  create_and_load('mytable', 'col0', 'col1 *', 'col2 i', 'col3 i*');
#
sub create_and_load {
  my $tablename = shift;
  my @cols = @_;

  my $sql = "CREATE TABLE $tablename ( ";

  my @col_defs;
  my @idx_defs;
  my @col_names;

  foreach my $col (@cols) {
    my ($name, $type) = split(/\s+/,$col);
    push @col_names, $name;

    if(defined($type) && $type =~ /i/) {
      push @col_defs, "$name INT";
    } else {
      push @col_defs, "$name VARCHAR(255)";
    }

    if(defined($type) && $type =~ /\*/) {
      push @idx_defs, "KEY ${name}_idx($name)";
    }
  }

  my $create_cols = join( ",\n", @col_defs, @idx_defs);


  $sql .= $create_cols.")";

  $dbVar->do( $sql );

  load( $tablename, @col_names );
}



sub debug {
  print STDERR @_, "\n";
}


#
# prints number of rows in a given table, used for debugging
#

sub count_rows {
  my $tablename = shift;

  my ($count) = $dbVar->selectall_arrayref
                    ("SELECT count(*) FROM $tablename")->[0]->[0];

  print STDERR "table $tablename has $count rows\n";
}


#
# reverse compliments nucleotide sequence
#

sub reverse_comp {
  my $seqref = shift;

  $$seqref = reverse( $$seqref );
  $$seqref =~
    tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;

  return;
}

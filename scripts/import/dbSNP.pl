#!/usr/local/ensembl/bin/perl -w
# this is the experimental script to fill the new variation 
# schema with data from dbSNP
# we use the local mysql copy of dbSNP at the sanger center

use strict;
use warnings;
use DBI;
use DBH;
use Getopt::Long;
use Benchmark;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use ImportUtils qw(dumpSQL debug create_and_load load );

my ($TAX_ID, $LIMIT_SQL, $CONTIG_SQL, $TMP_DIR, $TMP_FILE);

my $dbSNP;
my $dbVar;
my $dbCore;

{
  my($dshost, $dsuser, $dspass, $dsport, $dsdbname, # dbSNP db
     $chost, $cuser, $cpass, $cport, $cdbname,      # ensembl core db
     $vhost, $vuser, $vpass, $vport, $vdbname,      # ensembl variation db
     $limit);

  GetOptions('dshost=s'   => \$dshost,
             'dsuser=s'   => \$dsuser,
             'dspass=s'   => \$dspass,
             'dsport=i'   => \$dsport,
             'dsdbname=s' => \$dsdbname,
             'chost=s'   => \$chost,
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
             'limit=i'   => \$limit);

  $dshost   ||= 'cbi2.internal.sanger.ac.uk';
  $dsdbname ||= 'dbSNP_121';
  $dsuser   ||= 'dbsnpro';
  $dsport   ||= 3306;

  $chost    ||= 'ecs2';
  $cuser    ||= 'ensro';
  $cport    ||= 3364;

  $vport    ||= 3362;
  $vuser    ||= 'ensadmin';

  usage('-cdbname argument is required.') if(!$cdbname);
  usage('-vdbname argument is required.') if(!$vdbname);

  $TMP_DIR  = $ImportUtils::TMP_DIR;
  $TMP_FILE = $ImportUtils::TMP_FILE;


  $LIMIT_SQL = ($limit) ? " LIMIT $limit " : '';

  $dbSNP = DBH->connect
    ("DBI:mysql:host=$dshost;dbname=$dsdbname;port=$dsport",$dsuser, $dspass,
    {'RaiseError' => 1});
  die("Could not connect to dbSNP db: $!") if(!$dbSNP);

  $dbSNP->{mysql_auto_reconnect} = 1;

  $dbSNP->do("SET SESSION wait_timeout = 2678200");

  $dbVar = DBH->connect
    ("DBI:mysql:host=$vhost;dbname=$vdbname;port=$vport",$vuser, $vpass,
    {'RaiseError' => 1});
  die("Could not connect to variation database: $!") if(!$dbVar);

  $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host   => $chost,
     -user   => $cuser,
     -pass   => $cpass,
     -port   => $cport,
     -dbname => $cdbname);

  my $mc = $dbCore->get_MetaContainer();
  my $species = $mc->get_Species();

  throw("Unable to determine species from core database.") if(!$species);

  $TAX_ID = $mc->get_taxonomy_id();


  if($species->binomial() eq 'Homo sapiens') {
    $CONTIG_SQL = ' CONCAT(contig_acc, ".", contig_ver) ';
  } else {
    $CONTIG_SQL = ' contig_acc ';
  }
}

my $SPECIES_PREFIX = get_species_prefix($TAX_ID);



#source_table();
#population_table();
#individual_table();
#variation_table();
#individual_genotypes();
#population_genotypes();
#allele_table();
#flanking_sequence_table();
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

  dumpSQL($dbSNP,  qq{
           SELECT 1, concat( "rs", snp_id), validation_status, snp_id
           FROM SNP
           WHERE tax_id = $TAX_ID
           $LIMIT_SQL
          }
      );

  debug("Loading RefSNPs into variation table");

  load( $dbVar, "variation", "source_id", "name", "validation_status", "snp_id" );

  $dbVar->do( "ALTER TABLE variation ADD INDEX snpidx( snp_id )" );

  # create a temp table of subSNP info

  debug("Dumping SubSNPs");

  dump_subSNPs();

  create_and_load( $dbVar, "tmp_var_allele", "subsnp_id i*", "refsnp_id i*",
                   "pop_id i", "allele", "substrand_reversed_flag i");

  # load the synonym table with the subsnp identifiers

  debug("loading variation_synonym table with subsnps");

  $dbVar->do("ALTER TABLE variation_synonym add column subsnp_id int");
  $dbVar->do(qq{ALTER TABLE variation_synonym 
                add column substrand_reversed_flag tinyint});

  $dbVar->do( qq{INSERT INTO variation_synonym (variation_id, source_id, name,
                               subsnp_id, substrand_reversed_flag )
                 SELECT v.variation_id, 1,
                        CONCAT('ss',tv.subsnp_id), tv.subsnp_id,
                        tv.substrand_reversed_flag
                 FROM tmp_var_allele tv, variation v
                 WHERE tv.refsnp_id = v.snp_id
                 GROUP BY tv.subsnp_id
                });

  $dbVar->do("ALTER TABLE variation_synonym ADD INDEX subsnp_id(subsnp_id)");

  ### FIX: Not sure if all RefSNPs have subsnps, and if ones which do not
  ### should possibly be eliminated

  return;
}


#
# dumps subSNPs and associated allele information
#
sub dump_subSNPs {

  my $sth = $dbSNP->prepare
    (qq{SELECT subsnp.subsnp_id, subsnplink.snp_id, b.pop_id, ov.pattern,
             subsnplink.substrand_reversed_flag
        FROM SubSNP subsnp, SNPSubSNPLink subsnplink, ObsVariation ov,
             Batch b
        WHERE subsnp.batch_id = b.batch_id
        AND   subsnp.subsnp_id = subsnplink.subsnp_id
        AND   ov.var_id = subsnp.variation_id
        AND   b.tax_id = $TAX_ID
        $LIMIT_SQL
       } );

  $sth->execute();

  open ( FH, ">$TMP_DIR/$TMP_FILE" );

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

  dumpSQL($dbSNP, "SELECT pop_class, pop_class_id, pop_class_text FROM PopClassCode");

  load($dbVar, 'population', 'name', 'pop_class_id', 'description');

  $dbVar->do(qq{ALTER TABLE population ADD INDEX pop_class_id (pop_class_id)});

  debug("Dumping population data");

  # load Population data as populations

  dumpSQL($dbSNP, qq{SELECT concat(p.handle, ':', p.loc_pop_id),
                    p.pop_id, pc.pop_class_id
             FROM   Population p
             LEFT JOIN PopClass pc ON p.pop_id = pc.pop_id});

  debug("Loading population data");

  create_and_load( $dbVar, "tmp_pop", "name", "pop_id i*", "pop_class_id i*" );

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
  dumpSQL($dbSNP, qq{ SELECT si.pop_id, si.loc_ind_id, i.descrip, i.ind_id
             FROM   SubmittedIndividual si, Individual i
             WHERE  si.ind_id = i.ind_id
             AND    i.tax_id = $TAX_ID
             GROUP BY i.ind_id});

  create_and_load($dbVar, 'tmp_ind', 'pop_id i*', 'loc_ind_id', 'description',
                  'ind_id i*');

  # load pedigree into seperate tmp table because there are no
  # indexes on it in dbsnp and it makes the left join b/w tables v. slow
  # one individual has 2 (!) pedigree rows, thus the group by

  dumpSQL($dbSNP, qq{ SELECT ind_id, pa_ind_id, ma_ind_id, sex
              FROM PedigreeIndividual GROUP BY ind_id});

  create_and_load($dbVar, 'tmp_ped', 'ind_id i*', 'pa_ind_id i', 'ma_ind_id i', 'sex');

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
  dumpSQL($dbSNP, qq(SELECT a1.allele, a2.allele
             FROM Allele a1, Allele a2
             WHERE a1.rev_allele_id = a2.allele_id));

  create_and_load($dbVar, "tmp_rev_allele", "allele *", "rev_allele");

  # first load the allele data for alleles that we have population and
  # frequency data for

  dumpSQL($dbSNP, qq(SELECT afsp.subsnp_id, afsp.pop_id, a.allele_id, a.allele,
                    afsp.freq
             FROM   AlleleFreqBySsPop afsp, Allele a, SubSNP ss, Batch b
             WHERE  afsp.allele_id = a.allele_id
             AND    afsp.subsnp_id = ss.subsnp_id
             AND    ss.batch_id = b.batch_id
             AND    b.tax_id = $TAX_ID
             $LIMIT_SQL));

  debug("Loading allele frequency data");

  create_and_load($dbVar, "tmp_allele", "subsnp_id i*", "pop_id i*",
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
  $dbVar->do(qq{CREATE TABLE tmp_seq (variation_id int,
                                      subsnp_id int,
                                      line_num int,
                                      type enum ('5','3'),
                                      line varchar(255),
                                      revcom tinyint)
                MAX_ROWS = 100000000});

  # import both the 5prime and 3prime flanking sequence tables

  foreach my $type ('3','5') {

    debug("Dumping $type' flanking sequence");

    dumpSQL($dbSNP, qq{SELECT seq.subsnp_id, seq.line_num, seq.line
               FROM SubSNPSeq$type seq, Batch b, SubSNP ss
               WHERE ss.subsnp_id = seq.subsnp_id
               AND   ss.batch_id = b.batch_id
               AND   b.tax_id = $TAX_ID
               $LIMIT_SQL});


    $dbVar->do(qq{CREATE TABLE tmp_seq_$type (
                     subsnp_id int,
                     line_num int,
                     line varchar(255),
                     KEY subsnp_id_idx(subsnp_id))
                  MAX_ROWS = 100000000 });

    load($dbVar, "tmp_seq_$type", "subsnp_id", "line_num", "line");

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
                               ORDER BY ts.subsnp_id, ts.type, ts.line_num},{mysql_use_result => 1});

  $sth->execute();

  my ($vid, $ssid, $type, $line, $revcom);

  $sth->bind_columns(\$vid, \$ssid, \$type, \$line, \$revcom);

  open(FH, ">$TMP_DIR/$TMP_FILE");

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
  print FH join("\t", $cur_vid, $longest_up, $longest_dn), "\n" if ($longest_up ne '');
  #only to test when there is 1 ssID in the temp table
  print FH join("\t", $cur_vid, $upstream, $dnstream), "\n" if ($longest_up eq '');

  $sth->finish();

  close FH;

  debug("Loading flanking sequence data");

  # import the generated data
  $dbVar->do(qq{LOAD DATA LOCAL INFILE '$TMP_DIR/$TMP_FILE' INTO TABLE flanking_sequence });

  unlink("$TMP_DIR/$TMP_FILE");
  $dbVar->do("DROP TABLE tmp_seq_3");
  $dbVar->do("DROP TABLE tmp_seq_5");
  $dbVar->do("DROP TABLE tmp_seq");

  return;
}



sub variation_feature {

  ### TBD not sure if variations with map_weight > 1 or 2 should be
  ### imported.

  debug("Dumping seq_region data");

  dumpSQL($dbCore, qq{SELECT sr.seq_region_id, sr.name
              FROM   seq_region sr});

  debug("Loading seq_region data");
  create_and_load($dbVar, "tmp_seq_region", "seq_region_id", "name *");

  debug("Dumping SNPLoc data");

  my $tablename = $SPECIES_PREFIX . 'SNPContigLoc';

  dumpSQL($dbSNP, qq{SELECT snp_id, $CONTIG_SQL,
                     asn_from, 
                     IF(loc_type = 3,  asn_from - 1, asn_to), # 3 = between
                     IF(orientation, -1, 1)
              FROM   $tablename
              $LIMIT_SQL});


  debug("Loading SNPLoc data");

  create_and_load($dbVar, "tmp_contig_loc", "snp_id i*", "contig *", "start i", 
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

  dumpSQL($dbSNP, qq{SELECT  CONCAT(hs.handle, ':', hs.hapset_name),
                    hs.hapset_id, hssl.subsnp_id
             FROM HapSet hs, HapSetSnpList hssl, SubSNP ss, Batch b
             WHERE hs.hapset_id = hssl.hapset_id
             AND   hssl.subsnp_id = ss.subsnp_id
             AND   ss.batch_id = b.batch_id
             AND   b.tax_id = $TAX_ID});

  create_and_load($dbVar, 'tmp_var_grp', 'name', 'hapset_id i*', 'subsnp_id i*');

  $dbVar->do("ALTER TABLE variation_group add column hapset_id int");

  debug("Loading variation_group");

  $dbVar->do(qq{INSERT INTO variation_group (name, source_id, type, hapset_id)
                SELECT name, 1, 'haplotype', hapset_id
                FROM tmp_var_grp
                GROUP BY hapset_id});

  $dbVar->do("ALTER TABLE variation_group ADD INDEX hapset_id(hapset_id)");

  debug("Loading variation_group_variation");

  # there are a few hapsets in dbSNP which have two subsnps which have been
  # merged into the same refsnp.  Thus the group by clause.

  $dbVar->do(qq{INSERT INTO variation_group_variation
                     (variation_group_id, variation_id)
                SELECT vg.variation_group_id, vs.variation_id
                FROM   variation_group vg, variation_synonym vs,
                       tmp_var_grp tvg
                WHERE  tvg.hapset_id = vg.hapset_id
                AND    tvg.subsnp_id = vs.subsnp_id
                GROUP BY variation_group_id, variation_id});

  $dbVar->do("DROP TABLE tmp_var_grp");
}

#
# loads allele_group table
#
sub allele_group {
  debug("Dumping Hap data");

  dumpSQL($dbSNP, qq{SELECT  h.hap_id, h.hapset_id, h.loc_hap_id,
                    hsa.snp_allele, hsa.subsnp_id
             FROM   Hap h, HapSnpAllele hsa, SubSNP ss, Batch b
             WHERE  hsa.hap_id = h.hap_id
             AND    hsa.subsnp_id = ss.subsnp_id
             AND    ss.batch_id = b.batch_id
             AND    b.tax_id = $TAX_ID});

  create_and_load($dbVar, 'tmp_allele_group_allele','hap_id i*','hapset_id i*',
                  'name','snp_allele', 'subsnp_id i*');

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

  # there are a few haps in dbSNP which have two subsnps which have been
  # merged into the same refsnp.  Thus the group by clause.

  $dbVar->do(qq{INSERT INTO allele_group_allele (allele_group_id,
                                                 variation_id, allele)
                SELECT ag.allele_group_id, vs.variation_id, taga.snp_allele
                FROM   allele_group ag, tmp_allele_group_allele taga,
                       variation_synonym vs
                WHERE  ag.hap_id = taga.hap_id
                AND    vs.subsnp_id = taga.subsnp_id
                GROUP BY ag.allele_group_id, vs.variation_id});

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
  dumpSQL($dbSNP, qq{SELECT si.subsnp_id, sind.ind_id, og.obs
             FROM   SubInd si, ObsGenotype og, SubmittedIndividual sind
             WHERE  og.gty_id = si.gty_id
             AND    sind.submitted_ind_id = si.submitted_ind_id
             $LIMIT_SQL});

  create_and_load($dbVar, "tmp_gty", 'subsnp_id i*', 'ind_id i', 'genotype');

  # dump to file and split apart the genotype strings
  my $sth = $dbVar->prepare(qq{SELECT vs.variation_id, tg.ind_id, tg.genotype
                               FROM   tmp_gty tg, variation_synonym vs
                               WHERE  tg.subsnp_id = vs.subsnp_id});

  $sth->execute();

  open ( FH, ">$TMP_DIR/$TMP_FILE" );

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

  load($dbVar, "individual_genotype", 'variation_id', 'individual_id',
       'allele_1', 'allele_2');


  $dbVar->do("DROP TABLE tmp_gty");

  return;
}


#
# loads population genotypes into the population_genotype table
#
sub population_genotypes {
  debug("Dumping GtyFreqBySsPop and UniGty data");

  dumpSQL($dbSNP,qq{SELECT gtfsp.subsnp_id, gtfsp.pop_id, gtfsp.freq,
                    a1.allele, a2.allele
             FROM   GtyFreqBySsPop gtfsp, UniGty ug, Allele a1, Allele a2
             WHERE  gtfsp.unigty_id = ug.unigty_id
             AND    ug.allele_id_1 = a1.allele_id
             AND    ug.allele_id_2 = a2.allele_id
             $LIMIT_SQL});

  debug("loading genotype data");

  create_and_load($dbVar, "tmp_gty", 'subsnp_id i*', 'pop_id i*', 'freq',
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
  $dbVar->do(qq{ALTER TABLE variation_synonym
                DROP COLUMN substrand_reversed_flag});
  $dbVar->do('ALTER TABLE population DROP COLUMN pop_class_id');
  $dbVar->do('ALTER TABLE population DROP COLUMN pop_id');
  $dbVar->do('ALTER TABLE variation_group DROP COLUMN hapset_id');
  $dbVar->do('ALTER TABLE allele_group DROP COLUMN hap_id');
  $dbVar->do("DROP TABLE tmp_seq_region");
}



sub get_species_prefix {
  my $tax_id = shift;

  my $arr = $dbSNP->selectall_arrayref
    (qq{SELECT ot.prefix
        FROM   OrganismTax ot
        WHERE  ot.tax_id = $tax_id});

  if(@$arr) {
    return $arr->[0]->[0];
  }

  warn("tax_id=$tax_id not found in OrganismTax table." .
       "Assuming no species prefix");
  return '';
}




sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl dbSNP.pl <options>

options:
    -dshost <hostname>   hostname of dbSNP MySQL database (default = cbi2.internal.sanger.ac.uk)
    -dsuser <user>       username of dbSNP MySQL database (default = dbsnpro)
    -dspass <pass>       password of dbSNP MySQL database
    -dsport <port>       TCP port of dbSNP MySQL database (default = 3306)
    -dsdbname <dbname>   dbname of dbSNP MySQL database   (default = dbSNP_121)
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
    -limit <number>      limit the number of rows transfered for testing
    -tmpdir <dir>        temporary directory to use (with lots of space!)
    -tmpfile <filename>  temporary filename to use
EOF

  die("\n$msg\n\n");
}

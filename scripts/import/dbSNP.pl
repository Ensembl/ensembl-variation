# this is the experimental script to fill the new variation 
# schema with data from dbSNP
# we use the local mysql copy of dbSNP at the sanger center

use strict;
use DBI;

my $tmp_dir = "/ecs2/scratch5/ensembl/mcvicker/dbSNP";

my $dbSNP = DBI->connect( "DBI:mysql:host=cbi2.internal.sanger.ac.uk;dbname=dbSNP_120", "dbsnpro" );

my $dbVar = DBI->connect( "DBI:mysql:host=ecs4.internal.sanger.ac.uk;dbname=mcvicker_variation;port=3352","ensadmin", "ensembl" );

my $dbCore = DBI->connect( "DBI:mysql:host=ecs2.internal.sanger.ac.uk;dbname=homo_sapiens_core_22_34d;port=3364","ensro" );


source_table();
variation_table();
population_table();
allele_table();
flanking_sequence_table();
variation_feature();
cleanup();


sub source_table {
  $dbVar->do(qq(INSERT INTO source SET source_id = 1, name = "dbSNP"));
}


# filling of the variation table from SubSNP and SNP
# creating of a link table variation_id --> subsnp_id
sub variation_table {

  ### TBD make sure only the correct SNPs for the species come across
  ### by looking at the tax_id

  $dbVar->do( "ALTER TABLE variation add column snp_id int" );
  $dbVar->do( "ALTER TABLE variation add column subsnp_id int" );

  # load refSNPs into the variation table

  debug("Dumping RefSNPs");

  dumpSQL( qq{
           SELECT 1,concat( "rs", snp_id), snp_id
           FROM SNP
           LIMIT 10000
          }
      );

  debug("Loading RefSNPs into variation table");

  load( "variation", "source_id", "name", "snp_id" );

  debug("Adding index to variation table");

  $dbVar->do( "ALTER TABLE variation ADD INDEX snpidx( snp_id )" );

  # create a temp table of subSNP info 
  # containing RefSNP id, SubSNP id and validation status

  debug("Dumping SubSNPs");

  dumpSQL( qq{
           SELECT subsnp.subsnp_id, subsnplink.snp_id,
           if ( subsnp.validation_status > 0, "VALIDATED", "NOT_VALIDATED" )
           FROM SubSNP subsnp, SNPSubSNPLink subsnplink
           WHERE subsnp.subsnp_id = subsnplink.subsnp_id
           LIMIT 10000
          } );

  create_and_load( "tmp_var", "subsnp", "refsnp", "valid" );

  $dbVar->do( qq{
                 ALTER TABLE tmp_var MODIFY subsnp int 
                } );
  $dbVar->do( qq{
                 ALTER TABLE tmp_var add INDEX subsnp_idx( subsnp )
                } );

  debug("Building SubSNP parent ids");

  # create a second temp table containing the parent_id of the subsnps
  $dbVar->do( qq{
                 CREATE TABLE tmp_var2
                 SELECT tv.subsnp, v.variation_id, tv.valid
                 FROM tmp_var tv, variation v
                 WHERE tv.refsnp = v.snp_id
                });


  debug("Loading SubSNPs into variation table");
  # load the SubSNPs into the variation table
  $dbVar->do( qq{
                 INSERT INTO variation( source_id, name, parent_variation_id,
                                        validation_status, subsnp_id )
                 SELECT 1, concat("ss",subsnp), variation_id, valid, subsnp
                 FROM tmp_var2
                } );

  $dbVar->do('DROP table tmp_var');
  $dbVar->do('DROP table tmp_var2');

  return;
}


#
# loads the population table
#
sub population_table {

  ### TBD create new parent pops for the few populations that have multiple
  ### PopClassCodes in dbSNP, either that or just use the MULTINATIONAL
  ### class

  ### TBD add description

  $dbVar->do("ALTER TABLE population ADD column pop_id int");
  $dbVar->do("ALTER TABLE population ADD column pop_class_id int");

  # load PopClassCode data into tmp table

  debug("Dumping population class data");

  dumpSQL("SELECT pop_class, pop_class_id FROM PopClassCode");

  debug("Loading population class data");

  create_and_load( "tmp_pop_class", "name", "pop_class_id");

  $dbVar->do("ALTER TABLE tmp_pop_class MODIFY pop_class_id int");
  $dbVar->do(qq{ALTER TABLE tmp_pop_class 
                ADD INDEX pop_class_id (pop_class_id)});

  debug("Dumping population data");

  # load Population data into tmp table

  dumpSQL(qq{SELECT concat(p.handle, ':', p.loc_pop_id),
                    p.pop_id, pc.pop_class_id
             FROM   Population p, PopClass pc
             WHERE  p.pop_id = pc.pop_id });

  debug("Loading population data");

  create_and_load( "tmp_pop", "name", "pop_id", "pop_class_id");

  $dbVar->do("ALTER TABLE tmp_pop MODIFY pop_class_id int");
  $dbVar->do("ALTER TABLE tmp_pop ADD INDEX pop_class_id(pop_class_id)");
  $dbVar->do("ALTER TABLE tmp_pop MODIFY pop_id int");
  $dbVar->do("ALTER TABLE tmp_pop ADD INDEX pop_id(pop_id)");


  debug("Creating population table");

  # load PopClasses as parent populations
  $dbVar->do(qq{INSERT INTO population (name, pop_class_id)
                SELECT name, pop_class_id from tmp_pop_class});

  # set the parent pop_id of the other populations

  $dbVar->do(qq(CREATE table tmp_pop2
                SELECT tp.name as name, tp.pop_id as pop_id,
                       p.population_id as parent_population_id
                FROM   tmp_pop tp, population p
                WHERE  tp.pop_class_id = p.pop_class_id));

  $dbVar->do(qq{INSERT INTO population (name, parent_population_id, pop_id)
                SELECT name, parent_population_id, pop_id
                FROM   tmp_pop2});

  $dbVar->do("DROP TABLE tmp_pop_class");
  $dbVar->do("DROP TABLE tmp_pop");
  $dbVar->do("DROP TABLE tmp_pop2");

  return;
}



#
# loads the allele table
#
sub allele_table {
  debug("Dumping allele data");

  dumpSQL(qq(SELECT afsp.subsnp_id, afsp.pop_id, a.allele_id, a.allele,
                    afsp.freq
             FROM   AlleleFreqBySsPop afsp, Allele a
             WHERE  afsp.allele_id = a.allele_id
             LIMIT  10000));

  debug("Loading allele data");

  create_and_load("tmp_allele", "subsnp_id", "pop_id",
                  "allele_id", "allele", "freq");

  $dbVar->do("ALTER TABLE tmp_allele MODIFY subsnp_id int");
  $dbVar->do("ALTER TABLE tmp_allele MODIFY pop_id    int");
  $dbVar->do("ALTER TABLE tmp_allele MODIFY allele_id int");

  $dbVar->do("ALTER TABLE tmp_allele ADD INDEX subsnp_id(subsnp_id)");
  $dbVar->do("ALTER TABLE tmp_allele ADD INDEX pop_id(pop_id)");

  debug("Creating allele table");

  $dbVar->do(qq(INSERT INTO allele (variation_id, allele,
                                    frequency, population_id)
                SELECT v.variation_id, ta.allele, ta.freq,
                       p.population_id
                FROM   tmp_allele ta, variation v, population p
                WHERE  ta.subsnp_id = v.subsnp_id
                AND    ta.pop_id    = p.pop_id));

  $dbVar->do("DROP TABLE tmp_allele");
}




#
# loads the flanking sequence table
#
sub flanking_sequence_table {
  $dbVar->do(qq{CREATE TABLE tmp_seq (subsnp_id int,
                                      line_num int,
                                      type enum ('5','3'),
                                      line varchar(255))});

  # import both the 5prime and 3prime flanking sequence tables

  foreach my $type ('3','5') {

    debug("Dumping $type' flanking sequence");

    dumpSQL(qq{SELECT subsnp_id, line_num, line
               FROM SubSNPSeq$type
               LIMIT 10000});
    create_and_load("tmp_seq_$type", "subsnp_id", "line_num", "line");
    $dbVar->do("ALTER TABLE tmp_seq_$type MODIFY subsnp_id int");
    $dbVar->do("ALTER TABLE tmp_seq_$type MODIFY line_num int");

    # merge the tables into a single tmp table
    $dbVar->do(qq{INSERT INTO tmp_seq (subsnp_id, line_num, type, line)
                  SELECT v.variation_id, ts.line_num, '$type', ts.line
                  FROM   tmp_seq_$type ts, variation v
                  WHERE  v.subsnp_id = ts.subsnp_id});
  }

  $dbVar->do("ALTER TABLE tmp_seq ADD INDEX idx (subsnp_id, type, line_num)");

  my $sth = $dbVar->prepare(qq{SELECT v.variation_id, ts.type, ts.line
                               FROM   tmp_seq ts, variation v
                               WHERE  v.subsnp_id = ts.subsnp_id
                               ORDER BY v.variation_id, ts.type, ts.line_num});

  $sth->execute();

  my ($vid, $type, $line);

  $sth->bind_columns(\$vid, \$type, \$line);

  open(FH, ">$tmp_dir/flankingdump.txt");

  my $upstream = '';
  my $dnstream = '';
  my $cur_vid;

  debug("Rearranging flanking sequence data");

  # dump sequences to file that can be imported all at once
  while($sth->fetch()) {
    if(defined($cur_vid) && $cur_vid != $vid) {
        $upstream = '\N' if(!$upstream); # null
        $dnstream = '\N' if(!$dnstream);
        print FH join("\t", $cur_vid, $upstream, $dnstream), "\n";
        $upstream = '';
        $dnstream = '';
    }
    $cur_vid  = $vid;

    if($type == 5) {
      $dnstream .= $line;
    } else {
      $upstream .= $line;
    }
  }
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

  debug("Dumping seq_region data");

  dumpSQL( qq{SELECT sr.seq_region_id, sr.name
              FROM   seq_region sr},
           $dbCore);

  debug("Loading seq_region data");
  create_and_load("tmp_seq_region", "seq_region_id", "name");

  $dbVar->do("ALTER TABLE tmp_seq_region ADD INDEX name(name)");


  debug("Dumping SNPLoc data");
  dumpSQL( qq{SELECT snp_id, CONCAT(contig_acc, '.', contig_ver),
                     asn_from, asn_to, orientation
              FROM   SNPContigLoc
              LIMIT 10000});


  debug("Loading SNPLoc data");

  create_and_load("tmp_contig_loc", "snp_id", "contig", "start", "end",
                  "strand");

  $dbVar->do("ALTER TABLE tmp_contig_loc ADD INDEX contig_idx(contig)");

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

  # need to fix strands of zero, push snps to toplevel and check
  # that flanking sequence seems right

}





# cleans up some of the necessary temporary data structures after the
# import is complete
sub cleanup {
  $dbVar->do('ALTER TABLE variation  DROP COLUMN snp_id');
  $dbVar->do('ALTER TABLE variation  DROP COLUMN subsnp_id');
  $dbVar->do('ALTER TABLE population DROP COLUMN pop_class_id');
  $dbVar->do('ALTER TABLE population DROP COLUMN pop_id');
}



# successive dumping and loading of tables is typical for this process
# dump does effectively a select into outfile without server file system access
sub dumpSQL {
  my $sql = shift;
  my $db  = shift;

  $db ||= $dbSNP;

  local *FH;

  open ( FH, ">$tmp_dir/tabledump.txt" );

  my $sth = $db->prepare( $sql, { mysql_use_result => 1 });
  $sth->execute();
  my $first;
  while ( my $aref = $sth->fetchrow_arrayref() ) {
    $first = 1;
    for my $col ( @$aref ) {
      if ( $first ) {
        $first = 0;
      } else {
        print FH "\t";
      }
      if ( defined $col ) {
        print FH $col;
      } else {
        print FH "\\N";
      }
    }
    print FH "\n";
  }

  $sth->finish();
}


# load imports a table, optionally not all columns
# if table doesnt exist, create a varchar(255) for each column
sub load {
  my $tablename = shift;
  my @colnames = @_;

  my $cols = join( ",", @colnames );

  local *FH;
  open( FH, "<$tmp_dir/tabledump.txt" );
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


sub create_and_load {
  my $tablename = shift;
  my @cols = @_;

  my $sql = "CREATE TABLE $tablename ( ";
  my $create_cols = join( ",\n", map { "$_ varchar(255)" } @cols );
  $sql .= $create_cols.")";
  $dbVar->do( $sql );

  load( $tablename, @cols );
}



sub debug {
  print STDERR @_, "\n";
}

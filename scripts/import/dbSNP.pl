# this is the experimental script to fill the new variation 
# schema with data from dbSNP
# we use the local mysql copy of dbSNP at the sanger center



use strict;
use DBI;

my $tmp_dir = "/acari/scratch1/ensembl/arne/db";

my $dbSNP = DBI->connect( "DBI:mysql:host=cbi2.internal.sanger.ac.uk;dbname=dbSNP_120", "dbsnpro" );

my $dbVar = DBI->connect( "DBI:mysql:host=ecs4.internal.sanger.ac.uk;dbname=arne_variation;port=3352","ensadmin", "ensembl" );


variation_table();




# filling of the variation table from SubSNP and SNP
# creating of a link table variation_id --> subsnp_id
sub variation_table {
  
  $dbVar->do( "ALTER TABLE variation add column snp_id int" );
  $dbVar->do( "ALTER TABLE variation add column subsnp_id int" );

  &::dump( qq{
	      SELECT 1,concat( "rs", snp_id), snp_id
	      FROM SNP
	     }
	 );

  load( "variation", "source_id", "name", "snp_id" );

  $dbVar->do( "ALTER TABLE variation ADD INDEX snpidx( snp_id )" );

  &::dump( qq{
	SELECT subsnp.subsnp_id, subsnplink.snp_id, 
               if( subsnp.validation_status > 0, "VALIDATED", "NOT_VALIDATED" )
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

  $dbVar->do( qq{
		 CREATE TABLE tmp_var2 
                 SELECT tv.subsnp, v.variation_id, tv.valid
		 FROM tmp_var tv, variation v
                 WHERE tv.refsnp = v.snp_id
		});
  $dbVar->do( qq{
		 INSERT INTO variation( source_id, name, parent_variation_id, validation_status, subsnp_id ) 
		 SELECT 1, concat("ss",subsnp), variation_id, valid, subsnp
                 FROM tmp_var2
		} );
}



# successive dumping and loading of tables is typical for this process
# dump does effectively a select into outfile fithout server file system access
sub dump {
  my $sql = shift;
  
  local *FH;

  open ( FH, ">$tmp_dir/tabledump.txt" );

  my $sth = $dbSNP->prepare( $sql, { mysql_use_result => 1 });
  $sth->execute();
  my $first;
  while( my $aref = $sth->fetchrow_arrayref() ) {
    $first = 1;
    for my $col ( @$aref ) {
      if( $first ) {
	$first = 0;
      } else {
	print FH "\t";
      }
      if( defined $col ) {
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

  if( @colnames ) {

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

  

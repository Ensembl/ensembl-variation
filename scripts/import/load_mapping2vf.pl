#!/usr/local/ensembl/bin/perl -w
# this script read mapping data ALL_DIFF into variation_feature table

use strict;
use warnings;
use DBI;
use DBH;
use Getopt::Long;
use Benchmark;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use ImportUtils qw(dumpSQL debug create_and_load load );

my ($TAX_ID, $LIMIT_SQL, $CONTIG_SQL, $TMP_DIR, $TMP_FILE);

my $dbSNP;
my $dbVar;
my $dbCore;
my %rec;

{
  my($dshost, $dsuser, $dspass, $dsport, $dsdbname, # dbSNP db
     $chost, $cuser, $cpass, $cport, $cdbname,      # ensembl core db
     $vhost, $vuser, $vpass, $vport, $vdbname,      # ensembl variation db
     $alldiff_file, $limit);
  
  GetOptions('vuser=s'   => \$vuser,
             'vhost=s'   => \$vhost,
             'vpass=s'   => \$vpass,
             'vport=i'   => \$vport,
             'vdbname=s' => \$vdbname,
	     'chost=s'   => \$chost,
             'cuser=s'   => \$cuser,
             'cpass=s'   => \$cpass,
             'cport=i'   => \$cport,
             'cdbname=s' => \$cdbname,
	     'alldiff=s' => \$alldiff_file,
             'tmpdir=s'  => \$ImportUtils::TMP_DIR,
             'tmpfile=s' => \$ImportUtils::TMP_FILE,
             'limit=i'   => \$limit);
  
  $vhost    ||= 'ecs2';
  $vport    ||= 3361;
  $vuser    ||= 'ensadmin';
  
  $chost    ||= 'ecs2';
  $cuser    ||= 'ensro';
  $cport    ||= 3365;

  usage('-cdbname argument is required.') if(!$cdbname);
  usage('-vdbname argument is required.') if(!$vdbname);
  usage('-alldiff argument is required.') if(!$alldiff_file);

  $TMP_DIR  = $ImportUtils::TMP_DIR;
  $TMP_FILE = $ImportUtils::TMP_FILE;


  $LIMIT_SQL = ($limit) ? " LIMIT $limit " : '';

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


  variation_feature($alldiff_file);
}

sub variation_feature {
  
  my $alldiff_file = shift;

  debug("Dumping Variation");
  
  my $sth = $dbVar->prepare (qq{SELECT variation_id, name
                                FROM variation});
  $sth->execute();

  while(my ($variation_id, $name) = $sth->fetchrow_array()) {
    $rec{$name} = $variation_id;
  }

  $sth->finish();

  open (IN, "$alldiff_file");
  open ( FH, ">$TMP_DIR/$TMP_FILE" );

  while (<IN>) {
    s/^MORE_HITS//;
    my ($ref_id, $slice_name, $ver, $start, $end, $exact, $snp_type, $strand, $score, $ratio) =split;
    next if $ratio <0.5;

    ##make start < end if it is a between type
    if ($exact =~ /between/i) {
      $end = $start-1;
    }

    my ($coord_sys,$assembly,$seq_region_name,$seq_region_start,$seq_region_end,$version) = split /\:/, $slice_name
      if $slice_name =~ /\:/;
    
    my $sth = $dbCore->prepare (qq{SELECT seq_region_id from seq_region where name = ?});
    $sth->execute("$seq_region_name");
    
    my $seq_region_id = $sth->fetchrow_array();

    my $new_seq_region_start = $seq_region_start + $start -1 if ($seq_region_start);
    my $new_seq_region_end = $seq_region_start + $end -1 if ($seq_region_start);
    
    $ref_id = "rs$ref_id";
    print FH "$seq_region_id\t$new_seq_region_start\t$new_seq_region_end\t$strand\t$rec{$ref_id}\t$ref_id\n";
  }
  
  $sth->finish();
  
  close IN;
  close FH;
  
  load($dbVar, "variation_feature","seq_region_id","seq_region_start","seq_region_end",
       "seq_region_strand","variation_id","variation_name");
}

###
###example of usage###
###/creat_variation_feature.pl -cdbname danio_rerio_core_24_4 -vdbname yuan_zfish_snp -alldiff /ecs2/scratch4/yuan/zfish/ds_chNotOndir/ALL_DIFF -vpass ensembl

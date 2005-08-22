#!/usr/local/ensembl/bin/perl -w
# this script read mapping data ALL_DIFF into variation_feature table

use strict;
use warnings;
use DBI;
use DBH;
use Getopt::Long;
use Benchmark;
use Bio::EnsEMBL::DBSQL::DBConnection;
use ImportUtils qw(dumpSQL debug create_and_load load );

my ($TAX_ID, $LIMIT_SQL, $CONTIG_SQL, $TMP_DIR, $TMP_FILE);

my $dbSNP;
my $dbVar;
my $dbCore;

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
	     'mapping_file=s' => \$mapping_file,
             'tmpdir=s'  => \$ImportUtils::TMP_DIR,
             'tmpfile=s' => \$ImportUtils::TMP_FILE,
             'limit=i'   => \$limit);
  
  $vhost    ||= 'ecs2';
  $vport    ||= 3366;
  $vuser    ||= 'ensadmin';
  
  $chost    ||= 'ecs2';
  $cuser    ||= 'ensro';
  $cport    ||= 3364;

  usage('-cdbname argument is required.') if(!$cdbname);
  usage('-vdbname argument is required.') if(!$vdbname);
  usage('-mapping_file argument is required.') if(!$mapping_file);

  $TMP_DIR  = $ImportUtils::TMP_DIR;
  $TMP_FILE = $ImportUtils::TMP_FILE;


  $LIMIT_SQL = ($limit) ? " LIMIT $limit " : '';

  $dbVar = DBH->connect
    ("DBI:mysql:host=$vhost;dbname=$vdbname;port=$vport",$vuser, $vpass,
    {'RaiseError' => 1});
  die("Could not connect to variation database: $!") if(!$dbVar);

  $dbCore = Bio::EnsEMBL::DBSQL::DBConnection->new
    (-host   => $chost,
     -user   => $cuser,
     -pass   => $cpass,
     -port   => $cport,
     -dbname => $cdbname);


  variation_feature($mapping_file);
}

sub variation_feature {
  
  my $mapping_file = shift;

  my (%rec, %source, %status, %rec_pos, %rec_line, %rec_seq_region);

  debug("Dumping Variation");
  
  my $sth = $dbVar->prepare (qq{SELECT variation_id, name, source_id, validation_status
                                  FROM variation});
  $sth->execute();

  while(my ($variation_id, $name, $source_id, $validation_status) = $sth->fetchrow_array()) {
    $rec{$name} = $variation_id;
    $source{$name} = $source_id;
    $status{$name} = defined $validation_status ? $validation_status : '\N';
    #$status{$name} = '\N' if ! defined $validation_status;
  }

  $sth->finish();

  open (IN, "$mapping_file");
  open ( FH, ">$TMP_DIR/$TMP_FILE" );

  while (<IN>) {
    chomp;
    next if /^more |^PARSING/;
    s/^MORE_HITS\s+//;
    my ($ref_id, $slice_name, $start, $end, $strand, $ratio) =split;
    push @{$rec_line{$ref_id}}, $_;
  }
  
  foreach my $key (keys %rec_line) {
    next if @{$rec_line{$key}} >3;
    foreach my $line (@{$rec_line{$key}}) {
      my ($ref_id, $slice_name, $start, $end, $strand, $ratio) =split /\s+/,$line;
      next if $ratio <0.7;

      my ($coord_sys,$assembly,$seq_region_name,$seq_region_start,$seq_region_end,$version);

      ($seq_region_name,$seq_region_start,$seq_region_end) = split /\-/, $slice_name
	if $slice_name =~ /\-/;

      ($coord_sys,$assembly,$seq_region_name,$seq_region_start,$seq_region_end,$version) = split /\:/, $slice_name
	if $slice_name =~ /\:/;

      my $sth;

      if (!$rec_seq_region{$seq_region_name}) {
	$sth = $dbCore->prepare (qq{SELECT seq_region_id from seq_region where name = ?});
	$sth->execute("$seq_region_name");
     
	my $seq_region_id = $sth->fetchrow_array();
	$rec_seq_region{$seq_region_name}=$seq_region_id;
      
	if (!$seq_region_id) {
	  warn "There is no seq_region_id for $ref_id\n";
	  next;
	}
      }
      my $seq_region_id = $rec_seq_region{$seq_region_name};
      my $new_seq_region_start = $seq_region_start + $start -1 if ($seq_region_start);
      my $new_seq_region_end = $seq_region_start + $end -1 if ($seq_region_start);
    
      if (!$rec_pos{$ref_id}{$seq_region_id}{$new_seq_region_start}{$new_seq_region_end}) {
	print FH "$seq_region_id\t$new_seq_region_start\t$new_seq_region_end\t$strand\t$rec{$ref_id}\t$ref_id\t$source{$ref_id}\t",$status{$ref_id} ? $status{$ref_id} : 0,"\n";
	$rec_pos{$ref_id}{$seq_region_id}{$new_seq_region_start}{$new_seq_region_end}=1;
      }
    }
  }
  
  $sth->finish();
  
  close IN;
  close FH;
  
  load($dbVar, "variation_feature","seq_region_id","seq_region_start","seq_region_end",
       "seq_region_strand","variation_id","variation_name","source_id", "validation_status");
}

###
###example of usage###
###./load_mapping2vf.pl -cdbname danio_rerio_core_24_4 -vdbname yuan_zfish_snp -alldiff /ecs2/scratch4/yuan/zfish/ds_chNotOndir/ALL_DIFF -vpass xxxx

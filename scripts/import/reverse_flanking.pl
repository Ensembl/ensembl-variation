#! /usr/local/bin/perl

#

use strict;
#use DBH;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Data::Dumper;
use Bio::SeqIO;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug create_and_load load);

our ($species,$input_dir, $output_dir,$TMP_DIR,$TMP_FILE, $registry_file);

GetOptions('species=s'          => \$species,
	   'tmpdir=s'           => \$ImportUtils::TMP_DIR,
	   'tmpfile=s'          => \$ImportUtils::TMP_FILE,
	  );
$registry_file ||= $Bin . "/ensembl.registry";

#usage('-species argument is required') if(!$species);

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $dbCore = $cdb->dbc;
my $dbVar = $vdb->dbc;

my $sth=$dbVar->prepare(qq{select * from flanking_sequence where seq_region_strand=-1});
$sth->execute();
my ($variation_id,$up_seq,$down_seq,$up_seq_start,$up_seq_end,$down_seq_start,$down_seq_end,$seq_region_id,$seq_region_strand);

open OUT, ">$TMP_DIR/$TMP_FILE" or die "can't open reversed_flank:$!";
$sth->bind_columns(\$variation_id,\$up_seq,\$down_seq,\$up_seq_start,\$up_seq_end,\$down_seq_start,\$down_seq_end,\$seq_region_id,\$seq_region_strand);
while($sth->fetch()) {
  #print "$variation_id,$up_seq,$down_seq,$up_seq_start,$up_seq_end,$down_seq_start,$down_seq_end,$seq_region_id,$seq_region_strand\n";
  if ($up_seq and $down_seq) {
    ($up_seq, $down_seq) = ($down_seq, $up_seq);
    reverse_comp(\$up_seq);
    reverse_comp(\$down_seq);
    $up_seq_start = $up_seq_end = $down_seq_start = $down_seq_end = '\N';
  }
  elsif (! $up_seq and ! $down_seq) {
    my $tmp_seq_start = $up_seq_start;
    my $tmp_seq_end = $up_seq_end;
    ($up_seq_start, $up_seq_end) = ($down_seq_start, $down_seq_end);
    ($down_seq_start, $down_seq_end) = ($tmp_seq_start, $tmp_seq_end);
    $up_seq = $down_seq = '\N';
  }
  elsif ($up_seq and ! $down_seq) {
    $down_seq = $up_seq;
    reverse_comp(\$down_seq);
    $up_seq = '\N';
    ($up_seq_start, $up_seq_end) = ($down_seq_start, $down_seq_end);
    $down_seq_start = '\N';
    $down_seq_end = '\N';
  }
  elsif (! $up_seq and $down_seq) {
    $up_seq = $down_seq;
    reverse_comp(\$up_seq);
    $down_seq = '\N';
    ($down_seq_start, $down_seq_end) = ($up_seq_start, $up_seq_end);
    $up_seq_start = '\N';
    $up_seq_end = '\N';
  }
  $seq_region_strand = 1 if ($seq_region_strand = -1);
  print OUT "$variation_id\t$up_seq\t$down_seq\t$up_seq_start\t$up_seq_end\t$down_seq_start\t$down_seq_end\t$seq_region_id\t$seq_region_strand\n";
}

$sth->finish;

$dbVar->do(qq{CREATE TABLE reversed_flank like flanking_sequence});
load($dbVar,"reversed_flank","variation_id","up_seq","down_seq","up_seq_region_start","up_seq_region_end","down_seq_region_start","down_seq_region_end","seq_region_id","seq_region_strand");

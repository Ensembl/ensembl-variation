#!/usr/bin/env perl
#

use strict;
Use DBH;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use Bio::SeqIO;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug create_and_load load);

our ($species, $seq_region_id, $chr_name, $TMP_DIR, $TMP_FILE);

GetOptions('species=s'    => \$species,
	   'seq_region_id=i' => \$seq_region_id,
	   'tmpdir=s'     => \$ImportUtils::TMP_DIR,
	   'tmpfile=s'    => \$ImportUtils::TMP_FILE,
	  );
my $registry_file ||= $Bin . "/ensembl.registry";
$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;
$species ||= 'mouse';

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
print "core_db is ",$cdb->dbc->dbname, " and var_db is ",$vdb->dbc->dbname," and seq_region_id is $seq_region_id\n";

my $dbCore = $cdb->dbc->db_handle;
my $dbVar = $vdb->dbc->db_handle;
my %rec_repeat;
my $slice_adaptor = $cdb->get_SliceAdaptor();

&run();

sub run {

  my (%rec_name,%rec_start);
  my ($variation_name,$seq_start);
  my $sthv = $dbVar->prepare(qq{select variation_name,seq_region_start
                                from variation_feature
                                where seq_region_id = $seq_region_id
                                and source_id=2
  			       });
  $sthv->execute();
  $sthv->bind_columns(\$variation_name,\$seq_start);
  my $count;
  while($sthv->fetch()) {
    $count++;
    $rec_name{$variation_name}=$seq_start;
  }
  
  my $slice = $slice_adaptor->fetch_by_seq_region_id($seq_region_id);
  my $slice_size = 10000;
  my $start = 1;
  my $end = $start + $slice_size -1;
  $end = $slice->length if($end > $slice->length);
  while ($start <= $slice->length) {
    my $sub_slice = $slice->sub_Slice($start,$end);
    my @variation_features = @{$sub_slice->get_all_VariationFeatures};
    my @vfs;
    foreach my $vf (@variation_features) {
      if ($rec_name{$vf->variation_name}) {
	push @vfs, $vf;
      }
    }
    my %rec_starts = map {$_->start, $_->variation_name}  @vfs;
    
    my $repeat_masked_slice = $sub_slice->get_repeatmasked_seq();
    my $dna = $repeat_masked_slice->seq();
    
    my @dnas = split "", $dna;
    #warn `ps -p $$ -o vsz |tail -1`;
    foreach my $var_start (sort {$a<=>$b} keys %rec_starts) {
      if ($dnas[$var_start-1] eq 'N') {
	my $variation_name = $rec_starts{$var_start};
	$rec_repeat{$variation_name}=1;
      }
    }
    $start = $end + 1;
    $end = $start + $slice_size -1;
    $end = $slice->length if($end > $slice->length);
  }
}


open OUT, ">$TMP_DIR/$TMP_FILE\_$seq_region_id";
foreach my $variation_name (keys %rec_repeat) {
  print OUT "$variation_name\n";
}

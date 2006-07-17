#! /usr/local/ensembl/bin/perl
#
#./bsub_ssaha2.pl -species human -input_dir /ecs4/scratch2/yuan/hum/mapping_36 -output_dir /ecs4/scratch2/yuan/hum/mapping_36 -target_dir /ecs4/scratch2/yuan/hum/mapping_36/target_dna -start 1 -end 9

use strict;
#use DBH;
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

# try to use named options here and write a sub usage() function
# eg -host -user -pass -port -snp_dbname -core_dbname etc
# optional chromosome name or genomic sequence file
# optional more than one genomic sequence file
# optional a directory or sequence files (for unknown placing)


our ($species, $output_dir, $input_dir, $target_dir, $TMP_DIR, $TMP_FILE, $start, $end);

GetOptions('species=s'    => \$species,
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

my $dbCore = $cdb->dbc->db_handle;
my $dbVar = $vdb->dbc->db_handle;

$start ||=1;

&run();

sub run {

  my $slice_adaptor = $cdb->get_SliceAdaptor();

  my (%seq_region_names);

  #my $sthc = $dbCore->prepare(qq{select sr.seq_region_id from seq_region_attrib sra, attrib_type at, seq_region sr where sra.attrib_type_id=at.attrib_type_id and at.code="toplevel" and sr.seq_region_id = sra.seq_region_id });
  my $sthc = $dbCore->prepare(qq{select sr.seq_region_id from seq_region sr, coord_system cs where cs.coord_system_id=sr.coord_system_id and cs.name="chromosome"});
  $sthc->execute();
  while (my ($seq_region_id) = $sthc->fetchrow_array()) {
    my $call = "bsub -q normal -o mapping_out /usr/local/bin/perl map_variation_feature.pl -species $species -seq_region_id $seq_region_id -tmpdir $TMP_DIR -tmpfile $TMP_FILE";
    system($call);
    print "submitting job for $seq_region_id\n";
  }
}


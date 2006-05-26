#! /usr/local/bin/perl
#

use strict;
#use DBH;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use Bio::SeqIO;
use FindBin qw( $Bin );
use Getopt::Long;
use ImportUtils qw(dumpSQL debug create_and_load load);
use DBI qw(:sql_types);
our ($species, $input_file, $seq_region_id, $TMP_DIR, $TMP_FILE);

GetOptions('species=s'         => \$species,
	   'seq_region_id=i'   => \$seq_region_id,
           'tmpdir=s'          => \$ImportUtils::TMP_DIR,
           'tmpfile=s'         => \$ImportUtils::TMP_FILE,
          );
my $registry_file ||= $Bin . "/ensembl.registry";

usage('-species argument is required') if(!$species);

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;
$TMP_FILE .= "\.$seq_region_id";
Bio::EnsEMBL::Registry->load_all( $registry_file );
my $cdb1 = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core1');
my $cdb2 = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core2');
my $vdb1 = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation1');
my $vdb2 = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation2');
my $dbCore1 = $cdb1->dbc->db_handle;
my $dbCore2 = $cdb2->dbc->db_handle;
my $dbVar1 = $vdb1->dbc->db_handle;
my $dbVar2 = $vdb2->dbc->db_handle;

my $sa1 = $cdb1->get_SliceAdaptor();
my $sa2 = $cdb2->get_SliceAdaptor();

open OUT, ">$TMP_DIR/$TMP_FILE" or die "can't open tmpfile:$!";
my ($seq_region_start,$seq_region_end,$seq_region_strand,$variation_id,$allele_string,$variation_name,$map_weight,$flags,$source_id,$validation_status);

print "seq_region_id is $seq_region_id\n";
my $sthv1 = $dbVar1->prepare(qq{select seq_region_start,seq_region_end,seq_region_strand,variation_id,allele_string,variation_name,map_weight,flags,source_id,validation_status
                                from variation_feature
                                where seq_region_id = ?
                                #limit 10
                                }
                              );
$sthv1->bind_param(1,$seq_region_id,SQL_INTEGER);
$sthv1->execute();
$sthv1->bind_columns(\$seq_region_start,\$seq_region_end,\$seq_region_strand,\$variation_id,\$allele_string,\$variation_name,\$map_weight,\$flags,\$source_id,\$validation_status);

while ($sthv1->fetch()) {
  my $old_slice = $sa1->fetch_by_seq_region_id($seq_region_id);
  #print "seq_region_id is $seq_region_id and old_slice is $old_slice\n";
  my $slice_newdb_oldasm = $sa2->fetch_by_region('chromosome', $old_slice->seq_region_name, undef, undef, undef, 'NCBIm35');
  my $feat = new Bio::EnsEMBL::SimpleFeature (
					      -START => $seq_region_start,
					      -END => $seq_region_end,
					      -STRAND => $seq_region_strand,
					      -SLICE => $slice_newdb_oldasm,
					     );

  #print $feat->start," ",$feat->strand," ",$feat->slice,"\n";
  my @segments;
  eval {
    @segments = @{$feat->feature_Slice->project('chromosome','NCBIM36')};
  };
  if ($@) {
    print "name is ",$feat->feature_Slice->name,"\n";
  }
  my @slices_newdb_newasm = map { $_->to_Slice } @segments;

  # new feature spans start of first projection segment to end of last segment
  #print "old_seq_region_name is ",$old_slice->get_seq_region_id," old_start is ",$feat->start," old_end is ",$feat->end," new_seq_region_name is ",$slices_newdb_newasm[0]->get_seq_region_id," new_start is ",$slices_newdb_newasm[0]->start," new_end is ",$slices_newdb_newasm[-1]->end,"\n" if (@slices_newdb_newasm);
  if (@slices_newdb_newasm) {
    print OUT join "\t",$slices_newdb_newasm[0]->get_seq_region_id,$slices_newdb_newasm[0]->start,$slices_newdb_newasm[-1]->end,$slices_newdb_newasm[0]->strand,"$variation_id\t$allele_string\t$variation_name\t$map_weight\t$flags\t$source_id\t$validation_status\n";
  }
  else {
    print "there is no mapping for $seq_region_id\t$seq_region_start\t$seq_region_end\n";
  }
}

load ($dbVar2,"variation_feature","seq_region_id","seq_region_start","seq_region_end","seq_region_strand","variation_id","allele_string","variation_name","map_weight","flags","source_id","validation_status");

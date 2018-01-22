#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

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
our ($species, $input_file, $seq_region_id, $source_name, $TMP_DIR, $TMP_FILE);

GetOptions('species=s'         => \$species,
	   'seq_region_id=i'   => \$seq_region_id,
	   'source_name=s'     => \$source_name,
           'tmpdir=s'          => \$ImportUtils::TMP_DIR,
           'tmpfile=s'         => \$ImportUtils::TMP_FILE,
          );
my $registry_file ||= $Bin . "/ensembl.registry";

usage('-species argument is required') if(!$species);

$ImportUtils::TMP_FILE .= "_$seq_region_id";
$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

Bio::EnsEMBL::Registry->load_all( $registry_file );
my $cdb1 = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core1');
my $cdb2 = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core2');
my $vdb1 = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation2');
my $vdb2 = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation2');
my $dbCore1 = $cdb1->dbc->db_handle;
my $dbCore2 = $cdb2->dbc->db_handle;
my $dbVar1 = $vdb1->dbc->db_handle;
my $dbVar2 = $vdb2->dbc->db_handle;

my $sa1 = $cdb1->get_SliceAdaptor();
my $sa2 = $cdb2->get_SliceAdaptor();
my $buffer = {};

#open OUT, ">$TMP_DIR/$TMP_FILE\_$seq_region_id" or die "can't open tmpfile:$!";
#open ERR, ">$TMP_DIR/ERR\_$seq_region_id" or die "can't open errfile:$!";

my ($seq_region_start,$seq_region_end,$seq_region_strand,$variation_id,$allele_string,$variation_name,$map_weight,$flags,$source_id,$validation_status);

print "seq_region_id is $seq_region_id\n";
my $sthv1 = $dbVar1->prepare(qq{select seq_region_start,seq_region_end,seq_region_strand,variation_id,allele_string,variation_name,map_weight,flags,source_id,validation_status
                                from variation_feature
                                #from variation_feature_venter_map1_no_hap
                                where seq_region_id = ?
                                #and variation_id=7073660
                                #limit 10
                                }
                              );
$sthv1->bind_param(1,$seq_region_id,SQL_INTEGER);
$sthv1->execute();
$sthv1->bind_columns(\$seq_region_start,\$seq_region_end,\$seq_region_strand,\$variation_id,\$allele_string,\$variation_name,\$map_weight,\$flags,\$source_id,\$validation_status);

while ($sthv1->fetch()) {
  my $old_slice = $sa1->fetch_by_seq_region_id($seq_region_id);
  #print "seq_region_id is $seq_region_id and old_slice is $old_slice\n";
  my $slice_newdb_oldasm = $sa2->fetch_by_region('chromosome', $old_slice->seq_region_name, undef, undef, undef, 'NCBI36');
  if (!$slice_newdb_oldasm) {
    print_buffered($buffer,"$TMP_DIR/ERR\_$seq_region_id\_$source_name",join ("\t",$variation_id,$seq_region_id,$seq_region_start,$seq_region_end) . "\n");
    next;
  }

  my $indel;
  #for indels, change start and end
  if ($seq_region_start == $seq_region_end+1) {
    ($seq_region_start,$seq_region_end) = ($seq_region_end,$seq_region_start);
    $indel = 1;
  }

  my $feat = new Bio::EnsEMBL::SimpleFeature (
					      -START => $seq_region_start,
					      -END => $seq_region_end,
					      -STRAND => $seq_region_strand,
					      -SLICE => $slice_newdb_oldasm,
					     );

  #print $feat->start," ",$feat->strand," ",$feat->slice,"\n";
  my @segments;

  eval {
    @segments = @{$feat->feature_Slice->project('chromosome','GRCh37')};
  };
  if ($@) {
    print "no segments found for $variation_id\n" if $@;
    print_buffered($buffer,"$TMP_DIR/ERR\_$seq_region_id\_$source_name",join ("\t",$variation_id,$seq_region_id,$seq_region_start,$seq_region_end) . "\n");
    next;
  }  

  my @slices_newdb_newasm = map { $_->to_Slice }  @segments;
  @slices_newdb_newasm = sort {$a->start<=>$b->start} @slices_newdb_newasm;
  
  $validation_status ||= '\N';

  #my $new_source_ref = $dbVar2->selectall_arrayref(qq{SELECT source_id FROM source WHERE name = "$source_name"});
  #my $new_source_id = ($new_source_ref) ? $new_source_ref->[0][0] : $source_id;
  my $new_source_id = $source_id;#for celera/venter/watson dbsnp variations

  #make new strand is same as old strand 
  #my $new_strand = $seq_region_strand;
  my ($new_start,$new_end);

  if (@slices_newdb_newasm) {
    if ($indel) {
      $new_start = $slices_newdb_newasm[0]->start +1;
      $new_end = $slices_newdb_newasm[-1]->end -1;
    }
    else {
      $new_start = $slices_newdb_newasm[0]->start;
      $new_end = $slices_newdb_newasm[-1]->end;
    }
   
    my $new_strand = $slices_newdb_newasm[0]->strand;
    # new feature spans start of first projection segment to end of last segment
    #print "old_seq_region_name is ",$old_slice->get_seq_region_id," old_start is ",$feat->start," old_end is ",$feat->end," new_seq_region_name is ",$slices_newdb_newasm[0]->get_seq_region_id," new_start is ",$new_start," new_end is ",$new_end,"\n" if (@slices_newdb_newasm);

    print_buffered($buffer,"$TMP_DIR/$TMP_FILE\_$seq_region_id\_$source_name",join "\t",$slices_newdb_newasm[0]->get_seq_region_id,$new_start,$new_end,$new_strand,"$variation_id\t$allele_string\t$variation_name\t$map_weight\t$flags\t$new_source_id\t$validation_status\n");
    #print OUT join "\t",$slices_newdb_newasm[0]->get_seq_region_id,$new_start,$new_end,$slices_newdb_newasm[0]->strand,"$variation_id\t$allele_string\t$variation_name\t$map_weight\t$flags\t$source_id\t$validation_status\n";
  }
  else {
    print_buffered($buffer,"$TMP_DIR/ERR\_$seq_region_id\_$source_name",join ("\t",$variation_id,$seq_region_id,$seq_region_start,$seq_region_end) . "\n");
    #print ERR "$variation_id\t$seq_region_id\t$seq_region_start\t$seq_region_end\n";
  }
}

print_buffered($buffer);

debug("Loading mapping data...");
if (-e "$TMP_DIR/$TMP_FILE\_$seq_region_id\_$source_name") {
  system("mv $TMP_DIR/$TMP_FILE\_$seq_region_id\_$source_name  $TMP_DIR/$TMP_FILE") ;
  my $mapping_ref = $dbVar2->selectall_arrayref(qq{show tables like "variation_feature_GRCh37_MAPPING_PATRIC\_$source_name"});
  if (!$mapping_ref->[0][0]) {
    create_and_load ($dbVar2,"variation_feature_GRCh37_MAPPING_PATRIC\_$source_name","seq_region_id i*","seq_region_start i","seq_region_end i","seq_region_strand i","variation_id i*","allele_string","variation_name","map_weight","flags","source_id i","validation_status");
  }
  else {
    load ($dbVar2,"variation_feature_GRCh37_MAPPING_PATRIC\_$source_name","seq_region_id","seq_region_start","seq_region_end","seq_region_strand","variation_id","allele_string","variation_name","map_weight","flags","source_id","validation_status");
  }
}

debug("Loading ERROR Mapping data...");
if (-e "$TMP_DIR/ERR\_$seq_region_id\_$source_name") {
  system("mv $TMP_DIR/ERR\_$seq_region_id\_$source_name  $TMP_DIR/$TMP_FILE");
  my $failed_ref = $dbVar2->selectall_arrayref(qq{show tables like "failed_mapping_vf\_$source_name"});
  if (!$failed_ref->[0][0]) {
    create_and_load($dbVar2,"failed_mapping_vf\_$source_name","variation_id i*","seq_region_id i","seq_region_start i","seq_region_end i");
  }
  else {
    load($dbVar2,"failed_mapping_vf\_$source_name","variation_id","seq_region_id","seq_region_start","seq_region_end");
  }
}

sub print_buffered {
    my $buffer = shift;
    my $filename = shift;
    my $text = shift;

    local *FH;

    if( ! $filename ) {
	# flush the buffer
	foreach my $file (keys %{$buffer}){
	    open( FH, ">>$file" ) or die;
	    print FH $buffer->{ $file };
	    close FH;
	}
	%{$buffer} = ();

    } else {
	$buffer->{ $filename } .= $text;
	if( length( $buffer->{ $filename } ) > 10_000 ) {
	    open( FH, ">>$filename" ) or die;
	    print FH $buffer->{ $filename };
	    close FH;
	    $buffer->{ $filename } = '';
	}
    }
}

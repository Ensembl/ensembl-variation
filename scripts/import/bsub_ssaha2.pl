#!/usr/bin/env perl
# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
#./bsub_ssaha2.pl -species human -input_dir /ecs4/scratch2/yuan/hum/mapping_36 -output_dir /ecs4/scratch2/yuan/hum/mapping_36 -target_dir /ecs4/scratch2/yuan/hum/mapping_36/target_dna -start 1 -end 9 -map_all

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


our ($species, $map_by_chr, $map_all, $split, $rerun, $output_dir, $input_dir, $target_file, $run, $parse, $TMP_DIR, $TMP_FILE, $start, $end);

GetOptions('species=s'    => \$species,
	   'output_dir=s' => \$output_dir,
	   'input_dir=s'  => \$input_dir,
           'target_file=s' => \$target_file,
           'start=i'      => \$start,
           'end=i'        => \$end,
	   'run'          => \$run,
	   'parse'        => \$parse,
	   'split'        => \$split, #option to split reads in groups to reduce run time
	   'rerun=i'      => \$rerun, #to rerun it if some jobs failed i is the split start number
	   'map_by_chr=i' => \$map_by_chr,
	   'map_all'      => \$map_all,
	   'tmpdir=s'     => \$ImportUtils::TMP_DIR,
	   'tmpfile=s'    => \$ImportUtils::TMP_FILE,
	  );
my $registry_file ||= $Bin . "/ensembl.registry";

my $run_parse;
$run_parse = "-run" if $run;
$run_parse = "-parse" if $parse;

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $dbCore = $cdb->dbc->db_handle;
#my $dbVar = $vdb->dbc->db_handle;

#$start ||=1;

&run();
#get_read_flank_seq($cdb, $vdb, $read_file);

sub run {

  my $slice_adaptor = $cdb->get_SliceAdaptor();

  my (%seq_region_names);

  my $sthc = $dbCore->prepare(qq{select sr.name from seq_region_attrib sra, attrib_type at, seq_region sr where sra.attrib_type_id=at.attrib_type_id and at.code="toplevel" and sr.seq_region_id = sra.seq_region_id });
  #my $sthc = $dbCore->prepare(qq{select sr.name from seq_region sr, coord_system cs where sr.coord_system_id = cs.coord_system_id and cs.name='chromosome'});
  $sthc->execute();

  while (my ($seq_region_name) = $sthc->fetchrow_array()) {
    if ($seq_region_name =~ /^\d+|NT|^A|^B|^C|scaffold|random|^M|^Zv/g) {
      $seq_region_name  = "random";
    }
    $seq_region_names{$seq_region_name}=1;
  }

  if ($map_by_chr) {
    foreach my $seq_region_name (keys %seq_region_names) {
      my $input_dir_new = "$input_dir/$seq_region_name";
      system(mkdir "$output_dir/$seq_region_name");
      my $output_dir_new = "$output_dir/$seq_region_name";
      print "input_dir is $input_dir_new and output_dir is $output_dir_new\n";
      my $num_seq = `ls $input_dir_new/*seq |wc`;
      ($num_seq) = $num_seq =~ /\s+(\d+)\s+/;
      print "number of seq in $input_dir is #$num_seq#\n";
      next if ($num_seq ==1 and -z "$input_dir_new/1_query_seq");
      my $call = "./run_ssaha2.pl -start $start -end $num_seq -chr $seq_region_name -input_dir $input_dir_new -output_dir $output_dir_new -target_file $target_file";
      system($call);
      print "submit job for $seq_region_name\n";
    }
  }
  elsif ($map_all) {
    print "input_dir is $input_dir and output_dir is $output_dir\n";
    my $call = "./run_ssaha2.pl ";
    $call .= ($start) ? "-start $start -end $end" : '' ;
    $call .= " -split $split " if $split ;
    $call .= " -rerun $rerun" if $rerun;
    $call .= " -input_dir $input_dir -output_dir $output_dir -target_file $target_file $run_parse";
    print "call is $call\n";
    system($call);
    print "submit job for $start and $end\n";
  }
}


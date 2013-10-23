#!/usr/bin/env perl
# Copyright 2013 Ensembl
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
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

#
#./bsub_mapping.pl -species mouse -job flank -tmpdir [tmp_dir] -tmpfile patric_mapping

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


our ($species, $job, $source_name, $output_dir, $input_dir, $target_dir, $TMP_DIR, $TMP_FILE, $start, $end);

GetOptions('species=s'    => \$species,
	   'job=s'        => \$job,
	   'source_name=s'=> \$source_name,
           'tmpdir=s'     => \$ImportUtils::TMP_DIR,
	   'tmpfile=s'    => \$ImportUtils::TMP_FILE,
	  );
my $registry_file ||= $Bin . "/ensembl.registry";

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

my $script_name;
if ($job =~/flank/) {
  $script_name = "map_flanking_sequence.pl";
}
elsif($job =~ /var/) {
  $script_name = "map_variation_feature.pl";
}


Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core1');
#my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');

my $dbCore = $cdb->dbc->db_handle;
#my $dbVar = $vdb->dbc->db_handle;

$start ||=1;

run();

sub run {

  my $slice_adaptor = $cdb->get_SliceAdaptor();

  my ($count,%seq_region_names);

  #my $sthc = $dbCore->prepare(qq{select sr.seq_region_id from seq_region_attrib sra, attrib_type at, seq_region sr where sra.attrib_type_id=at.attrib_type_id and at.code="toplevel" and sr.seq_region_id = sra.seq_region_id });
  my $sthc = $dbCore->prepare(qq{select sr.seq_region_id from seq_region_attrib sra, attrib_type at, seq_region sr where sra.attrib_type_id=at.attrib_type_id and at.code="toplevel" and sr.seq_region_id = sra.seq_region_id});
  $sthc->execute();
  while (my ($seq_region_id) = $sthc->fetchrow_array()) {
    my $call = "bsub -q normal -o $TMP_DIR/mapping_out\_$seq_region_id\_$job /usr/bin/env perl $script_name -species $species -seq_region_id $seq_region_id -source_name $source_name -tmpdir $TMP_DIR -tmpfile $TMP_FILE";
    system($call);
    $count++;
    print "submitting job for $seq_region_id count $count\n";
  }
}


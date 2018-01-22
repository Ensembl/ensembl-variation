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
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(expand reverse_comp);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Data::Dumper;
use Bio::SeqIO;
use FindBin qw( $Bin );
use Getopt::Long;
#use ImportUtils qw(dumpSQL debug create_and_load load);

# try to use named options here and write a sub usage() function
# eg -host -user -pass -port -snp_dbname -core_dbname etc
# optional chromosome name or genomic sequence file
# optional more than one genomic sequence file
# optional a directory or sequence files (for unknown placing)


our ($species, $target_dir, $start,$end, $reverse_comp, $file,$seq_region_id,$coord_sys_name, $seq_region_name, $TMP_DIR, $TMP_FILE, $repeat_masked, $chr_name,$read_file, $registry_file);

GetOptions('species=s'         => \$species,
	   'target_dir=s'      => \$target_dir,
	   'repeat_masked'        => \$repeat_masked,
	   'start=i'           => \$start,
	   'end=i'             => \$end,
	   'reverse_comp'      => \$reverse_comp,
	   'seq_region_id=i' => \$seq_region_id,
	   'file=s'            => \$file,
	   'seq_region_name=s' => \$seq_region_name,
	   'coord_sys_name=s'  => \$coord_sys_name,
        'registry_file=s' => \$registry_file,		
	  );

$registry_file ||= $Bin . "/ensembl.registry";

$coord_sys_name  ||= 'toplevel';

die('-species argument is required') if(!$species);
#usage('-target_dir argument is required') if(!$target_dir);

#$TMP_DIR  = $ImportUtils::TMP_DIR;
#$TMP_FILE = $ImportUtils::TMP_FILE;

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $dbCore = $cdb->dbc->db_handle;


#&get_chrom($cdb,$start,$end);
#&get_contigs_old_schema($dba);
#&get_olap_slice($cdb);
&fetch_all($cdb,$coord_sys_name,$target_dir);
#reverse_file($file);

sub get_chrom {

  if (!$coord_sys_name and !$seq_region_name) {
    print STDERR "$coord_sys_name or $seq_region_name are undefined\n";
    $coord_sys_name = "toplevel";
  }

  my ($dba) = @_;

  my $sa = $dba->get_SliceAdaptor();
  my $slice;
  if ($seq_region_id) {
    $slice = $sa->fetch_by_seq_region_id($seq_region_id,$start,$end);
  }
  elsif ($seq_region_name) {
    $slice = $sa->fetch_by_region( $coord_sys_name, $seq_region_name,$start,$end);
  }
  my $output_name = "$target_dir/".$slice->seq_region_name().'-'.$slice->start.'.fa';
  $output_name .= '_rev' if $reverse_comp;
  #open OUT, ">$target_dir/".$slice->seq_region_name().'-'.$slice->start.'_rev.fa';
  open OUT, ">$output_name";
  print OUT ">".$slice->seq_region_name().'-'.$slice->start().'-'.$slice->end()."\n";
  my $seq = $slice->seq;
  reverse_comp(\$seq) if ($seq and $reverse_comp);
  fasta($seq);

}

sub reverse_file {
  my $file = shift;
  my $seq;

  open IN, "$target_dir/$file";
  open OUT, ">$target_dir/$file\_rev";

  while (<IN>) {
    if (/^>/) {
      print OUT;
    }
    else {
      chomp;
      $seq .= $_;
    }
  }
  reverse_comp(\$seq);
  fasta($seq);
}

sub get_contigs_old_schema {

  my $dba = shift;
  my $sth = $dba->prepare(qq{
			     SELECT a.contig_id
			     FROM assembly a
			     , contig g
			     , clone c
			     WHERE a.contig_id = g.contig_id
			     AND g.clone_id = c.clone_id
			     #AND a.chromosome_id = \"$chr_name\"
			     AND a.type = 'RGSC3.1'
			    });
  $sth->execute;

  my $out = Bio::SeqIO->new(
			    -FORMAT => 'fasta',
			   );

  my $cla = $dba->get_RawContigAdaptor;
  while (my ($cid) = $sth->fetchrow) {
    my $contig = $cla->fetch_by_dbID($cid);
    $out->write_seq($contig);
  }
  $sth->finish;

}

sub get_olap_slice {

  ####change to use ensembl-4 api
  my $dba = shift;
  my $slice_adaptor = $dba->get_SliceAdaptor();
  #my $slice = $slice_adaptor->fetch_by_chr_name($chr_name);
  my $slice = $slice_adaptor->fetch_by_region ('chromosome',4,68800000,69500000);
  my $SIZE=10000;
  my $OVERLAP = 0;

  my $start = 1;
  my $end   = $start + $SIZE -1;
  $end = $slice->length if($end > $slice->length);

  #open OUT, ">$target_dir/".$slice->seq_region_name().'_olap.fa' or die "can't open fasta files $! \n";

  while($start+$OVERLAP <= $slice->length) {
    open OUT, ">$target_dir/olap_".$start.'.fa' or die "can't open fasta files $! \n";
    print OUT ">". ($slice->start() + $start - 1) . '-' . ($slice->start() + $end - 1) . "\n";
    my $seq = $slice->subseq($start, $end);

    fasta($seq);

    $start = $end - $OVERLAP + 1;
    $end = $start + $SIZE -1;
    $end = $slice->length if($end > $slice->length);
  }
}

sub fetch_all {

  my ($dba,$coord_sys_name,$target_dir) = @_;

  if (!$coord_sys_name) {
    print STDERR "$coord_sys_name are undefined using toplevel instead\n";
    $coord_sys_name = "toplevel";
  }
  
  my $slice_adaptor = $dba->get_SliceAdaptor();
  my @names;
  
  my @slices = @{$slice_adaptor->fetch_all($coord_sys_name,undef,1,1)};
  print "repeat_masked is $repeat_masked\n";
  if (@slices) {
    foreach my $slice (@slices) {
      if ($repeat_masked) { 
	$slice = $slice->get_repeatmasked_seq();
      }
      if ($slice->seq_region_name() =~ /NT|^A|^B|^C|^c|^Zv7/) {
      #print "seq_region_name is ",$slice->seq_region_name(),"\n";
      #if ($slice->seq_region_name() =~ /random|M|NT|^Contig/) {
	#print "name is ",$slice->seq_region_name(),"\n";
	open OUT, ">>$target_dir/random.fa" or die "can't open NT_files $!\n ";
      }
      else {
	open OUT, ">$target_dir/".$slice->seq_region_name().'.fa' or die "can't open fasta files $! \n";
	#open OUT, ">>$target_dir/contig.fa" or die "can't open fasta files $! \n";
      }
      print OUT ">".$slice->seq_region_name().'-'.$slice->start().'-'.$slice->end()."\n";
      my $seq = $slice->seq;
      fasta($seq);
    }
  }
}



sub fasta {
  my $seq = shift;

  my $len = length($seq);
  my $start = 0;
  while($start < $len) {
    print OUT substr($seq, $start, 80), "\n";
    $start += 80;
  }

}

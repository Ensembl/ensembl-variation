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

#example: ./fetch_base_qual.pl -i [index_file] -r [reads_file] 19866807248411 -1 11

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::Index::Fastq;
use Bio::EnsEMBL::Utils::Sequence qw(expand reverse_comp);
use Getopt::Long;

our ($index_file, $reads_file);

GetOptions(
    'index_file=s' => \$index_file,
	'reads_file=s'  => \$reads_file
);

my ($query_name,$query_strand,$snp_pos,$new) = @ARGV;
my $flank;
$flank = get_flanking_seq($query_name,$query_strand,$snp_pos) if !$new;
$flank = get_flanking_seq1($query_name,$query_strand,$snp_pos) if $new;

my ($flank_5_qual,$flank_3_qual,$snp_base,$snp_qual) = @$flank;
print "flanking_5_qual is @$flank_5_qual and flanking_3_qual is @$flank_3_qual, snp_base is $snp_base and snp_qual is $snp_qual\n";


sub get_pfetch_sequence_and_quality {
    my $query_name = shift;
    
    my $fastq_index = Bio::Index::Fastq->new(-filename => $index_file);

    my $seq = $fastq_index->fetch($query_name);

    return $seq;
}

sub get_flanking_seq1 {
  my ($query_name,$query_strand,$snp_pos) = @_;
  my ($qual_5,$qual_3,$snp_base,$snp_qual,$start,$end);
  
  my ($REC_SEQ,$REC_QUAL) = make_seq_qual_hash($reads_file);
  my %REC_SEQ = %$REC_SEQ;
  my %REC_QUAL = %$REC_QUAL;
  
  die "hash %REC_SEQ is empty" unless (keys %REC_SEQ >1);

  my $seq = $REC_SEQ{"$query_name"};
  my $qual = $REC_QUAL{"$query_name"};
  print "$seq\n";
  print "$qual\n";
  my $len = length($seq);

  $start = $snp_pos-5;
  $end = $snp_pos+5;
  $start = 1 if $snp_pos -5 <= 0;
  $end = $len if $snp_pos+5 > $len;

  $snp_base = substr($seq,$snp_pos-1,1);
  my $flanking_bases = substr($seq,$start-1,11);
  print "old_flanking_bases are $flanking_bases\n";
  $snp_qual = substr($qual,$snp_pos-1,1);
  $qual_5 = substr($qual,$start-1,5);
  $qual_3 = substr($qual,$snp_pos,5);
  if ($query_strand ==-1) {
    $snp_base =~ tr/ACGTacgt/TGCAtgca/;
    ($qual_5,$qual_3) = ($qual_3,$qual_5);
    print "new_snp_base is $snp_base\n";
  }
  my $flanking_quals = substr($qual,$start-1,11);
  print "unpacked qual is $flanking_quals\n";
  my @flanking_quals = split '',$flanking_quals;
  @flanking_quals = map{unpack("C",$_)-33} @flanking_quals;
  print "old_flanking_quals are @flanking_quals\n";
  my @qual_5 = split '',$qual_5;
  my @qual_3 = split '',$qual_3;
    
  if ($query_strand==-1) {
    @qual_5 = reverse @qual_5;
    @qual_3 = reverse @qual_3;
  }
  @qual_5 = map{unpack("C",$_)-33} @qual_5;
  @qual_3 = map{unpack("C",$_)-33} @qual_3;
  $snp_qual = unpack("C",$snp_qual)-33;

  $qual_5 = [ @qual_5 ];
  $qual_3 = [ @qual_3 ];
                                  
  return ([$qual_5,$qual_3,$snp_base,$snp_qual]);
}

sub make_seq_qual_hash {
  my $reads_file = shift;
  my (%REC_SEQ,%REC_QUAL,$name);
  #print "file is $reads_file\n";
  open FASTQ, $reads_file or die "Can't open $reads_file";
  while (<FASTQ>) {
    chomp;
    if (/^\@(.*)/) {
      $name = $1;
      ($name) = split /\s+/,$name;
    }
    elsif (/^\+(.*)/) {
      $name = $1;
      ($name) = split /\s+/,$name;
    }
    elsif ($name and /^(\!.*)$/) {
      $REC_QUAL{$name} = $1;
      undef $name;
    }
    elsif ($name and /[^?%*]/) {
      $REC_SEQ{$name} = $_;
      undef $name;
    }
  }
  close FASTQ;
  return (\%REC_SEQ,\%REC_QUAL);
}
                                                              

sub get_flanking_seq {
  my ($query_name,$query_strand,$snp_pos) = @_;
  my ($qual_5,$qual_3,$snp_base,$snp_qual,$start,$end);

  #my $t0 = [gettimeofday];

  local($\) = undef;
  my $seqobj = &get_pfetch_sequence_and_quality($query_name);
  my ($seqstr, @qualarray);
  my $len = $seqobj->length;
  print $seqobj->seq,"\n";
  print @{$seqobj->qual},"\n";
  print "old_base ",substr($seqobj->seq,$snp_pos-6,11),"\n";
  print "old_qual ",@{$seqobj->subqual($snp_pos-5,$snp_pos+5)},"\n";

  $start = $snp_pos-5;
  $end = $snp_pos+5;
  $start = 1 if $snp_pos -5 < 0;
  $end = $len if $snp_pos+5 > $len;

  $seqstr = $seqobj->baseat($snp_pos);
  $snp_qual = $seqobj->qualat($snp_pos);
  $qual_5 = $seqobj->subqual($start,$snp_pos-1);
  $qual_3 = $seqobj->subqual($snp_pos+1,$end);

  if ($query_strand ==-1) {
    $seqstr =~ tr/ACTGactg/TGACtgac/;
    ($qual_5,$qual_3) = ($qual_3,$qual_5);
    my @qual_5 = reverse @$qual_5;
    my @qual_3 = reverse @$qual_3;
    #print "qual_5 is @qual_5 and qual_3 is @qual_3\n";
    $qual_5 = [@qual_5];
    $qual_3 = [@qual_3];
  } 
                        
  #my $t1 = [gettimeofday];
  #print "Time to run pfetch: ",tv_interval($t0,$t1),"\n";

  return ([$qual_5,$qual_3,$seqstr,$snp_qual]);
  
}

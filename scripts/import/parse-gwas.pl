#!/software/bin/env/usr/bin/env perl
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


use strict;
use Data::Dumper;

# define the names of the columns we want, as the order may change
our @COLS = (
  'PUBMEDID',
  'Disease/Trait',
  'Reported Gene(s)',
  'Strongest SNP-Risk Allele',
  'SNPs',
  'Risk Allele Frequency',
  'p-Value',
);

parse_file();
process_file();

sub parse_file {
  die "You must provide a filename to $0 to be parsed" unless @ARGV;
  
  open IN, $ARGV[0];
  
  # get header
  my $header = <IN>;
  chomp $header;
  my @headers = split /\t/, $header;
  
  my %hash;
  $hash{$headers[$_]} = $_ for (0..$#headers);
  
  open OUT, ">initial_data_new";
  
  # read file
  while(<IN>) {
	chomp;
	
	my @data = split /\t/, $_;
	
	# check number of cols
	next unless scalar @data == scalar @headers;
	
	print OUT (join "\t", map {$data[$hash{$_}]} @COLS);
	print OUT "\n";
  }
}

sub process_file {# join pubmed_id together and for each rsname in name list,put it in separate entry

  open IN, "initial_data_new";
  open OUT, ">initial_data2";

  while (<IN>) {
    next if /^PubMedID/;
    chomp;
    my @all=split /\t/, $_;
    $all[0] = "pubmed/$all[0]";
    my @rsnames = split /\,/, $all[4] if $all[4];
    foreach my $n (@rsnames) {
      print OUT join ("\t", $n,@all)."\n";
    }
  }
}

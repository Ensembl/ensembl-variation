#!/software/bin/env/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

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

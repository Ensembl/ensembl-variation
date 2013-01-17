#!/software/bin/env/usr/bin/env perl

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

#!/usr/bin/env perl

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
use Bio::EnsEMBL::Registry;
use Getopt::Long;
use FileHandle;

# COMMAND LINE OPTIONS
######################

my $config = {};

GetOptions(
	$config,
	'help',      # displays help message
	'input=s',   # input file
	'sample=s',  # sample name
);

# print usage message if requested or no args supplied
if(defined($config->{help})) {
	&usage;
	exit(0);
}

# defaults
$config->{sample} ||= 'DEFAULT_SAMPLE_NAME';

my $in_file_handle = new FileHandle;

if(defined($config->{input})) {
	
	# check defined input file exists
	die("ERROR: Could not find input file ", $config->{input}, "\n") unless -e $config->{input};
	
	# open file
	if($config->{input} =~ /\.gz$/){
		$in_file_handle->open("zcat ". $config->{input} . " | " ) or die("ERROR: Could not read from input file ", $config->{input}, "\n");
	}
	else {
		$in_file_handle->open( $config->{input} ) or die("ERROR: Could not read from input file ", $config->{input}, "\n");
	}
}
else {
	$in_file_handle = 'STDIN';
}

# header
print join "\t", (
	'#CHROM',
	'POS',
	'ID',
	'REF',
	'ALT',
	'QUAL',
	'FILTER',
	'INFO',
	'FORMAT',
	$config->{sample}."\n"
);

while(<$in_file_handle>) {
	chomp;
	
	my @data = split /\s+/, $_;
		
	# some files have an extra col at the beginning
	shift @data if $data[0] =~ /^\d+$/ && $data[1] =~ /chr/;
	
	# trim off characters from chrom name
	$data[0] =~ s/chr(om)?//ig;
	
	# alleles
	my @alleles = split /\//, $data[3];
	unshift @alleles, 'N' if scalar @alleles == 1;
	
	# get genotype
	my @gt = split /\//, $data[3];
	push @gt, $gt[0] if scalar @gt == 1;
	
	# make codes
	my %codes = map {$alleles[$_] => $_} 0..$#alleles;
	
	# make coded genotype
	my $code = join '/', map {$codes{$_}} @gt;
	
	# print data
	print join "\t", (
		$data[0],
		$data[1] + 1,
		'.',
		$alleles[0],
		$alleles[1],
		'.',
		'.',
		'.',
		'GT',
		$code."\n"
	);
}
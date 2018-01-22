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


use strict;

use Getopt::Long;
use FileHandle;

my $config = {};

GetOptions(
	$config,
	'input|i=s',
	'samples|s=s',
);

my $in_file_handle = new FileHandle;

if(defined($config->{input})) {
	
	# check defined input file exists
	die("ERROR: Could not find input file ", $config->{input}, "\n") unless -e $config->{input};
	
	if($config->{input} =~ /\.gz$/){
		$in_file_handle->open("zcat ". $config->{input} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
	}
	else {
		$in_file_handle->open( $config->{input} ) or die("ERROR: Could not read from input file ", $config->{input}, "\n");
	}
}

# no file specified - try to read data off command line
else {
	$in_file_handle = 'STDIN';
}

my %include_samples;

if(defined($config->{samples})) {
	open IN, $config->{samples} or die("ERROR: Could not open samples file ".$config->{samples}."\n");
	
	while(<IN>) {
		chomp;
		$include_samples{(split)[0]} = 1;
	}
	
	close IN;
}


my %headers;

while(<$in_file_handle>) {
	next if /^##/;
	
	my @split = split /\t/;
	my $data = {};
	$data->{line} = $_;
	
	# column definition line
	if(/^#/) {
		%headers = %{parse_header($config, \@split)};
	}
	
	# data
	else {
		$data->{$_} = $split[$headers{$_}] for keys %headers;
		
		# don't want any with more than one alt
		next if $data->{ALT} =~ /\,/;
		
		# don't want any without rs
		next unless $data->{ID} =~ /rs/;
		
		my $line = join "\t", (
			$data->{'#CHROM'},
			$data->{POS},
			$data->{POS},
			1,
		);
		
		for(my $i = $config->{first_sample_col}; $i <= $#split; $i++) {
			my $sample_id = $config->{individuals}->[$i - $config->{first_sample_col}];
			next if scalar keys %include_samples and !$include_samples{$sample_id};
			
			my $gt = (split /\:/, $split[$i])[0];
			
			# skip missing data
			next if $gt =~ /\./;
			
			$gt =~ s/\||\\|\///g;
			$gt =~ tr/01/Aa/;
			
			$sample_id =~ s/\D//g;
			$sample_id =~ s/^0+//g;
			
			print "$line\t$sample_id\t$gt\n";
		}
	}
}
	

sub parse_header {
	my $config     = shift;
	my $split_ref  = shift;
	
	my @split = @$split_ref;
	
	my %headers;
	$headers{$split[$_]} = $_ for(0..$#split);
	
	# set location of first sample col
	if(defined($headers{FORMAT})) {
		$config->{first_sample_col} = $headers{FORMAT} + 1;
		
		# splice @split hash to get just individual IDs
		splice(@split, 0, $config->{first_sample_col});
		
		$config->{individuals} = \@split;
	}
	
	return \%headers;
}
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
use warnings;

use Getopt::Long;
use FileHandle;

my ($freq_file, $ld_file, $out_file, $pop, $chr, $min_r2);

GetOptions(
	'freqs=s' => \$freq_file,
	'ld=s' => \$ld_file,
	'output=s' => \$out_file,
	'population=s' => \$pop,
	'chr=s' => \$chr,
	'r2=s' => \$min_r2,
);

# check args
die("ERROR: Frequencies file $freq_file not found\n") unless -e $freq_file;
#die("ERROR: No output file specified\n") unless defined($out_file);
die("ERROR: No population specified\n") unless defined($pop);
$min_r2 ||= 0.99; #default value for r2

# load MAFs
open IN, "<$freq_file" or die("ERROR: Could not open frequencies file $freq_file\n");


my $maf = {};

while(<IN>) {
    chomp;
    /^\CHROM/ && next;
	
    my ($chr, $pos, $n_a, $n_c, $a1, $a2) = split;
	
	$a1 =~ s/^.+://g;
	$a2 =~ s/^.+://g;
	
	# maf is lowest, not always first/last
	$maf->{$pos} = (sort {$a <=> $b} ($a1, $a2))[0];
}
close IN;

# load LD
my $ld_file_handle = new FileHandle;

if(defined($ld_file)) {
	$ld_file_handle->open("<$ld_file") or die("ERROR: Could not open LD file $ld_file\n");
}
else {
	$ld_file_handle = 'STDIN';
}


my $ld = {};

while(<$ld_file_handle>) {
	chomp;
	
	my @split = split;
	
	my ($p1, $p2, $r2) = @split[2..4];
	
	next unless $r2 > $min_r2;
	
	push @{$ld->{$p1}}, $p2;
}

$ld_file_handle->close;

#do the algorithm
my $remove_snps = {}; #hash containing the snps that must be removed from the entry, they have been ruled out
my $tagged_by = {};

foreach my $pos (sort {$maf->{$b} <=> $maf->{$a} || $a <=> $b} keys %{$maf}){
    if (!defined $remove_snps->{$pos}){
		#add the SNPs that should be removed in future iterations
		#and delete from the hash with the MAF_snps, the ones that have a r2 greater than r2 with $variation_feature_id
		#map {$remove_snps->{$_}++;delete $MAF_snps->{$_}} @{$LD_values->{$seq_region_start}};
		if(defined($ld->{$pos})) {
			foreach(@{$ld->{$pos}}) {
				$remove_snps->{$_}++;
				delete $maf->{$_};
				push @{$tagged_by->{$pos}}, $_;
			}
		}
		
		# also delete the ones that don't exist in $LD_values
		# since they can't be tag SNPs
		else {
			$remove_snps->{$pos}++;
			delete $maf->{$pos};
		}
    }
}

my $out_file_handle = new FileHandle;

if(defined($out_file)) {
	$out_file_handle->open(">$out_file") or die ("ERROR: Could not open output file $out_file\n");
}
else {
	$out_file_handle = *STDOUT;
}


foreach my $p1 (keys %{$maf}){
	
	foreach my $p2 (@{$tagged_by->{$p1}}) {
		print $out_file_handle join("\t", $chr, $p1, $p2, $pop);
		print $out_file_handle "\n";
	}
}
$out_file_handle->close or die ("Could not close output file with tagged SNPs");

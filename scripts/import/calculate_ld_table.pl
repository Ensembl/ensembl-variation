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

my ($TMP_DIR, $TMP_FILE, $ld_file);


GetOptions('tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	   'ldfile=s'  => \$ld_file);



$TMP_DIR = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

#added default options

#first of all, sort the file by position

`sort -n -k 1 -o $ld_file $ld_file` unless $ld_file =~ /gz$/; #order snps by position

#now, read the file, and convert it to the AA,Aa or aa format
my $seq_region_start;
my $individual_id;
my $population_id;
my $seq_region_id;
my $allele_1;
my $allele_2;
my $previous_seq_region_start = 0;

$ld_file =~ /dump_data_(\d+)\_(\d+)/;
$population_id = $1;
$seq_region_id = $2;
my %individual_information = (); #hash containing relation snps->individuals
my %alleles_information = (); #hash containing a record of alleles in the variation. A will be the major and a the minor. When more than
my $buffer = {};
#2 alleles per variation, genotype will be discarded
#get the seq_region_id and population_id from the file name
if($ld_file =~ /gz$/) {
	open IN, "gzip -dc $ld_file |";
}
else {
	open IN,"<$ld_file" or die "Could not open input file: $ld_file\n";
}

while (<IN>){
    chomp;
    ($seq_region_start,$individual_id, $allele_1,$allele_2) = split; #get all the fields in the file
	
	next unless $seq_region_start and $individual_id and $population_id and $allele_1 and $allele_2;
	
    if ($previous_seq_region_start == 0 or $seq_region_start == $previous_seq_region_start){
	$previous_seq_region_start = $seq_region_start;
	if ($allele_1 ne 'N' and $allele_2 ne 'N'){
	    $alleles_information{$population_id}{$allele_1}++;
	    $alleles_information{$population_id}{$allele_2}++;
	    
	    $individual_information{$population_id}{$individual_id}{allele_1} = $allele_1;
	    $individual_information{$population_id}{$individual_id}{allele_2} = $allele_2;    
	}
    }
    else{
	#at this point, we have seen all variations for individual, time to convert genotype and print to file
	foreach my $population (keys %alleles_information){
	    next if (keys %{$alleles_information{$population}} > 2); #skip variations with 3 alleles
	    convert_genotype($alleles_information{$population},$individual_information{$population});
	    #print all individuals to a file
	    map {print_buffered ($buffer,"$TMP_DIR/$TMP_FILE\.$population_id\.$seq_region_id",
				 join("\t",$seq_region_id,$previous_seq_region_start, $previous_seq_region_start,
				      $population, $_,
				      $individual_information{$population}{$_}{genotype})."\n")} keys %{$individual_information{$population}};
	}
	#and finally, empty structures
	%alleles_information = ();
	%individual_information = ();
	$previous_seq_region_start = $seq_region_start;
	if ($allele_1 ne 'N' and $allele_2 ne 'N'){
	    $alleles_information{$population_id}{$allele_1}++;
	    $alleles_information{$population_id}{$allele_2}++;
	    
	    $individual_information{$population_id}{$individual_id}{allele_1} = $allele_1;
	    $individual_information{$population_id}{$individual_id}{allele_2} = $allele_2;
	}
    }
}
close IN or die "Could not close input file";
#print remaining region, if has 2 alleles
foreach my $population (keys %alleles_information){
    if (keys %{$alleles_information{$population}} <= 2){
	convert_genotype($alleles_information{$population},$individual_information{$population});
	#print all individuals to a file
	map {print_buffered ($buffer,"$TMP_DIR/$TMP_FILE\.$population_id\.$seq_region_id",
			     join("\t",$seq_region_id,$previous_seq_region_start, $previous_seq_region_start,
				  $population, $_,
				  $individual_information{$population}{$_}{genotype})."\n")} keys %{$individual_information{$population}};
    }
}
print_buffered($buffer);
#and run ld calculation
my $file = "$TMP_DIR/$TMP_FILE\.$population_id\.$seq_region_id"; #file containing the genotype in the AA format
#once the files are created, we have to calculate the ld
my $call .= "calc_genotypes $file $file\.out";
#print $call,"\n";
system($call);
unlink("$TMP_DIR/$TMP_FILE\.$population_id\.$seq_region_id");

#
# Converts the genotype into the required format for the calculation of the pairwise_ld value: AA, Aa or aa
# From the Allele table, will select the alleles and compare to the alleles in the genotype
#

sub convert_genotype{
    my $alleles_information = shift; #reference to the hash containing the alleles for the variation present in the genotypes
    my $individual_information = shift; #reference to a hash containing the values to be written to the file
    my @alleles_ordered; #the array will contain the alleles ordered by apparitions in the genotypes (only 2 values possible)
    
    @alleles_ordered = sort({$alleles_information->{$b} <=> $alleles_information->{$a}} keys %{$alleles_information});
    
    #let's convert the allele_1 allele_2 to a genotype in the AA, Aa or aa format, where A corresponds to the major allele and a to the minor
    foreach my $individual_id (keys %{$individual_information}){
	    #if both alleles are different, this is the Aa genotype
	    if ($individual_information->{$individual_id}{allele_1} ne $individual_information->{$individual_id}{allele_2}){
		$individual_information->{$individual_id}{genotype} = 'Aa';
	    }
	    #when they are the same, must find out which is the major
	    else{	    
		if ($alleles_ordered[0] eq $individual_information->{$individual_id}{allele_1}){
		    #it is the major allele
		    $individual_information->{$individual_id}{genotype} = 'AA';
		}
		else{
		    $individual_information->{$individual_id}{genotype} = 'aa';
		}
		
	    }
	}
}

sub print_buffered {
    my $buffer = shift;
    my $filename = shift;
    my $text = shift;

    local *FH;

    if( ! $filename ) {
	# flush the buffer
	foreach my $file (keys %{$buffer}){
	    open( FH, ">>$file" ) or die;
	    print FH $buffer->{ $file };
	    close FH;
	}

    } else {
	$buffer->{ $filename } .= $text;
	if( length( $buffer->{ $filename } ) > 10_000 ) {
	    open( FH, ">>$filename" ) or die;
	    print FH $buffer->{ $filename };
	    close FH;
	    $buffer->{ $filename } = '';
	}
    }
}

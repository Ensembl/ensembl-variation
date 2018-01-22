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
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor;
use FindBin qw( $Bin );
use Data::Dumper;

my ($TMP_DIR, $population_id, $species, $registry_file);


GetOptions('species=s' => \$species,
	   'tmpdir=s'  => \$TMP_DIR,
	   'registry_file=s' => \$registry_file);

warn("Make sure you have a updated ensembl.registry file!\n");

$registry_file ||= $Bin . "/ensembl.registry";

#added default options
$species ||= 'human';
Bio::EnsEMBL::Registry->load_all( $registry_file );

my $dbVariation = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');

#first, calculate the MAF for the SNPs in the chromosome
my $file = "$ARGV[0]" if (defined @ARGV);
die "Not possible to calculate SNP tagging without file with SNPs" if (!defined @ARGV);

open IN, "<$file" or die("Could not open input file: $file\n");

my $r2 = 0.99; #default value for r2

my $previous_seq_region_start = 0;
my ($allele_1, $allele_2, $seq_region_start, $ind_id);
my $genotypes_snp = {};
my $last_snp = 0; #to indicate if there is a last snp that needs the MAF calculation
my $MAF_snps = {};
while (<IN>){
    chomp;
    /^\#/ && next;
    ($seq_region_start, $ind_id, $allele_1, $allele_2) = split;
	next unless $allele_1 and $allele_2 and $seq_region_start;
	
    #when we change variation_feature, calculate the MAF
    if (($previous_seq_region_start != $seq_region_start) && $previous_seq_region_start != 0){
		$MAF_snps->{$previous_seq_region_start} = &calculate_MAF($genotypes_snp);
		
		$genotypes_snp = {};
		$previous_seq_region_start = $seq_region_start;
		$last_snp = 1;
    }
    if ($previous_seq_region_start == 0){ #initialize the variables
	$previous_seq_region_start = $seq_region_start;
    }
    $genotypes_snp->{$allele_1}++ if ($allele_1 ne 'N');
    $genotypes_snp->{$allele_2}++ if ($allele_2 ne 'N');    
    $last_snp = 0;
}
close IN;

#calculate the MAF for the last SNP, if necessary
$MAF_snps->{$previous_seq_region_start} = &calculate_MAF($genotypes_snp)  if ($last_snp == 0);

#get LD values for the chromosome
$file =~/.*_(\d+)_(\d+)\.txt/; #extract the seq_region_id from the name of the file
$population_id = $1;
my $seq_region_id = $2;

my $pos_to_vf = {};
my $host = `hostname`;
chop $host;
my $LD_values = &get_LD_chromosome($dbVariation,$seq_region_id,$r2,$population_id, $pos_to_vf);
#do the algorithm
my $remove_snps = {}; #hash containing the snps that must be removed from the entry, they have been ruled out
my $tagged_by = {};
foreach $seq_region_start (sort {$MAF_snps->{$b} <=> $MAF_snps->{$a} || $a <=> $b} keys %{$MAF_snps}){
    if (!defined $remove_snps->{$seq_region_start}){
		#add the SNPs that should be removed in future iterations
		#and delete from the hash with the MAF_snps, the ones that have a r2 greater than r2 with $variation_feature_id
		#map {$remove_snps->{$_}++;delete $MAF_snps->{$_}} @{$LD_values->{$seq_region_start}};
		if(defined($LD_values->{$seq_region_start})) {
			foreach(@{$LD_values->{$seq_region_start}}) {
				$remove_snps->{$_}++;
				delete $MAF_snps->{$_};
				push @{$tagged_by->{$seq_region_start}}, $_;
			}
		}
		
		# also delete the ones that don't exist in $LD_values
		# since they can't be tag SNPs
		else {
			$remove_snps->{$seq_region_start}++;
			delete $MAF_snps->{$seq_region_start};
		}
    }
}
my $genotype_without_vf = 0;
my $genotype_with_vf = 0;
open OUT, ">$TMP_DIR/snps_tagged_$population_id\_$host\-$$\.txt" or die ("Could not open output file");
foreach my $position_vf (keys %{$MAF_snps}){
	# do pos_to_vf on the tagged snps
	foreach my $tagged(@{$tagged_by->{$position_vf}}) {
		$pos_to_vf->{$tagged} ||= &get_vf_id_from_position($dbVariation, $seq_region_id, $tagged);
	}
	
    if (! defined $pos_to_vf->{$position_vf}){ #some variations might not have LD, get dbID from database
		#get it from the database
		$pos_to_vf->{$position_vf} = &get_vf_id_from_position($dbVariation,$seq_region_id,$position_vf);
    }
    if ($pos_to_vf->{$position_vf} ne ''){
		print OUT join("\t",$pos_to_vf->{$position_vf},$pos_to_vf->{$_},$population_id),"\n" for grep {defined($pos_to_vf->{$_})} @{$tagged_by->{$position_vf}};
		$genotype_with_vf++;
    }
    else{
		$genotype_without_vf++;
    }
}
close OUT or die ("Could not close output file with tagged SNPs");
unlink($file);

#for a given position retrieve the vf_id from the database
sub get_vf_id_from_position{
	my $dbVar = shift;
    my $seq_region_id = shift;
    my $seq_region_start = shift;
	
	my $sth = $dbVar->dbc->prepare(qq{
	  SELECT variation_feature_id, source_id
	  FROM variation_feature
	  WHERE seq_region_id = ?
	  AND seq_region_start = ?
	  AND seq_region_end = seq_region_start
	});
	
	$sth->execute($seq_region_id, $seq_region_start);
	
	my ($vf_id, $source);
	$sth->bind_columns(\$vf_id, \$source);
	
	my %by_source;
	
	push @{$by_source{$source}}, $vf_id while $sth->fetch;
	$sth->finish;
	
	if(scalar keys %by_source) {
		foreach my $s(sort {$a <=> $b} keys %by_source) {
			return shift @{$by_source{$s}};
		}
	}
	
	return '';
}



#for a list of genotypes, get the MAF
sub calculate_MAF{
    my $genotypes_snp = shift;
    my $MAF;
    my $total = 0;
    my $allele_freq;

    if (keys %{$genotypes_snp} == 2){
		$total += $_ for values %{$genotypes_snp};
		
		foreach (values %{$genotypes_snp}){
			$allele_freq = $_ / $total;
			last;
		}
		return ($allele_freq < 0.5 ? $allele_freq : 1 - $allele_freq);
    }    
    return 0;
#    die "genotype with more than 2 alleles!!";
}

#creates a hash with all the variation_features in the chromosome with a r2 greater than r2
sub get_LD_chromosome{
    my $dbVariation = shift;
    my $seq_region_id = shift;
    my $r2 = shift;
    my $population_id = shift;
    
    my $variation_features = {};
    my $sth = $dbVariation->dbc->prepare(qq{SELECT seq_region_start,seq_region_end
						FROM pairwise_ld
						WHERE seq_region_id = ?
						AND r2 > ?
						AND sample_id = ?
					    },{mysql_use_result =>1});
    $sth->execute($seq_region_id,$r2,$population_id);
    my ($seq_region_start,$seq_region_end);
    $sth->bind_columns(\$seq_region_start, \$seq_region_end);
    while ($sth->fetch()){
	push @{$variation_features->{$seq_region_start}}, $seq_region_end;
	push @{$variation_features->{$seq_region_end}}, $seq_region_start;
    }
    $sth->finish();
    return $variation_features;
}

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
use FindBin qw( $Bin );
use Bio::EnsEMBL::Registry;
use Data::Dumper;

my ($species, $snpassay_file, $snpind_file,$population_name);



#
# bsub -q bigmem -W4:00 -R"select[mem>3500] rusage[mem=3500]" -M3500000
#

GetOptions('species=s' => \$species,
	   'snpassay_file=s' => \$snpassay_file,
	   'snpind_file=s' => \$snpind_file,
	   'population_name=s' => \$population_name
	   );
warn("Make sure you have a updated ensembl.registry file!\n");

my $registry_file ||= $Bin . "/ensembl.registry";
if ($population_name ne 'CELERA_STRAIN:SD' and $population_name ne 'ENSEMBL:STAR'){
    die "Only allow STAR and CELERA data for the moment\n\n";
}
#added default options
$species ||= 'rat';
Bio::EnsEMBL::Registry->load_all( $registry_file );

my $dbVariation = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');

open SNP, ">$snpassay_file" or die "could not open output file for snp assay info: $!\n";
#open IND, ">$snpind_file" or die "could not open output file for ind assay info: $!\n";

#first of all, write header for files
print_snp_headers();
print_ind_headers();
#get data
my $var_adaptor = $dbVariation->get_VariationAdaptor();
my $pop_adaptor = $dbVariation->get_PopulationAdaptor();
my $population = $pop_adaptor->fetch_by_name($population_name);
my $variations = $var_adaptor->fetch_all_by_Population($population);
my $vf_adaptor = $dbVariation->get_VariationFeatureAdaptor();
print_snp_data($variations,$vf_adaptor);


close SNP or die "Could not close snp assay info file: $!\n";
#close IND or die "Could not close ind assay info file: $!\n";


#will print the snpassay file header
sub print_snp_headers{
  
    print_cont_section(SNP); #contacts
    print_pub_section(SNP); #publications
    print_method_section(IND); #methods
    print_pop_section(SNP); #print population section

    #snpassay
    print SNP "TYPE:\tSNPASSAY\n";
    print SNP "HANDLE:\tENSEMBL\n";
    print SNP "BATCH:\t2007\n";  #use year of submission for batch ??
    print SNP "MOLTYPE:\tGenomic\n";
    print SNP "METHOD:\tEnsembl-SSAHA\n";
    print SNP "SAMPLESIZE:\t2\n"; #samplesize ??
    print SNP "ORGANISM:\tRattus norvegicus\n";
    print SNP "CITATION:\t\n";
    print SNP "POPULATION:\t",$population_name,"\n";
    print SNP "COMMENT:\t\n"; #any comment ??
    print SNP "||\n";

}


#for a list of variation objects, print the necessary information for dbSNP
sub print_snp_data{
    my $variations = shift;
    my $vf_adaptor = shift;


    foreach my $variation (@{$variations}){
	my $vf = shift @{$vf_adaptor->fetch_all_by_Variation($variation)};
	print SNP "SNP:\t",$variation->name,"\n";
	print SNP "ACCESSION:\t",$vf->slice->accession_number,"\n"; #is this the info they want ??
	print SNP "SAMPLESIZE:\t2\n"; #number of chromosomes ??
	print SNP "LENGTH:\t", length($variation->five_prime_flanking_seq) + length($variation->three_prime_flanking_seq) + 1,"\n";
	print SNP "5\'_FLANK:\t",$variation->five_prime_flanking_seq,"\n";
	print SNP "OBSERVED:\t", $vf->allele_string,"\n";
	print SNP "3\'_FLANK:\t",$variation->three_prime_flanking_seq,"\n";
	print SNP "||\n";
    }
}


#prints contacts section
sub print_cont_section{
    my $fh  = shift;

  #contact details
    print $fh "TYPE:\tCONT\n";
    print $fh "HANDLE:\tENSEMBL\n"; #should this be the source ??
    print $fh "NAME:\tDaniel Rios\n"; #my name ??
    print $fh "FAX:\t00441223494468\n";
    print $fh "TEL:\t00441223494684\n";
    print $fh "EMAIL:\tdani\@ebi.ac.uk\n";
    print $fh "LAB:\tEnsembl project\n";
    print $fh "INST:\tEuopearn Bioinformatics Institute\n";
    print $fh "ADDR:\tEMBL-EBI,Wellcome Trust Genome Campus,Hinxton,CB10 1SD Cambridge, UK\n";
    print $fh "||\n";

}

#prints pub section

sub print_pub_section{
    my $fh = shift;

   #publications
    print $fh "TYPE:\tPUB:\n";
    print $fh "||\n";

}

#print population section

sub print_pop_section{
    my $fh = shift;

    #population
    print $fh "TYPE:\tPOPULATION\n";
    print $fh "HANDLE:\tENSEMBL\n";
    print $fh "ID:\t",$population_name,"\n";
    print $fh "POP_CLASS:\tUNKNOWN\n";
    print $fh "POPULATION:\t"; #add population description ??
    print $fh "||\n";

}

sub print_method_section{
    my $fh = shift;

    #method
    print $fh "TYPE:\tMETHOD\n";
    print $fh "HANDLE:\tENSEMBL\n";
    print $fh "ID:\tEnsembl-SSAHA\n"; #which method ??
    print $fh "METHOD_CLASS:\tComputation\n";
    print $fh "SEQ_BOTH_STRANDS:\tNA\n"; #both strands ??
    print $fh "TEMPLATE_TYPE:\tUNKNOWN\n";
    print $fh "MULT_PCR_AMPLIFICATION:\tNA\n";
    print $fh "MULT_CLONES_TESTED:\tNA\n";
    print $fh "METHOD:\tComputationally discovered SNPs usng SSAHA\n"; #another comment ??
    print $fh "||\n";

}

#prints the headers for the snpinduse file
sub print_ind_headers{

    print_cont_section(IND); #contacts section
    print_pub_section(IND); #publications section
    print_pop_section(IND); #population section
    #print some individual specific section
    print_method_section(IND); #method section
    
}

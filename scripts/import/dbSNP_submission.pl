#!/usr/local/ensembl/bin/perl
use strict;
use warnings;

use Getopt::Long;
use FindBin qw( $Bin );
use Bio::EnsEMBL::Registry;
use Data::Dumper;

my ($species, $snpassay_file, $snpind_file,$population_name);


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
    #contact details
    print SNP "TYPE:\tCONT\n";
    print SNP "HANDLE:\tENSEMBL\n"; #should this be the source ??
    print SNP "NAME:\tDaniel Rios\n"; #my name ??
    print SNP "FAX:\t00441223494468\n";
    print SNP "TEL:\t00441223494684\n";
    print SNP "EMAIL:\tdani\@ebi.ac.uk\n";
    print SNP "LAB:\tEnsembl project\n";
    print SNP "INST:\tEuopearn Bioinformatics Institute\n";
    print SNP "ADDR:\tEMBL-EBI,Wellcome Trust Genome Campus,Hinxton,CB10 1SD Cambridge, UK\n";
    print SNP "||\n";
    #publications
    print SNP "TYPE:\tPUB:\n";
    print SNP "||\n";
    #method
    print SNP "TYPE:\tMETHOD\n";
    print SNP "HANDLE:\tENSEMBL\n";
    print SNP "ID:\tEnsembl-SSAHA\n"; #which method ??
    print SNP "METHOD_CLASS:\tComputation\n";
    print SNP "SEQ_BOTH_STRANDS:\tNA\n"; #both strands ??
    print SNP "TEMPLATE_TYPE:\tUNKNOWN\n";
    print SNP "MULT_PCR_AMPLIFICATION:\tNA\n";
    print SNP "MULT_CLONES_TESTED:\tNA\n";
    print SNP "METHOD:\tComputationally discovered SNPs usng SSAHA\n"; #another comment ??
    print SNP "||\n";
    #population
    print SNP "TYPE:\tPOPULATION\n";
    print SNP "HANDLE:\tENSEMBL\n";
    print SNP "ID:\t",$population_name,"\n";
    print SNP "POP_CLASS:\tUNKNOWN\n";
    print SNP "POPULATION:\t"; #add population description ??
    print SNP "||\n";
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

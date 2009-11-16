#!/usr/local/ensembl/bin/perl
use lib '/nfs/users/nfs_y/yuan/ensembl/src/ensembl-variation-55/modules';

use strict;
use warnings;

use Getopt::Long;
use FindBin qw( $Bin );
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use ImportUtils qw(dumpSQL load debug);

my ($species, $snpassay_file, $snpind_file,$population_name,$seq_region_id,$TMP_DIR,$TMP_FILE);



#
# bsub -q normal -W4:00 -R"select[mem>3500] rusage[mem=3500]" -M3500000
#

GetOptions('species=s' => \$species,
	   'tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	   'seq_region_id=i' => \$seq_region_id,
           'population_name=s' => \$population_name
	   );
warn("Make sure you have a updated ensembl.registry file!\n");

my $registry_file ||= $Bin . "/ensembl.registry";

#need change for different submisson
my $sample_size = 8;
#my $sample_size = 932;
my $pop_class ='Unknown';
my $tax_id = 10116;
my $batch = '2009-11_STAR_4_strain';
#my $batch = '2009-11_STAR-genotype';
my $orgainism_name = "Rattus Norvegicus";
my $seq_source = "STAR";
my $source_type = 'submitter'; # choice of repository or curator or submitter 
my $pop_source = "NA";
my $place_source = "NA";
my $breed = "NA";
my $sex = "Unknown";
my $method = 'RAT_STRAIN-READS_SNPS_200712'; #'RAT_STRAIN-GENOTYPES_200712''Ensembl-SSAHA';

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $dbVariation = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');

$TMP_DIR = $ImportUtils::TMP_DIR;
$ImportUtils::TMP_FILE .= $seq_region_id ? "\_$seq_region_id" : '';
$TMP_FILE = $ImportUtils::TMP_FILE;
$seq_region_id ||='';

my $new_pop_name = $population_name;
$new_pop_name =~ s/\:/\_/;#dbSNP more like '_' then ':'

print "new_pop_name is $new_pop_name\n";

open SNP, ">$TMP_DIR/snpassay_file\_$seq_region_id" or die "could not open output file for snp assay info: $!\n";
open IND, ">$TMP_DIR/snpind_file\_$seq_region_id" or die "could not open output file for ind assay info: $!\n";

#get data
my $var_adaptor = $dbVariation->get_VariationAdaptor();
my $pop_adaptor = $dbVariation->get_PopulationAdaptor();
#my $population = $pop_adaptor->fetch_by_name($population_name);
my $ind_adaptor = $dbVariation->get_IndividualAdaptor();
#my $variations = $var_adaptor->fetch_all_by_Population($population);
my $vf_adaptor = $dbVariation->get_VariationFeatureAdaptor();
my $ind_gtype_adaptor = $dbVariation->get_IndividualGenotypeAdaptor();
my $slice_adaptor = $dbCore->get_SliceAdaptor();

#first of all, write header for files
print_snp_headers();
print_ind_headers();
#print_snp_data($variations,$vf_adaptor,$ind_gtype_adaptor);#this is for small amonut of SNPs
print_snp_data_whole($dbVariation);#this is for big amount of SNPs

close SNP or die "Could not close snp assay info file: $!\n";
close IND or die "Could not close ind assay info file: $!\n";


#will print the snpassay file header
sub print_snp_headers{

    #print_cont_section(); #contacts
    #print_pub_section(); #publications
    #print_method_section(); #methods
    print_pop_section(); #print population section
    print_individual_section($pop_adaptor,$ind_adaptor); #print individual section


    #snpassay
    print SNP "TYPE:SNPASSAY\n";
    print SNP "HANDLE:ENSEMBL\n";
    print SNP "BATCH:$batch\n";  #use year of submission for batch ?? #need change for different submisson
    print SNP "MOLTYPE:Genomic\n";
    #print SNP "METHOD:Ensembl-SSAHA\n";
    print SNP "METHOD:$method\n";
    print SNP "SAMPLESIZE:$sample_size\n"; #samplesize ??
    print SNP "ORGANISM:$orgainism_name\n";#need change for different submisson
    print SNP "CITATION:\n";
    print SNP "POPULATION:",$new_pop_name,"\n";
    #need change for different submisson
    #print SNP "COMMENT:Mouse strain : A/J,129X1/SvJ,C3HeB/FeJ,129S1/SvImJ,DBA/2J,NOD/DIL,MSM/Ms\n"; #any comment ??
    #print SNP "COMMENT: Genomics sequences are from three individuals : gsc-plamid,wibr-plasmid and gsc-BAC\n";
    print SNP "||\n";

}

sub print_snp_data_whole {
  my $dbVar = shift;
  my $LIMIT = $seq_region_id ? "AND vf.seq_region_id = $seq_region_id" : '';

  my $sql=qq(SELECT vf.variation_name,vf.allele_string,tg.allele_1,tg.allele_2,s.name,f.seq_region_id,f.up_seq_region_start,f.up_seq_region_end,f.down_seq_region_start,f.down_seq_region_end
              FROM variation_feature vf,sample s, flanking_sequence f, tmp_individual_genotype_single_bp tg
              WHERE vf.variation_id=f.variation_id
              AND vf.variation_id=tg.variation_id
              AND tg.sample_id=s.sample_id 
              $LIMIT
              );

  dumpSQL($dbVar->dbc,$sql);

  system("sort $TMP_DIR/$TMP_FILE >$TMP_DIR/$TMP_FILE\.s");
  system("mv $TMP_DIR/$TMP_FILE\.s $TMP_DIR/$TMP_FILE");

  open FH, "$TMP_DIR/$TMP_FILE" or die "can't open $TMP_FILE file";

  my ($var,$pre_var_name,$up_seq,$down_seq);

  while (<FH>) {

    my ($variation_name,$allele_string,$allele_1,$allele_2,$ind_name,$seq_region_id,$up_seq_region_start,$up_seq_region_end,$down_seq_region_start,$down_seq_region_end) = split;
    next if ($variation_name !~ /^ENS/);
    #print "$variation_name,$allele_string,$allele_1,$allele_2,$ind_name,$seq_region_id,$up_seq_region_start,$up_seq_region_end,$down_seq_region_start,$down_seq_region_end\n";
    if ($pre_var_name and $variation_name ne $pre_var_name) {
      print SNP "SNP:",$var->{'var_name'},"\n";
      print IND "SNP:ENSEMBL|",$var->{'var_name'},"|SS_STRAND_FWD\n";
      print SNP "SAMPLESIZE:$sample_size\n"; #number of chromosomes ??
      print SNP "LENGTH:", length($var->{'up_seq'}) + length($var->{'down_seq'}) + 1,"\n";
      print SNP "5\'_FLANK:",$var->{'up_seq'},"\n";
      print SNP "OBSERVED:", $var->{'allele_string'},"\n";
      print SNP "3\'_FLANK:",$var->{'down_seq'},"\n";
      print SNP "||\n";
      foreach my $ind_na (@{$var->{'ind_name'}}) {
	next if !$var->{$ind_na}->{'allele_1'};
	my $allele1 = $var->{$ind_na}{'allele_1'};
	my $allele2 = $var->{$ind_na}{'allele_2'};
	print IND "ID:ENSEMBL|$new_pop_name:$ind_na:$allele1/$allele2\n";
      }
      print IND "||\n";
      undef $var;
    }

    ($up_seq,$down_seq) = get_flanking_seq($seq_region_id,$up_seq_region_start-300,$up_seq_region_end,$down_seq_region_start,$down_seq_region_end+300) if !$var->{'up_seq'};

    $pre_var_name = $variation_name;
    $var->{'var_name'} = $variation_name;
    $var->{'allele_string'} = $allele_string;
    $var->{$ind_name}{'allele_1'} = $allele_1;
    $var->{$ind_name}{'allele_2'} = $allele_2;
    push @{$var->{'ind_name'}},$ind_name;
    $var->{'up_seq'} = $$up_seq if $up_seq;
    $var->{'down_seq'} = $$down_seq if $down_seq;
  }

  if ($var) {
    #print last variation
    print SNP "SNP:",$var->{'var_name'},"\n";
    print IND "SNP:ENSEMBL|",$var->{'var_name'},"|SS_STRAND_FWD\n";
    print SNP "SAMPLESIZE:$sample_size\n"; #number of chromosomes ??
    print SNP "LENGTH:", length($var->{'up_seq'}) + length($var->{'down_seq'}) + 1,"\n";
    print SNP "5\'_FLANK:",$var->{'up_seq'},"\n";
    print SNP "OBSERVED:", $var->{'allele_string'},"\n";
    print SNP "3\'_FLANK:",$var->{'down_seq'},"\n";
    print SNP "||\n";
    foreach my $ind_name (@{$var->{'ind_name'}}) {
      next if !$var->{$ind_name}->{'allele_1'};
      my $allele_1 = $var->{$ind_name}->{'allele_1'};
      my $allele_2 = $var->{$ind_name}->{'allele_2'};
      print IND "ID:ENSEMBL|$new_pop_name:$ind_name:$allele_1/$allele_2\n" if ($allele_1 and $allele_2);
    }
    print IND "||\n";
  }
}

sub get_flanking_seq {
 my ($seq_region_id,$up_seq_region_start,$up_seq_region_end,$down_seq_region_start,$down_seq_region_end) = @_;

 my $tmp_slice = $slice_adaptor->fetch_by_seq_region_id($seq_region_id);

 if ($up_seq_region_start<1) {
   $up_seq_region_start=1;
 }
 elsif ($down_seq_region_end>$tmp_slice->length) {
   $down_seq_region_end = $tmp_slice->length;
 }

 my $up_slice = $tmp_slice->sub_Slice($up_seq_region_start,$up_seq_region_end,1);
 my $down_slice = $tmp_slice->sub_Slice($down_seq_region_start,$down_seq_region_end,1);

 if (!$up_slice or !$down_slice) {
   debug("Can't make up_slice or down slice with this line : $seq_region_id,$up_seq_region_start,$up_seq_region_end,$down_seq_region_start,$down_seq_region_end");
 }

 my $up_seq = $up_slice->seq;
 my $down_seq = $down_slice->seq;
 return (\$up_seq,\$down_seq);
}

#for a list of variation objects, print the necessary information for dbSNP
sub print_snp_data{
    my $variations = shift;
    my $vf_adaptor = shift;
    my $ind_gtype_adaptor = shift;

    foreach my $variation (@{$variations}){
	my $vf = shift @{$vf_adaptor->fetch_all_by_Variation($variation)};
	
	print SNP "SNP:",$variation->name,"\n";
	print IND "SNP:ENSEMBL|",$variation->name,"|SS_STRAND_FWD\n";
	print SNP "ACCESSION:",$vf->slice->accession_number,"\n"; #is this the info they want ??
	print SNP "SAMPLESIZE:$sample_size\n"; #number of chromosomes ??
	print SNP "LENGTH:", length($variation->five_prime_flanking_seq) + length($variation->three_prime_flanking_seq) + 1,"\n";
	print SNP "5\'_FLANK:",$variation->five_prime_flanking_seq,"\n";
	print SNP "OBSERVED:", $vf->allele_string,"\n";
	print SNP "3\'_FLANK:",$variation->three_prime_flanking_seq,"\n";
	print SNP "||\n";

	foreach my $ig (@{$ind_gtype_adaptor->fetch_all_by_Variation($variation)}){
	  my $ind_name = $ig->individual()->name;
	  print IND "ID:ENSEMBL|$new_pop_name:$ind_name",$ig->allele1,"/",$ig->allele2,"\n";
	}
	print IND "||\n";
    }
}


#prints contacts section
sub print_cont_section{

  #contact details
    print SNP "TYPE:CONT\n";
    print SNP "HANDLE:ENSEMBL\n"; #should this be the source ??
    print SNP "NAME:Yuan Chen\n"; #my name ??
    print SNP "FAX:00441223494468\n";
    print SNP "TEL:00441223494684\n";
    print SNP "EMAIL:yuan\@ebi.ac.uk\n";
    print SNP "LAB:Ensembl project\n";
    print SNP "INST:European Bioinformatics Institute\n";
    print SNP "ADDR:EMBL-EBI,Wellcome Trust Genome Campus,Hinxton,CB10 1SD Cambridge, UK\n";
    print SNP "||\n";

}

#prints pub section

sub print_pub_section{

   #publications
    print SNP "TYPE:PUB:\n";
    print SNP "||\n";

}

#print population section

sub print_pop_section{

    #population
    print SNP "TYPE:POPULATION\n";
    print SNP "HANDLE:ENSEMBL\n";
    print SNP "ID:",$new_pop_name,"\n";
    print SNP "POP_CLASS:$pop_class\n";
    print SNP "POPULATION:\n"; #add population description ??
    print SNP "||\n";

}

sub print_individual_section{

  my $pop_adaptor = shift;
  my $ind_adaptor = shift;

  my $pop = $pop_adaptor->fetch_by_name($population_name);
  my @ind_names = map {$_->name} @{$ind_adaptor->fetch_all_by_Population($pop)};

    print "individual_names are @ind_names\n";
    #individual
    foreach my $ind_name (@ind_names) {
      if ($ind_name =~/gsc/i) {
	$seq_source = "Genoscope";
	$pop_source = "NA";
	$place_source = "NA";
      }
      elsif ($ind_name =~ /wibr/i) {
	$seq_source = "Broad Institute (MIT)";
	$pop_source = "NA";
	$place_source = "NA";
      }
      elsif ($ind_name =~ /abeii/) {
	$seq_source = "abelii";
	$pop_source = "NA";
	$place_source = "NA";
      }
      print SNP "TYPE:INDIVIDUAL\n";
      #print SNP "HANDLE:ENSEMBL\|$new_pop_name\|$ind_name\|$tax_id\|$sex\|$breed\|$pop_source\n";
      print SNP "IND:ENSEMBL\|$new_pop_name\|$ind_name\|$tax_id\|$sex\|$breed\|$pop_source\n"; #dbSNP suggest use IND rather than HANDLE in the begining
      print SNP "SOURCE:$source_type\|$seq_source|$ind_name\|$place_source\n";
      #print SNP "PEDIGREE:NA\|NA\|NA\|NA\|NA\n";delete suggested by dbSNP
      print SNP "||\n";
    }
}

sub print_method_section{

    #method
    print SNP "TYPE:METHOD\n";
    print SNP "HANDLE:ENSEMBL\n";
    print SNP "ID:Ensembl-SSAHA\n"; #which method ??
    print SNP "METHOD_CLASS:Computation\n";
    print SNP "SEQ_BOTH_STRANDS:NA\n"; #both strands ??
    print SNP "TEMPLATE_TYPE:UNKNOWN\n";
    print SNP "MULT_PCR_AMPLIFICATION:NA\n";
    print SNP "MULT_CLONES_TESTED:NA\n";
    #need change for different submisson
    #print SNP "METHOD:Computationally discovered SNPs using ssahaSNP (see ssahaSNP description: http://www.sanger.ac.uk/Software/analysis/ssahaSNP/). The sequencing reads were from three individuals, two from gsc plasmid and BAC clones and one from wibr plasmid clone. These reads are component of reference sequence. The default parameters for ssahaSNP were used.\n";
    #print SNP "METHOD:Computationally discovered SNPs using ssahaSNP (see ssahaSNP description: http://www.sanger.ac.uk/Software/analysis/ssahaSNP/). The sequencing reads were from one individual abelii. These reads are component of reference sequence. The default parameters for ssahaSNP were used.\n"; #another comment ??
    print SNP "METHOD:Using ssaha2SNP, mouse sequencing reads from several sources were aligned to the top_level coordinates associated with build NCBI37 of the mouse genome. The read set includes whole genome sequencing reads from the DBA/2J, 129S1/ImJ, 129X1/SvJ, A/J strains of Mus musculus, C3HeB/FeJ strain of Mus musculus, BAC-end sequence reads from the NOD strain of Mus musculus, and additional sequencing reads from the MSM/Ms strain of Mus musculus molossinus. The reads for C3HeB/FeJ and NOD were generated at the Wellcome Trust Sanger Institute, the MSM-Ms reads were generated by Riken. Reads for the other strains were generated by Celera\n";
    print SNP "||\n";


}

#prints the headers for the snpinduse file
sub print_ind_headers{

  #contact details
    print IND "TYPE:CONT\n";
    print IND "HANDLE:ENSEMBL\n"; #should this be the source ??
    print IND "NAME:Yuan Chen\n"; #my name ??
    print IND "FAX:00441223494468\n";
    print IND "TEL:00441223494684\n";
    print IND "EMAIL:yuan\@ebi.ac.uk\n";
    print IND "LAB:Ensembl project\n";
    print IND "INST:Euopearn Bioinformatics Institute\n";
    print IND "ADDR:EMBL-EBI,Wellcome Trust Genome Campus,Hinxton,CB10 1SD Cambridge, UK\n";
    print IND "||\n";
 
    print IND "TYPE:SNPINDUSE\n";
    print IND "HANDLE:ENSEMBL\n";
    print IND "BATCH:$batch\n";
    #print IND "METHOD:Ensembl-SSAHA\n";
    print IND "METHOD:$method\n";
    print IND "||\n";

    #print_cont_section(IND); #contacts section
    #print_pub_section(IND); #publications section
    #print_pop_section(IND); #population section
    #print some individual specific section
    #print_method_section(IND); #method section

}

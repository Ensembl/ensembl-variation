#!/usr/bin/env perl
# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


use strict;
use warnings;

use Getopt::Long;
use FindBin qw( $Bin );
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use ImportUtils qw(dumpSQL load debug);

my ($species, $snpassay_file, $snpind_file,$population_name,$seq_region_id,$config_file,$TMP_DIR,$TMP_FILE);

=head
The example of config file :
#need change for different submisson
#The method section is already generated in a separated file, otherwise should be put inside config file 
#If unsure about a particular option, look at dbSNP web page : http://www.ncbi.nlm.nih.gov/SNP/how_to_submit.html
$sample_size = 4
$pop_class = Unknown
$tax_id = 9258
$batch = 2009-12_platypus_ssaha
$organism_name = Ornithorhynchus anatinus
$seq_source = WUGSC/NISC
$source_type = repository
$pop_source = NA
$place_source = NA
$breed = NA
$sex = Unknown
$method = PLATYPUS-READS_SNPS_200801

=cut

#
# bsub -q normal  -R"select[mem>3500] rusage[mem=3500]" -M3500000
#

GetOptions('species=s' => \$species,
	   'tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	   'seq_region_id=i' => \$seq_region_id,
           'population_name=s' => \$population_name,
	   'config_file=s' => \$config_file,
	   );
warn("Make sure you have a updated ensembl.registry file!\n");

my $registry_file ||= $Bin . "/ensembl.registry";

my %config;
if(defined $config_file) {
  open IN, $config_file or die "Could not read from config file $config_file\n";

  while(<IN>) {
    chomp;
    my ($option,$value) = split /\=/, $_;
    $option =~ s/\s+|^\$//g; #incase there is a $ sign or spaces
    $value =~ s/\s+// if $value =~ /\s+/;
    $config{$option} = $value;
  }
}
#need change for different submisson, see config file
my $sample_size = $config{'sample_size'};
my $pop_class =$config{'pop_class'};
my $tax_id = $config{'tax_id'};
my $batch = $config{'batch'};
my $organism_name = $config{'organism_name'};
my $seq_source = $config{'seq_source'};
my $source_type = $config{'source_type'}; # choice of repository or curator or submitter 
my $pop_source = $config{'pop_source'};
my $place_source = $config{'place_source'};
my $breed = $config{'breed'};
my $sex = $config{'sex'};
my $method = $config{'method'}; #'RAT_STRAIN-GENOTYPES_200712''Ensembl-SSAHA';
my $method_description = $config{'method_description'};

print "Check is this correct ?? sample_size is $sample_size, pop_class is $pop_class, tax_id is $tax_id, batch is $batch, organism_name is $organism_name, seq_source is $seq_source, source_type is $source_type, pop_source is pop_source, place_source is $place_source, breed is $breed, sex is $sex, method is $method\n";
sleep(60);

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
my $ind_adaptor = $dbVariation->get_IndividualAdaptor();
my $vf_adaptor = $dbVariation->get_VariationFeatureAdaptor();
my $ind_gtype_adaptor = $dbVariation->get_IndividualGenotypeAdaptor();
my $slice_adaptor = $dbCore->get_SliceAdaptor();

#first of all, write header for files
print_snp_headers();
print_ind_headers();
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
    print SNP "METHOD:$method\n";
    print SNP "SAMPLESIZE:$sample_size\n"; #samplesize ??
    print SNP "ORGANISM:$organism_name\n";#need change for different submisson
    print SNP "CITATION:\n";
    print SNP "POPULATION:",$new_pop_name,"\n";
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
              AND vf.source_id=2
              $LIMIT
              );

  dumpSQL($dbVar->dbc,$sql);

  system("sort $TMP_DIR/$TMP_FILE >$TMP_DIR/$TMP_FILE\.s");
  system("mv $TMP_DIR/$TMP_FILE\.s $TMP_DIR/$TMP_FILE");

  open FH, "$TMP_DIR/$TMP_FILE" or die "can't open $TMP_FILE file";

  my ($var,$pre_var_name,$up_seq,$down_seq);

  while (<FH>) {

    my ($variation_name,$allele_string,$allele_1,$allele_2,$ind_name,$seq_region_id,$up_seq_region_start,$up_seq_region_end,$down_seq_region_start,$down_seq_region_end) = split /\t/;
    next if ($variation_name !~ /^ENS|PCAP/);
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

    my $num = 400 - ($up_seq_region_end - $up_seq_region_start +1);

    ($up_seq,$down_seq) = get_flanking_seq($seq_region_id,$up_seq_region_start-$num,$up_seq_region_end,$down_seq_region_start,$down_seq_region_end+$num) if !$var->{'up_seq'};

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
 if ($down_seq_region_end>$tmp_slice->length) {
   $down_seq_region_end = $tmp_slice->length;
 }

 my $up_slice = $tmp_slice->sub_Slice($up_seq_region_start,$up_seq_region_end,1);
 my $down_slice = $tmp_slice->sub_Slice($down_seq_region_start,$down_seq_region_end,1);

 if (!$up_slice or !$down_slice) {
   debug("Can't make up_slice or down slice with this line : $seq_region_id,$up_seq_region_start,$up_seq_region_end,$down_seq_region_start,$down_seq_region_end and the length of slice is ",$tmp_slice->length,"\n");
 }

 my $up_seq = $up_slice->seq;
 my $down_seq = $down_slice->seq;
 return (\$up_seq,\$down_seq);
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
      print SNP "TYPE:INDIVIDUAL\n";
      print SNP "IND:ENSEMBL\|$new_pop_name\|$ind_name\|$tax_id\|$sex\|$breed\|$pop_source\n"; #dbSNP suggest use IND rather than HANDLE in the begining
      print SNP "SOURCE:$source_type\|$seq_source|$ind_name\|$place_source\n";
      print SNP "||\n";
    }
}

sub print_method_section{

    #method
    print SNP "TYPE:METHOD\n";
    print SNP "HANDLE:ENSEMBL\n";
    print SNP "$method\n"; ID:Ensembl-SSAHA\n"; #which method ??
    print SNP "METHOD_CLASS:Computation\n";
    print SNP "SEQ_BOTH_STRANDS:NA\n"; #both strands ??
    print SNP "TEMPLATE_TYPE:UNKNOWN\n";
    print SNP "MULT_PCR_AMPLIFICATION:NA\n";
    print SNP "MULT_CLONES_TESTED:NA\n";
    #need change for different submisson
    print SNP "$method_description\n";
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

}

#!/usr/bin/env perl

use strict;
use warnings;
use Bio::EnsEMBL::Variation::Utils::QCUtils qw( get_reference_base );


#input from nextstrain_json2geno.pl
#nuc     20692   C20692T C       T       ORF8    V114I   Shanghai/SH0114/2020    unassigned      NODE_0000451

my %samples;
my %info;
my %alt;

die "\nUsage: nextstrain_geno2vcf.pl [variants file] [genome.fa] \n\n" unless scalar(@ARGV)==2;

## use sequence accession
my $seq_name = "MN908947.3";
## get ref from genome sequence
my $db = Bio::DB::Fasta->new( $ARGV[1] )||die "Failed to index fasta for checking: $!\n";


open my $in, $ARGV[0] ||die "Failed to open $ARGV[0] : $!\n";
open my $out, ">", "$ARGV[0]\.vcf" ||die "Failed to open file to write : $!\n";

while(<$in>){

  next if /\#/;
  chomp;
  my ($s, $pos, $name, $ref, $alt, $gene, $prot, $samp, $clade, $parent) = split/\t/;


  ### group all data by reference position only

  ## store alleles for sample name
  $samples{$samp}{$pos}{alt} = $alt;
  $samples{$samp}{$pos}{ref} = $ref;

  $alt{$pos}{$alt} = 1;
  ## the reference allele may be the alt wrt the reference assembly
  ## saving this for checking purposes
  $alt{$pos}{$ref} = 1;


  $info{$pos}{gene} = $gene if $gene ne ".";

  ## store DNA change descriptions
  $info{$pos}{nuc}{$name} = 1 ;

  ## protein changes by alt
  $info{$pos}{prot}{$alt} = $prot;

  ## could be on different banches
  $info{$pos}{clade}{$clade} = 1 if defined $clade; 
}

my $header = format_header(\%samples);
print $out $header;


## write line per ref position
foreach my $var( sort sort_num(keys %info)){

  print $out "$seq_name\t$var\t.\t";

  ## get expected ref from genome assembly 
  my $ref_seq = $db->seq($seq_name, $var, $var);


  ### deal with alleles - could be multi-allelic
  my $n =1; 
  my %al2code;
  my %code2al;

  ## fix genomic ref as 0
  $al2code{$ref_seq} = 0;
  $code2al{0} = $ref_seq;

  foreach my $alt (keys %{$alt{$var}} ){
    next if defined $al2code{$alt} && $al2code{$alt} == 0; ## skip ref
    $al2code{$alt} = $n ;
    $code2al{$n} = $alt;
    $n++;
  }

  ## collect alts
  my @ordered_alts;
  foreach ( my $n=1; $n < (keys %code2al); $n++){
    push @ordered_alts, $code2al{$n};
  }

  print $out "$ref_seq\t". join(",",@ordered_alts) . "\t.\tPASS\t";



  ## add extracted info  
  my @info;
  push @info, "Gene=$info{$var}{gene}" if defined $info{$var}{gene};


  #order protein changes by alt allele
  my @prot;
  foreach my $a(@ordered_alts){
    push @prot, $info{$var}{prot}{$a};
  }
  push @info, "Prot=".join(",",@prot) if scalar(@prot >0);;
   

  ## clades are NOT ordered
  my $clades = join(",", (keys %{$info{$var}{clade}})) if defined $info{$var}{clade};
  push @info, "clades=$clades" if defined $clades;


  ## names are NOT ordered
  my $names = join(",", (keys %{$info{$var}{nuc}})) ;
  push @info, "names=$names";
  

  ## count alleles and store genotypes
  my($ac, $an, $genos) =  format_genotypes(\%samples, \%al2code, $var);

  ## add allele counts to info
  push @info, "AN=$an";
  push @info, "AC=$ac";

  print $out join(";", @info) ;

  print $out "\tGT\t$genos\n";

}

## generic lines & 
## order samples
sub format_header{

  my $samples = shift;

  my $h =  "##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description=\"All variants set as passed\">
##fileDate=20200408
##source http://data.nextstrain.org/ncov.json
##contig=<ID=MN908947.3,assembly=GCA009858895.3,length=29903>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
#CHROM\tPOS\tID\tREF\tALT\tFILTER\tQUAL\tINFO\tFORMAT\t";

  foreach my $sam (sort (keys %{$samples})){
    $h .= $sam."\t";
  }
  $h .= "\n";

  return $h;
}


## encode genotypes & calculate alleles  
sub format_genotypes{

  my ($samples, $al2code, $var) = @_;

  my %counts;
  my $an;
  my $genos;
  my $swap = 0;

  foreach my $sam (sort(keys %$samples)){

    ## count total alleles
    $an++;

    if (defined $samples->{$sam}{$var}{ref} ){

      ### check ref swap
      if($al2code->{ $samples->{$sam}{$var}{ref} } ne "0"){
        $swap++;
      #  print "Ref swap:\t$var\t$sam\t$samples->{$sam}{$var}{ref}\t$samples->{$sam}{$var}{alt}\n";
      }
      ## save geno to write
      unless( defined $al2code->{ $samples->{$sam}{$var}{alt}} ){
        die "no code for $sam $var ($samples->{$sam}{$var}{alt}})\n";
      }
      $genos .= "$al2code->{ $samples->{$sam}{$var}{alt} }\t";

      ## allow for multi-allelic; skipping reference
      $counts{ $al2code->{ $samples->{$sam}{$var}{alt} } }++
        unless $al2code->{ $samples->{$sam}{$var}{alt} } eq "0";

    }
    else{
      ## assume reference unless stated otherwise
      $genos .= "0\t";
    }
  }
  print "$var $swap reference swaps\n" if $swap > 0;

  ## order alt allele counts;
  my @ac;
   for (my $n=1;$n<=scalar(keys %counts) ; $n++){
    push @ac, $counts{$n};
  }
  my $ac =  join(",", @ac);

  return ($ac, $an, $genos);
}

sub sort_num{
  return $a<=>$b;
}


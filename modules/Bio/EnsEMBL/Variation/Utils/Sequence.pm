# EnsEMBL module for Bio::EnsEMBL::Variation::Utils::Sequence
#
#

=head1 NAME

Bio::EnsEMBL::Variation::Utils::Sequence - Utility functions for sequences

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::Utils::Sequence qw(ambiguity_code variation_class);

  my $alleles = 'A|C';

  print "my alleles = $alleles\n";

  my $ambig_code = ambiguity_code($alleles);

  print "my ambiguity code is $ambig_code\n";

  print "my SNP class is = variation_class($alleles)";


=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Utils::Sequence;

use Bio::EnsEMBL::Utils::Exception qw(warning);
use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&ambiguity_code &variation_class &unambiguity_code &sequence_with_ambiguity);


=head2 ambiguity_code

  Arg[1]      : string $alleles
  Example     :  use Bio::EnsEMBL::Variation::Utils::Sequence qw(ambiguity_code)
                 my $alleles = 'A|C';
                 my $ambig_code = ambiguity_code($alleles);
                print "the ambiguity code for $alleles is: ",$ambig_code;
  Description : returns the ambiguity code for a SNP allele
  ReturnType  : String
                The ambiguity code
  Exceptions  : None
  Caller      : Variation, VariationFeature

=cut

sub ambiguity_code {
    my $alleles = shift;
    my %duplicates; #hash containing all alleles to remove duplicates
    map {$duplicates{$_}++} split /[\|\/\\]/, $alleles;
    $alleles = uc( join '', sort keys %duplicates );
    my %ambig = qw(AC M ACG V ACGT N ACT H AG R AGT D AT W CG S CGT B CT Y 
GT K C C A A T T G G - -); #we will need to decide what to do with alleles like -A. Is that possible??
    return $ambig{$alleles};
}

=head2 unambiguity_code

  Arg[1]      : string $alleles
  Example     :  use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code)
                 my $ambiguity_code = 'M';
                 my $alleles = unambiguity_code($ambiguity_code);
                print "the alleles for ambiguity code $ambiguity_code is: ",$alleles;
  Description : returns the alleles for an ambiguity code
  ReturnType  : String
                The Alleles, alphabetically sorted and in capital
  Exceptions  : None
  Caller      : Variation, VariationFeature

=cut

sub unambiguity_code {
    my $ambiguity_code = shift;
   
    my %unambig = qw(M AC V ACG N ACGT H ACT R AG D AGT W AT S CG B CGT Y CT K 
GT C CC A AA T TT G GG - --); #we will need to decide what to do with alleles like -A. Is that possible??
    return $unambig{$ambiguity_code};
}


=head2 variation_class

  Arg[1]      : string $alleles
  Example     : use Bio::EnsEMBL::Variation::Utils::Sequence qw (variation_class)
                my $alleles = 'A|C';    
                my $variation_class = variation_class($alleles);
                print "the variation class for the alleles $alleles is: ",$variation_class;
  Description : return the class of the alleles according to dbSNP classification(SNP,indel,mixed,mnp...)
  ReturnType  : String. The class of the alleles
  Exceptions  : none
  Caller      : Variation, VariationFeature

=cut

sub variation_class{
    my $alleles = shift;
    return 'snp' if $alleles =~ /^[ACGT]([\|\\\/][ACGT])+$/i;
    return 'cnv' if (($alleles eq 'cnv') || ($alleles eq 'CNV'));

    my @alleles = split /[\|\/\\]/, $alleles;

    if (@alleles == 1){       
	#(HETEROZYGOUS) 1 allele
	return 'het'
    }
    elsif(@alleles == 2){
	if ((($alleles[0] =~ tr/ACTG//)== length($alleles[0]) && ($alleles[1] =~ tr/-//) == 1) || (($alleles[0] =~ tr/-//) == 1 && ($alleles[1] =~ tr/ACTG//) == length($alleles[1]))){
	    #A/- 2 alleles
	    return 'in-del'
	    }
	elsif (($alleles[0] =~ /LARGE/) || ($alleles[1] =~ /LARGE/)){
	    #(LARGEDELETION) 2 alleles
	    return 'named'
	    }
	elsif (($alleles[0] =~ tr/ACTG//) > 1 || ($alleles[1] =~ tr/ACTG//) > 1){
	    #AA/GC 2 alleles
	    return 'mnp'
	    }
	else{
	    warning("not possible to determine class for  @alleles");
	    return '';
	}
    }
    elsif(@alleles > 2){
	if ($alleles[0] =~ /\d+/){ 
	    #(CA)14/15/16/17 > 2 alleles, all of them contain the number of repetitions of the allele
	    return 'microsat'
	    }
	
	elsif ((grep {/-/} @alleles) > 0){
	    #-/A/T/TTA > 2 alleles
	    return 'mixed'
	    }
	else{
	  #  warning("not possible to determine class of alleles " . @alleles);
	    return '';
	}
    }
    else{
	warning("no alleles available ");
	return '';
    }
}

=head2 sequence_with_ambiguity

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBAdaptor $dbCore
  Arg[2]      : Bio::EnsEMBL::Variation::DBSQL::DBAdaptor $dbVar
  Arg[3]      : string $chr
  Arg[4]      : int $start
  Arg[5]      : int $end
  Arg[6]      : int $strand
  Example     : use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity)
                my $slice = sequence_with_ambiguity($dbCore,$dbVar,1,100,200);
                print "the sequence with ambiguity code for your region is: ",$slice->seq()
  Description : given a region, returns a Bio::EnsEMBL::Slice object with 
                the sequence set with ambiguity codes
  ReturnType  : Bio::EnsEMBL::Slice object
  Exceptions  : none
  Caller      : general

=cut

sub sequence_with_ambiguity{
    my ($dbCore,$dbVar,$chr,$start,$end,$strand) = @_;

    my $slice;
    if (ref($dbCore) ne 'Bio::EnsEMBL::DBSQL::DBAdaptor'){
	warning('You need to provide a Bio::EnsEMBL::DBSQL::DBAdaptor as a first argument');
	return $slice;
    }
    if (ref($dbVar) ne 'Bio::EnsEMBL::Variation::DBSQL::DBAdaptor'){
	warning('You need to provide a Bio::EnsEMBL::Variation::DBSQL::DBAdaptor object as second argument');
	return $slice;
    }
    my $slice_adaptor = $dbCore->get_SliceAdaptor();
    my $vf_adaptor = $dbVar->get_VariationFeatureAdaptor;
    $slice = $slice_adaptor->fetch_by_region('chromosome',$chr,$start,$end,$strand); #get the slice
    my $seq = $slice->seq;
    foreach my $vf (@{$vf_adaptor->fetch_all_by_Slice($slice)}){
	substr($seq,$vf->start-1,1,$vf->ambig_code);
    }
    $slice->{'seq'} = $seq;
    return $slice;
}

1;

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

@EXPORT_OK = qw(&ambiguity_code &variation_class);


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
    $alleles = uc( join '', sort split /[\|\/\\]/, $alleles );
    my %ambig = qw(AC M ACG V ACGT N ACT H AG R AGT D AT W CG S CGT B CT Y 
GT K);
    return $ambig{$alleles};
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
    my @alleles = split /[\|\/\\]/, $alleles;

#    print STDERR "\n",@alleles," and number ",scalar(@alleles),"\n";
    if (@alleles == 1){       
	#(HETEROZYGOUS) 1 allele
	return 'het'
    }
    elsif(@alleles == 2){
	if (($alleles[0] =~ tr/ACTG//) == 1 && ($alleles[1] =~ tr/ACTG//) == 1){
	    #A/T 2 alleles
	    return 'snp';
	}
	elsif ((($alleles[0] =~ tr/ACTG//)== length($alleles[0]) && ($alleles[1] =~ tr/-//) == 1) || (($alleles[0] =~ tr/-//) == 1 && ($alleles[1] =~ tr/ACTG//) == length($alleles[1]))){
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
	    warning("not possible to determine class for " . @alleles);
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
	    warning("not possible to determine class of alleles " . @alleles);
	    return '';
	}
    }
    else{
	warning("no alleles available ");
	return '';
    }
}

1;

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

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp); 
use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&ambiguity_code &variation_class &unambiguity_code &sequence_with_ambiguity &hgvs_variant_notation);


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
	
	foreach my $a(split /[\|\/\\]/, $alleles) {
		# convert Ns
		my @a = ($a eq 'N' ? qw(A C G T) : ($a));
		map {$duplicates{$_}++} @a;
	}
    $alleles = uc( join '', sort keys %duplicates );
    #my %ambig = qw(AC M ACG V ACGT N ACT H AG R AGT D AT W CG S CGT B CT Y 
#GT K C C A A T T G G - - -A -A -C -C -G -G -T -T A- A- C- C- G- G- T- T-); #for now just make e.g. 'A-' -> 'A-'
	my %ambig = qw(AC M ACG V ACGT N ACT H AG R AGT D AT W CG S CGT B CT Y GT K C C A A T T G G - -);
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
   
    #my %unambig = qw(M AC V ACG N ACGT H ACT R AG D AGT W AT S CG B CGT Y CT K 
#GT C CC A AA T TT G GG - -- -A -A -C -C -G -G -T -T A- A- C- C- G- G- T- T-); #for now just make e.g. 'A-' -> 'A-'
	my %unambig = qw(M AC V ACG N ACGT H ACT R AG D AGT W AT S CG B CGT Y CT K GT C CC A AA T TT G GG - --);
    return $unambig{$ambiguity_code};
}


=head2 variation_class

  Arg[1]      : string $alleles
  Arg[2]      : boolean $is_somatic - flag that this variation is somatic
  Example     : use Bio::EnsEMBL::Variation::Utils::Sequence qw (variation_class)
                my $alleles = 'A|C';    
                my $variation_class = variation_class($alleles);
                print "the variation class for the alleles $alleles is: ",$variation_class;
  Description : return the class of the alleles according to dbSNP classification(SNP,indel,mixed,substitution...)
  ReturnType  : String. The class of the alleles
  Exceptions  : none
  Caller      : Variation, VariationFeature

=cut

sub variation_class{
    
    my ($alleles, $is_somatic) = @_;
     
    my $class;
    
    if ($alleles =~ /^[ACGTN]([\|\\\/][ACGTN])+$/i) {
        $class = 'snp';
    }
    elsif (($alleles eq 'cnv') || ($alleles eq 'CNV')) {
        $class = 'cnv';
    }
    elsif ($alleles =~ /CNV\_PROBE/i) {
        $class = 'cnv probe';
    }
    elsif ($alleles =~ /HGMD\_MUTATION/i) {
        $class = 'hgmd_mutation';
    }
    else {
        my @alleles = split /[\|\/\\]/, $alleles;

        if (@alleles == 1) {
            #(HETEROZYGOUS) 1 allele
	        $class =  'het';
        }
        elsif(@alleles == 2) {
	       if ((($alleles[0] =~ tr/ACTGN//)== length($alleles[0]) && ($alleles[1] =~ tr/-//) == 1) || 
	           (($alleles[0] =~ tr/-//) == 1 && ($alleles[1] =~ tr/ACTGN//) == length($alleles[1])) ){
	           #A/- 2 alleles
	           $class =  'in-del'
	       }
	       elsif (($alleles[0] =~ /LARGE|INS|DEL/) || ($alleles[1] =~ /LARGE|INS|DEL/)){
	           #(LARGEDELETION) 2 alleles
	           $class = 'named'
	       }
	       elsif (($alleles[0] =~ tr/ACTG//) > 1 || ($alleles[1] =~ tr/ACTG//) > 1){
	           #AA/GC 2 alleles
	           $class = 'substitution'
	       }
	       else {
	           warning("not possible to determine class for  @alleles");
	           $class = '';
	       }
        }
        elsif (@alleles > 2) {
	       
	       if ($alleles[0] =~ /\d+/) { 
	           #(CA)14/15/16/17 > 2 alleles, all of them contain the number of repetitions of the allele
	           $class = 'microsat'
	       }
	
	       elsif ((grep {/-/} @alleles) > 0) {
	           #-/A/T/TTA > 2 alleles
	           $class = 'mixed'
	       }
	       else {
	           #  warning("not possible to determine class of alleles " . @alleles);
	           $class = '';
	       }
        }
        else{
	       warning("no alleles available ");
	       $class = '';
        }
    }
    
    if ($is_somatic) {
        if ($class eq '') {
            # for undetermined classes just call it somatic
            $class = 'somatic';
        }
        else {       
            # somatic mutations aren't polymorphisms, so change SNPs to SNVs
            $class = 'snv' if $class eq 'snp'; 
 
            # and prefix the class with 'somatic' 
            $class = 'somatic_'.$class; 
        }
    }
    
    return $class;
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

=head2 hgvs_variant_notation

  Arg[1]      : string $alt_allele
  Arg[2]      : string $ref_sequence
  Arg[3]      : int $ref_start
  Arg[4]      : int $ref_end
  Arg[5]      : int $display_start (optional)
  Arg[6]      : int $display_end (optional)
  
  Example     : use Bio::EnsEMBL::Variation::Utils::Sequence qw (hgvs_variant_notation)
		my $alt_allele = 'A';
		my $ref_sequence = 'CCGTGATGTGC';
		my $ref_start = 4;
		my $ref_end = 4;
		my $ref_name = 'test_seq';
		my $ref_type = 'g';
		my $notation = hgvs_variant_notation($alt_allele,$ref_sequence,$ref_start,$ref_end);
                print "HGVS notation of your variant: $ref_name\:$ref_type\." . $notation->{'hgvs'};
		
  Description : Given an allele, a reference sequence and position of variant, returns a reference to a hash containing metadata and a 
		string with HGVS notation of the variant. Returns undef if reference and variant alleles are identical.
		The optional display_start and display_end, if specified, will be used in the notation instead of the ref_start and ref_end.
		This can be useful, e.g. if we want coordinates relative to chromosome but don't want to pass the entire chromosome sequence
		into the subroutine.
		The data fields in the returned hash are:
		'start'	-> Displayed start position of variant
		'end' -> Displayed end position of variant
		'ref' -> Reference allele
		'alt' -> Alternative allele
		'type' -> The variant class, e.g. ins, inv, >, delins
		'hgvs' -> A string with HGVS notation
  ReturnType  : reference to a hash
  Exceptions  : If the length of the interval to be displayed is different from the length of the reference allele
  Caller      : general

=cut
sub hgvs_variant_notation {
    my $alt_allele = shift;
    my $ref_sequence = shift;
    my $ref_start = shift;
    my $ref_end = shift;
    my $display_start = shift;
    my $display_end = shift;
    
    # If display_start and display_end were not specified, use ref_start and ref_end
    $display_start ||= $ref_start;
    $display_end ||= $ref_end;
    
    #ÊThrow an exception if the lengths of the display interval and reference interval are different
    throw("The coordinate interval for display is of different length than for the reference allele") if (($display_end - $display_start) != ($ref_end - $ref_start));
    
    # Length of the reference allele. Negative lengths make no sense
    my $ref_length = ($ref_end - $ref_start + 1);
    if ($ref_length < 0) {
	$ref_length = 0;
    }
    
    # Remove any gap characters in the alt allele
    $alt_allele =~ s/\-//g;
    
    # Length of alternative allele
    my $alt_length = length($alt_allele);
    
    # Get the reference allele
    my $ref_allele = substr($ref_sequence,($ref_start-1),$ref_length);
    
    # Check that the alleles are different, otherwise return undef
    return undef unless ($ref_allele ne $alt_allele);
    
    # Store the notation in a hash that will be returned
    my %notation;
    $notation{'start'} = $display_start;
    $notation{'end'} = $display_end;
    $notation{'ref'} = $ref_allele;
    $notation{'alt'} = $alt_allele;
    
    # The simplest case is a deletion
    if (!$alt_length) {
	$notation{'type'} = 'del';
        
        # Return the notation
        return \%notation;
    }
    
    # Another case is if the allele lengths are equal
    if ($ref_length == $alt_length) {
        
	# If length is 1 it's a single substitution
        if ($ref_length == 1) {
	    $notation{'type'} = '>';
	    return \%notation;
        }
        
	# Check if it's an inversion
        my $rev_ref = $ref_allele;
        reverse_comp(\$rev_ref);
        if ($alt_allele eq $rev_ref) {
	    $notation{'type'} = 'inv';
	    return \%notation;
        }
        
	$notation{'type'} = 'delins';
	
        return \%notation;
    }
    
    # If this is an insertion, we should check if the preceeding reference nucleotides match the insertion. In that case it should be annotated as a multiplication.
    if (!$ref_length) {
    
        # Get the same number of nucleotides preceding the insertion as the length of the insertion
        my $prev_str = substr($ref_sequence,($ref_end-$alt_length),$alt_length);
        
        # If they match, this is a duplication
        if ($prev_str eq $alt_allele) {
	    $notation{'start'} = ($display_end - $alt_length + 1);
	    $notation{'type'} = 'dup';
	    $notation{'ref'} = $prev_str;
            # Return the notation
	    return \%notation;
        }
        
        # If they didn't match it's a plain insertion
	$notation{'start'} = $display_end;
	$notation{'end'} = $display_start;
	$notation{'type'} = 'ins';
        
        return \%notation;
    }
    
    # Otherwise, the reference and allele are of different lengths. By default, this is a delins but
    # we need to check if the alt allele is a multiplication of the reference
    # Check if the length of the alt allele is a multiple of the reference allele
    if ($alt_length%$ref_length == 0) {
        my $multiple = ($alt_length / $ref_length);
        if ($alt_allele eq ($ref_allele x $multiple)) {
            if ($multiple == 2) {
		$notation{'type'} = 'dup';
            }
            else {
		$notation{'type'} = '[' . $multiple . ']';
            }
	    return \%notation;
        }
    }
    
    # Else, it's gotta be a delins
    $notation{'type'} = 'delins';
    
    return \%notation;
}

1;

=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

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
use Bio::EnsEMBL::Variation::Utils::Constants qw(:SO_class_terms);
use Bio::EnsEMBL::Utils::Scalar qw(wrap_array);
use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(
    &ambiguity_code 
    &variation_class 
    &unambiguity_code 
    &sequence_with_ambiguity 
    &hgvs_variant_notation 
    &format_hgvs_string
    &get_hgvs_alleles
    &SO_variation_class 
    &align_seqs 
    &strain_ambiguity_code
    &get_all_validation_states
    &get_validation_code
    &add_validation_state
    &revcomp_tandem
    %EVIDENCE_VALUES
);

# List of validation states. Order must match that of set in database
our @VALIDATION_STATES = qw(cluster freq submitter doublehit hapmap 1000Genome failed precious);

our %EVIDENCE_VALUES  = qw(Multiple_observations 1 Frequency 2 HapMap 2 1000Genomes 2 Cited 3 ESP 4);

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

=head2 strain_ambiguity_code

  Arg[1]      : string $alleles (separated by "/", "\" or "|")
  Example     :  use Bio::EnsEMBL::Variation::Utils::Sequence qw(strain_ambiguity_code)
                 my $alleles = 'A|C';
                 my $ambig_code = strain_ambiguity_code($alleles);
                print "the ambiguity code for $alleles is: ",$ambig_code;
  Description : returns the ambiguity code for a strain genotype
  ReturnType  : String
  Exceptions  : None
  Caller      : AlleleFeatureAdaptor

=cut

sub strain_ambiguity_code {
    my $alleles = shift;
	
	# return normal ambiguity code for a SNP
	return ambiguity_code($alleles) if($alleles =~ /^[ACGT][\|\/\\][ACGT]$/);
	
	# get alleles
	my ($a1, $a2) = split /[\|\/\\]/, $alleles;
	
	# pad
	if(length($a1) > length($a2)) {
		$a2 .= '-' x (length($a1) - length($a2));
	}
	else {
		$a1 .= '-' x (length($a2) - length($a1));
	}
	
	# build ambiguity code base by base
	my $ambig = '';
	
	for my $i(0..(length($a1) - 1)) {
		my $b1 = substr($a1, $i, 1);
		my $b2 = substr($a2, $i, 1);
		
		# -/- = -
		if($b1 eq '-' && $b2 eq '-') {
			$ambig .= '-';
		}
		
		# G/- = g
		elsif($b1 eq '-') {
			$ambig .= lc($b2);
		}
		
		# -/G = g
		elsif($b2 eq '-') {
			$ambig .= lc($b1);
		}
		
		# A/G = R
		else {
			$ambig .= ambiguity_code($b1.'|'.$b2);
		}
	}
	
	return $ambig;
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

=head2 SO_variation_class

  Arg[1]      : string $alleles - a slash ()'/') separated list of alleles, the first allele is 
                assumed to be the reference unless the $ref_correct argument is false
  Arg[2]      : boolean $ref_correct - flags that the first allele is not known to be the 
                reference sequence (so we can't call insertions or deletions and have to
                resort to 'indel')
  Example     : use Bio::EnsEMBL::Variation::Utils::Sequence qw (SO_variation_class)
                my $alleles = 'A/C';    
                my $SO_term = SO_variation_class($alleles);
                print "the SO term for the alleles $alleles is: ",$SO_term;
  Description : return the SO term for the class of the alleles
  ReturnType  : String. The SO term for the class of the alleles
  Exceptions  : none
  Caller      : Variation, VariationFeature

=cut

sub SO_variation_class {
    
    my $alleles     = shift;
    my $ref_correct = shift;
    
    $ref_correct = 1 unless defined $ref_correct;

    my $allele_class = '[A-Z]';

    # default to sequence_alteration
    my $class = SO_TERM_SEQUENCE_ALTERATION;

    if ($alleles =~ /^$allele_class(\/$allele_class)+$/) {
        # A/T, A/T/G
        $class = SO_TERM_SNV;
    }
    elsif ($alleles =~ /^$allele_class+(\/$allele_class+)+$/) {
        # AA/TTT
        $class = SO_TERM_SUBSTITUTION;
    }
    elsif ($alleles =~ /\)\d+/) {
        # (CAG)8/(CAG)9
        $class = SO_TERM_TANDEM_REPEAT;
    }
    else {
        my @alleles = split /\//, $alleles;

        if (@alleles > 1) {

            my $ref = shift @alleles;

            if ($ref eq '-') {

                if (@alleles == 1 && $alleles[0] =~ /DEL/) {
                    # -/(LARGEDELETION) (rather oddly!)
                    $class = $ref_correct ? SO_TERM_DELETION : SO_TERM_INDEL;
                }

                unless (grep { $_ !~ /^$allele_class+$|INS/ } @alleles) {
                    # -/ATT, -/(LARGEINSERTION)
                    $class = $ref_correct ? SO_TERM_INSERTION : SO_TERM_INDEL;
                }

                # else must be mixed insertion and deletion, so just called sequence_alteration
            }
            elsif ($ref =~ /^$allele_class+$/) {
                unless (grep { $_ !~ /-|DEL/ } @alleles) {
                    # A/-, A/(LARGEDELETION)
                    $class = $ref_correct ? SO_TERM_DELETION : SO_TERM_INDEL;
                }
            }
            elsif ($ref =~ /DEL/) {
                unless (grep { $_ !~ /-/ } @alleles) {
                    # (LARGEDELETION)/-, (2345 BP DELETION)/-
                    $class = $ref_correct ? SO_TERM_DELETION : SO_TERM_INDEL;
                }
            }
        }
        elsif (@alleles == 1) {
            if ($alleles[0] =~ /INS/) {
                # (LARGEINSERTION)
                $class = $ref_correct ? SO_TERM_INSERTION : SO_TERM_INDEL;
            }
            elsif($alleles[0] =~ /DEL/) {
                # (308 BP DELETION)
                $class = $ref_correct ? SO_TERM_DELETION : SO_TERM_INDEL;
            }
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
    my $var_name  = shift;
    
    # If display_start and display_end were not specified, use ref_start and ref_end
    $display_start ||= $ref_start;
    $display_end ||= $ref_end;
    
    #Throw an exception if the lengths of the display interval and reference interval are different
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
    if($ref_allele eq $alt_allele){
	warn "\nError in HGVS calculation for $var_name: alt allele ($alt_allele) is the same as the reference allele ($ref_allele) - potential strand or allele ordering problem - skipping\n";
	return undef ;
    }
    
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


=head2 format_hgvs_string

  Arg[1]      : string reference sequence name
  Arg[2]      : string strand
  Arg[3]      : hash of hgvs information
  Example     : 
  Description : Creates HGVS formatted string from input hash                
  ReturnType  : string in HGVS format
  Exceptions  : 
  Caller      : 

=cut

sub format_hgvs_string{
    ##### generic formatting routine for genomic and coding HGVS names
 
    my $hgvs_notation = shift;
  
    ### all start with refseq name & numbering type
    $hgvs_notation->{'hgvs'} = $hgvs_notation->{'ref_name'} . ":" . $hgvs_notation->{'numbering'} . ".";    

    my $coordinates;
    #### if single base event, list position only once
    if($hgvs_notation->{'start'} eq $hgvs_notation->{'end'}){
	$coordinates =  $hgvs_notation->{'start'};
    }
    else{
	$coordinates = $hgvs_notation->{'start'} . "_" . $hgvs_notation->{'end'};
    }

    ##### format rest of string according to type

    if($hgvs_notation->{'type'} eq 'del' ||  $hgvs_notation->{'type'} eq 'inv' || $hgvs_notation->{'type'} eq 'dup'){
	### inversion of reference bases => list ref not alt
	### deletion  of reference bases => list ref lost
	### duplication  of reference bases (eg ref = GAAA alt = GAAAGAAA) => list duplicated ref (dupGAAA)
	$hgvs_notation->{'hgvs'} .= $coordinates . $hgvs_notation->{'type'} . $hgvs_notation->{'ref'};      
    }

    elsif( $hgvs_notation->{'type'} eq '>'){
	### substitution - list both alleles
	$hgvs_notation->{'hgvs'} .= $hgvs_notation->{'start'} . $hgvs_notation->{'ref'} . $hgvs_notation->{'type'} . $hgvs_notation->{'alt'};
    }

    elsif( $hgvs_notation->{'type'} eq 'delins'){
	$hgvs_notation->{'hgvs'} .= $coordinates . 'del' . $hgvs_notation->{'ref'} . 'ins' . $hgvs_notation->{'alt'};
    }   

    elsif($hgvs_notation->{'type'} eq 'ins'){
	## reference not listed
	$hgvs_notation->{'hgvs'} .= $coordinates . $hgvs_notation->{'type'} . $hgvs_notation->{'alt'};
    }

    elsif($hgvs_notation->{'type'} =~ /\[\d+\]/){
	#### insertion described by string and number
	$hgvs_notation->{'hgvs'} .= $coordinates  . $hgvs_notation->{'type'} . $hgvs_notation->{'ref'};
    }

    else{
	warn "PROBLEM with generic HGVS formatter - type = ". $hgvs_notation->{'type'} ."\n";
    }
  
    return $hgvs_notation->{'hgvs'};

}

=head2 get_hgvs_alleles

  Arg[1]      : string - HGVS name for a variant
  Example     : ACTA1:c.1031delG
  Description : extracts reference and alternate alleles for an hgvs string
  ReturnType  : array of alleles (ref_allele_string , alt_allele_string)
  Exceptions  : if the type of variant is not recognisable from the HGVS name
  Caller      : 

=cut

sub get_hgvs_alleles{
    
    #### extract ref and alt alleles where possible from HGVS g/c/n string

  my ( $hgvs) = shift;

  my ($reference,$type,$description) = $hgvs =~ m/^([^\:]+)\:.*?([cgmnrp]?)\.?(.*?[\*\-0-9]+.*)$/i;
  
  ## remove protein no-change part of string if present
  $description =~ s/\(p\.=\)$//;  

  my ($ref_allele, $alt_allele) ;
    
  ### A single nt substitution, reference and alternative alleles are required
  if ($description =~ m/>/) {
    ($ref_allele,$alt_allele) = $description =~ m/([A-Z]+)>([A-Z]+)$/i;
  }
    
  # delins, the reference allele is optional
  elsif ($description =~ m/del.*ins/i) {
    ($ref_allele,$alt_allele) = $description =~ m/del(.*?)ins([A-Z]+)$/i;          
  }
    
  # A deletion, the reference allele is optional
  elsif ($description =~ m/del/i) {
    ($ref_allele) = $description =~ m/del([A-Z]*)$/i; 
    $alt_allele = '-';
  }
    
  # A duplication, the reference allele is optional
  elsif ($description =~ m/dup/i) {
     $ref_allele ="-";
     ($alt_allele) = $description =~ m/dup([A-Z]*)$/i;
  }
    
  # An inversion, the reference allele is optional
  elsif ($description =~ m/inv/i) {
    ($ref_allele) = $description =~ m/inv([A-Z]*)$/i;
    $alt_allele = $ref_allele;
    reverse_comp(\$alt_allele);
  }
    
  # An insertion, 
  elsif ($description =~ m/ins/i) {
    ($alt_allele) = $description =~ m/ins([A-Z]*)$/i;
    $ref_allele = '-';
  }
  ## A simple repeat (eg. ENST00000522587.1:c.-310+750[13]A => alt AAAAAAAAAAAAA)
  elsif ($description =~ m/\[/i) {    
  
    my ($number, $string) = $description =~ m/\[(\d+)\]([A-Z]*)$/i; 
    foreach my $n(1..$number){ $alt_allele .= $string;}
    $ref_allele = $string;
  }
  else {
    throw ("The variant class for HGVS notation '$hgvs' is unknown or could not be correctly recognized");
  }
  return ($ref_allele, $alt_allele) ;
}


=head2 align_seqs

  Arg[1]      : string $seq1
  Arg[2]      : string $seq2
  Example     : my $aligned_seqs = align_seqs($seq1, $seq2);
  Description : Does a simple NW align of two sequence strings. Best used on
                short (<1000bp) sequences, otherwise runtime will be long
  ReturnType  : arrayref to a pair of strings
  Exceptions  : none
  Caller      : web flanking sequence display

=cut

sub align_seqs {
	my $seq1 = shift;
	my $seq2 = shift;

	# align parameters
	my $match    = 10;
	my $mismatch = -10;
	my $gep      = -10;
	
	# split sequences into arrays
	my @split1 = split //, $seq1;
	my @split2 = split //, $seq2;
	
	# evaluate substitutions
	my $len1 = length($seq1);
	my $len2 = length($seq2);
	
	my (@smat, @tb);
	
	for (my $i=0; $i<=$len1; $i++) {
		$smat[$i][0] = $i * $gep;
		$tb[$i][0] = 1;
	}
	for (my $j=0; $j<=$len2; $j++) {
		$smat[0][$j] = $j * $gep;
		$tb[0][$j] = -1;
	}
	
	my ($s, $sub, $del, $ins);
	
	for (my $i=1; $i<=$len1; $i++) {
		for (my $j=1; $j<=$len2; $j++)	{
			
			# calculate score
			if($split1[$i-1] eq $split2[$j-1]) {
				$s = $match;
			}
			else {
				$s = $mismatch;
			}
			
			$sub = $smat[$i-1][$j-1] + $s;
			$del = $smat[$i][$j-1] + $gep;
			$ins = $smat[$i-1][$j] + $gep;
			
			if($sub > $del && $sub > $ins) {
				$smat[$i][$j] = $sub;
				$tb[$i][$j] = 0;
			}
			elsif($del > $ins) {
				$smat[$i][$j] = $del;
				$tb[$i][$j] = -1;
			}
			else {
				$smat[$i][$j] = $ins;
				$tb[$i][$j] = 1;
			}
		}
	}
	
	
	my $i = $len1;
	my $j = $len2;
	my $aln_len = 0;
	my (@aln1, @aln2);
	
	while(!($i == 0 && $j == 0)) {
		if($tb[$i][$j] == 0) {
			$aln1[$aln_len] = $split1[--$i];
			$aln2[$aln_len] = $split2[--$j];
		}
		elsif($tb[$i][$j] == -1) {
			$aln1[$aln_len] = '-';
			$aln2[$aln_len] = $split2[--$j];
		}
		elsif($tb[$i][$j] == 1) {
			$aln1[$aln_len] = $split1[--$i];
			$aln2[$aln_len] = '-';
		}
		
		$aln_len++;
	}
	
	return [(join "", reverse @aln1), (join "", reverse @aln2)];
}


=head2 array_to_bitval

  Arg[1]      : arrayref $arr
  Arg[2]      : arrayref $ref
  Example     : my $bitval = array_to_bitval(['hapmap','precious'],['cluster','freq','submitter','doublehit','hapmap','1000Genome','failed','precious']);
  Description : Takes a reference to an array as input and return a bit value representing the 
                combination of elements from a reference array. c.f. the SET datatype in MySQL 
  ReturnType  : bitvalue that represents the combination of elements in the reference array specified in the given array
  Exceptions  : none
  Caller      : get_validation_code

=cut

sub array_to_bitval {
    my $arr = shift;
    my $ref = shift;
    
    #Ensure that we have array references
    $arr = wrap_array($arr);
    $ref = wrap_array($ref);
    
    #Turn the reference array into a hash, the values will correspond to 2 raised to the power of the position in the array
    my $i=0;
    my %ref_hash = map {lc($_) => $i++;} @{$ref}; 
    
    #Set the bitval
    my $bitval = 0;
    foreach my $a (@{$arr}) {
        
        my $pos = $ref_hash{lc($a)};
        if (defined($pos)) {
            $bitval |= 2**$pos;    
        }
        # Warn if the element is not present in the reference array
        else {
            warning("$a is not a recognised element. Recognised elements are: " . join(",",@{$ref}));
        }
    }
    
    return $bitval;      
}

=head2 bitval_to_array

  Arg [1]    : int $bitval
  Arg [2]    : arrayref $ref
  Example    : my $arr = bitval_to_array(6,['cluster','freq','submitter','doublehit','hapmap','1000Genome','failed','precious']);
             : print join(",",@{$arr}); #Will print 'freq,submitter'
  Description: Returns an array with the combination of elements from the reference array specified by the supplied bitvalue. 
               c.f. the SET datatype in MySQL
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : get_all_validation_states

=cut

sub bitval_to_array {
    my $bitval = shift || 0;
    my $ref = shift;

    #Ensure that we have array references
    $ref = wrap_array($ref);
    
    # convert the bit value into an ordered array
    my @arr;
    for (my $i = 0; $i < @{$ref}; $i++) {
        push(@arr,$ref->[$i]) if ((1 << $i) & $bitval);
    }

    return \@arr;
}


=head2 add_validation_state

  Arg [1]    : string $state
  Example    : add_validation_state('cluster');
  Description: Adds a validation state to this variation.
  Returntype : none
  Exceptions : warning if validation state is not a recognised type
  Caller     : general
  Status     : At Risk

=cut

sub add_validation_state {
  my $obj = shift;
  my $state = shift;
  
  #Get the bitvalue for the new state 
  my $newbit = get_validation_code($state) || 0;
  
  #Bit-add it to the current validation_code
  my $oldbit = $obj->{'validation_code'} || 0;
  $newbit |= $oldbit;
  
  # Set the validation_code
  $obj->{'validation_code'} = $newbit;
  
  return;
}

=head2 get_all_validation_states

  Arg [1]    : int $bitval
  Example    : my @vstates = @{get_all_validation_states($var->{'validation_code'})};
  Description: Retrieves all validation states for a specified bit value.
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : general

=cut

sub get_all_validation_states {
    return bitval_to_array(shift,\@VALIDATION_STATES);
}

=head2 get_validation_code

  Arg [1]    : arrayref $validation_status
  Example    : $var->{'validation_code'} = get_validation_code(['submitter','precious']);
  Description: Retrieves the bit value for a combination of validation statuses.
  Returntype : int
  Exceptions : none
  Caller     : Variation::new

=cut

sub get_validation_code {
    return array_to_bitval(shift,\@VALIDATION_STATES);
}



=head2 revcomp_tandem

  Arg [1]    : string eg. "(AC)17/(AC)19"
  Example    : $rc_string = revcomp_tandem("(AC)17/(AC)19");
  Description: reverse compliments tandem repeat descriptions 
  Returntype : string eg. "(GT)17/(GT)19"
  Exceptions : none
  Caller     : Variation::Pipeline::VariantQC::VariantQC
  Status     : At Risk

=cut

sub revcomp_tandem{

    
  my $allele_string = shift;

  my $new_allele_string;

  my @parts = split/\//, $allele_string;
  
  foreach my $part (@parts){
  
    if( $part =~/\d+$/){

      my ($seq, $num ) = split/\)/,$part;
      $seq =~ s/\(//g;     
      reverse_comp(\$seq);

      $new_allele_string .= "(" . $seq .")" . $num . "/";
    }
    else{
      reverse_comp(\$part);

      $new_allele_string .= $part. "/";;
    }
  }
  $new_allele_string =~ s/\/$//;
  
  return $new_allele_string ;
}




1;

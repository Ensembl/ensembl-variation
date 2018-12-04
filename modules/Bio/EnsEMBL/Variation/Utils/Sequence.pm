=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

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

use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp); 
use Bio::EnsEMBL::Variation::Utils::Constants qw(:SO_class_terms);
use Bio::EnsEMBL::Utils::Scalar qw(wrap_array assert_ref);
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
    &get_3prime_seq_offset
    &SO_variation_class 
    &align_seqs 
    &strain_ambiguity_code
    &revcomp_tandem
    &get_matched_variant_alleles
    &trim_sequences
    &raw_freqs_from_gts
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
	return ambiguity_code($alleles) if($alleles =~ /^[ACGT]([\|\/\\][ACGT])*$/);
	
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

                elsif ($alleles =~ /^$allele_class+(\/$allele_class+)+$/) {
                    # AA/TT   => SO_TERM_SUBSTITUTION
                    # AA/TTTT => SO_TERM_INDEL
                    my $same_size = 1;
                    foreach (my $n =0; $n< scalar(@alleles); $n++){
                        $same_size = 0 unless length($alleles[$n]) eq length($ref);
                    }
                    $class = $same_size == 1 ? SO_TERM_SUBSTITUTION : SO_TERM_INDEL;
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
    #warn "\nError in HGVS calculation for $var_name: alt allele ($alt_allele) is the same as the reference allele ($ref_allele) - potential strand or allele ordering problem - skipping\n";
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
  ReturnType  : String in HGVS format
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

  if( $hgvs_notation->{'type'} eq '>' || ($hgvs_notation->{'type'} eq 'inv' && length($hgvs_notation->{'ref'}) ==1)  ){
    ### substitution - list both alleles
    $hgvs_notation->{'hgvs'} .= $hgvs_notation->{'start'} . $hgvs_notation->{'ref'} . '>' . $hgvs_notation->{'alt'};
  }

  elsif($hgvs_notation->{'type'} eq 'del' ||  $hgvs_notation->{'type'} eq 'inv' || $hgvs_notation->{'type'} eq 'dup'){
    ### inversion of reference bases => list ref not alt
    ### deletion  of reference bases => list ref lost
    ### duplication  of reference bases (eg ref = GAAA alt = GAAAGAAA) => list duplicated ref (dupGAAA)
    $hgvs_notation->{'hgvs'} .= $coordinates . $hgvs_notation->{'type'};
  }

  elsif( $hgvs_notation->{'type'} eq 'delins'){
    $hgvs_notation->{'hgvs'} .= $coordinates . 'delins' . $hgvs_notation->{'alt'};
  }   

  elsif($hgvs_notation->{'type'} eq 'ins'){
    ## reference not listed
    $hgvs_notation->{'hgvs'} .= $coordinates . $hgvs_notation->{'type'} . $hgvs_notation->{'alt'};
  }

  elsif($hgvs_notation->{'type'} =~ /\[\d+\]/){
    #### insertion described by string and number
    $hgvs_notation->{'hgvs'} .= $coordinates . $hgvs_notation->{'type'};
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
    
  # A deletion, the alleles should be defined by position alone
  # Previous HGVS versions listed the deleted sequence - captured here
  elsif ($description =~ m/del/i) {
    ($ref_allele) = $description =~ m/del([A-Z]*)$/i; 
    $alt_allele = '-';
  }
    
  # A duplication, the alleles should be defined by position alone
  # Previous HGVS versions listed the duplicated sequence - captured here
  elsif ($description =~ m/dup/i) {
    $ref_allele ="-";
    ($alt_allele) = $description =~ m/dup([A-Z]*)$/i;
  }
    
  # An inversion, the the alleles should be defined by position alone
  # Previous HGVS versions listed the duplicated sequence - captured here
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
  # no change
  elsif ($description =~ m/\=/i) {
    ($ref_allele) = $description =~ m/([A-Z]*)\=$/i;
     $alt_allele = $ref_allele;
  }

  else {
    throw ("The variant class for HGVS notation '$hgvs' is unknown or could not be correctly recognized");
  }
  ## get rid of: c.748_753del6insGGCCG
  undef $ref_allele  if defined $ref_allele && $ref_allele =~/\d+/;

  return ($ref_allele, $alt_allele) ;
}

=head2 get_3prime_seq_offset
  Arg[1]     : allele sequence
  Arg[2]     : downstream flank
  Description: Compare an allele to its 3' sequence to define the most 3'
               description of the change for HGVS
  Returntype : string or undef if this allele is not in the
  Exceptions : none
  Status     : Experimental
=cut
sub get_3prime_seq_offset{

  my $seq_to_check  = shift;
  my $post_seq      = shift;

  my $offset = 0;

  return ($seq_to_check, $offset)  unless defined  $post_seq;

  ## get length of pattern to check 
  my $check_length = length($post_seq) - length($seq_to_check);

  ## move along sequence after deletion looking for match to start of deletion
  for (my $n = 0; $n<= $check_length; $n++ ){

    ## check each position in deletion/ following seq for match
    my $check_next_al  = substr( $seq_to_check, 0, 1);
    my $check_next_ref = substr( $post_seq, $n, 1);

    ## stop if the sequences differ
    last if $check_next_al ne $check_next_ref;

    ## move position of deletion along
    $offset++;

    ## modify deleted sequence - remove start & append to end
    $seq_to_check = substr($seq_to_check,1);
    $seq_to_check .= $check_next_ref;
  }

  return ($seq_to_check, $offset);
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


=head2 trim_sequences

  Arg[1]      : string $ref
  Arg[2]      : string $alt
  Arg[3]      : (optional) int $start
  Arg[4]      : (optional) int $end
  Arg[5]      : (optional) bool $empty_to_dash
  Arg[6]      : (optional) bool $end_first
  Example     : my ($new_ref, $new_alt, $new_start) = @{trim_sequences($ref, $alt, $start)}
  Description : Takes a pair of reference and alternate sequences and trims common sequence
                from the start and then the end to give the minimal pair of alleles,
                adjusting $start and $end coordinates accordingly.
                
                $start is optional and will be given initial value of 0 if not supplied
                $end is optional and will be set to $start + (length($ref) - 1) if not supplied

                $ref or $alt may end up as an empty string - to comply with Ensembl Variation
                convention you may wish to set $empty_to_dash to a true value to
                return "-" instead of an empty string.

                By default common sequence is trimmed from the start of each sequence first;
                to trim from the end first, set $end_first to a true value.

                A boolean flag indicating if any change was made is returned.
  ReturnType  : arrayref:
                [
                  string $new_ref,
                  string $new_alt,
                  int $new_start,
                  int $new_end,
                  bool $changed
                ]
  Exceptions  : throws if no $ref and/or $alt supplied
  Caller      : VEP, general

=cut

sub trim_sequences {
  my ($ref, $alt, $start, $end, $empty_to_dash, $end_first) = @_;

  throw("Missing reference or alternate sequence") unless $ref && $alt;

  $start ||= 0;
  $end ||= $start + (length($ref) - 1);

  my $changed = 0;

  if($end_first) {
    # trim from right
    while($ref && $alt && substr($ref, -1, 1) eq substr($alt, -1, 1)) {
      $ref = substr($ref, 0, length($ref) - 1);
      $alt = substr($alt, 0, length($alt) - 1);
      $end--;
      $changed = 1;
    }

    # trim from left
    while($ref && $alt && substr($ref, 0, 1) eq substr($alt, 0, 1)) {
      $ref = substr($ref, 1);
      $alt = substr($alt, 1);
      $start++;
      $changed = 1;
    }
  }

  else {
    # trim from left
    while($ref && $alt && substr($ref, 0, 1) eq substr($alt, 0, 1)) {
      $ref = substr($ref, 1);
      $alt = substr($alt, 1);
      $start++;
      $changed = 1;
    }

    # trim from right
    while($ref && $alt && substr($ref, -1, 1) eq substr($alt, -1, 1)) {
      $ref = substr($ref, 0, length($ref) - 1);
      $alt = substr($alt, 0, length($alt) - 1);
      $end--;
      $changed = 1;
    }
  }

  if($empty_to_dash) {
    $ref = '-' if $ref eq '';
    $alt = '-' if $alt eq '';
  }

  return [$ref, $alt, $start, $end, $changed];
}


=head2 get_matched_variant_alleles

  Arg[1]      : hashref $var_a (see below for expected structure)
  Arg[2]      : hashref $var_b
  Example     : my @matches = @{get_matched_variant_alleles($var_a, $var_b)}
  Description : Takes a pair of simplified variant hashes and finds alternate
                alleles that match between the two.

                $var = {

                  # describe alleles as ref and alts
                  ref  => 'A',
                  alts => ['G', 'T'],

                  # OR as "/"-separated allele string, ref first (marginally less efficient)
                  allele_string => 'A/G/T',

                  pos => 12345,
                  strand => 1, # (optional, defaults to 1)
                }

                Within each variant, each reference (REF) and alternate (ALT)
                allele pair is considered separately, with identical sequence
                being trimmed from each end using trim_sequences(). Both
                directions are tested (trim from start then end, then trim from
                end then start).

                Matching pairs of alleles are returned as a hashref indicating
                the matched allele as it originally appeared as well as the index
                position of the matched allele in the array of ALTs:
  ReturnType  : arrayref of hashrefs:
                [
                  {
                    a_allele => $matched_allele_from_var_a,
                    a_index  => $index_of_matched_allele_in_b_alts,
                    b_allele => $matched_allele_from_var_b,
                    b_index  => $index_of_matched_allele_in_b_alts,
                  }
                ]
  Exceptions  : throws if
                - $var_a or $var_b are not hashrefs
                - no ref, alt and/or pos key supplied for either $var_a or $var_b
  Caller      : VEP, VCFCollection, DumpVEP pipeline

=cut

sub get_matched_variant_alleles {
  my ($a, $b) = @_;

  assert_ref($a, 'HASH');
  assert_ref($b, 'HASH');

  # convert allele_string key
  foreach my $var($a, $b) {
    if(my $as = $var->{allele_string}) {
      my @alleles = split('/', $as);
      $var->{ref}  ||= shift @alleles;
      $var->{alts}   = \@alleles unless exists($var->{alts});
    }
  }

  # check ref
  throw("Missing ref key in first variant") unless exists($a->{ref});
  throw("Missing ref key in second variant") unless exists($b->{ref});

  # check alts
  $a->{alts} ||= [$a->{alt}] if defined($a->{alt});
  $b->{alts} ||= [$b->{alt}] if defined($b->{alt});
  throw("Missing alt or alts key in first variant") unless exists($a->{alts});
  throw("Missing alt or alts key in second variant") unless exists($b->{alts});

  # check pos
  throw("Missing pos key in first variant") unless $a->{pos};
  throw("Missing pos key in second variant") unless $b->{pos};

  # munge in strand
  $a->{strand} = 1 unless exists($a->{strand});
  $b->{strand} = 1 unless exists($b->{strand});

  # reverse ref of $a if required
  # do it this way as typical use case will be $a is a VF and $b is from a VCF (always fwd)
  my $a_ref = $a->{ref};
  reverse_comp(\$a_ref) if $a->{strand} != $b->{strand};

  # we might need to minimise pairs of alleles if user has input weirdness
  # so we're going to create keys for each minimised ref/alt pair
  # also trim both ways in case?
  my (%minimised_a_alleles, %a_indexes);
  my $i = 0;

  foreach my $orig_a_alt(@{$a->{alts}}) {

    # we store the index position of each alt so we can look it up and return it later
    $a_indexes{$orig_a_alt} = $i++;

    my $rev_alt = $orig_a_alt;
    reverse_comp(\$rev_alt) if $a->{strand} != $b->{strand};

    # only need to trim both directions if length of either allele > 1
    foreach my $direction(@{_get_trim_directions($a_ref, $orig_a_alt)}) {
      my $pos = $a->{pos};
      my $ref = $a_ref;
      my $alt = $rev_alt;
      ($ref, $alt, $pos) = @{trim_sequences($ref, $alt, $pos, undef, 1, $direction)};

      # store the original alt as when we come to report results
      # we need to match back against the alt as it was in the user input
      $minimised_a_alleles{"$ref\_$alt\_$pos"} = $orig_a_alt;
    }
  }

  # use these as backups as we might modify them
  my $orig_b_pos = $b->{pos};
  my $orig_b_ref = $b->{ref};

  # we're going to return a list of hashrefs containing the matched allele indexes and alleles
  my @matches;

  $i = 0;
  foreach my $orig_b_alt(@{$b->{alts}}) {

    # we're going to minimise $b's alleles one by one and compare the pos and alleles to $a
    # only need to trim both directions if length of either allele > 1
    DIRECTION: foreach my $direction(@{_get_trim_directions($orig_b_ref, $orig_b_alt)}) {

      # first make copies of everything
      my $pos = $orig_b_pos;
      my $ref = $orig_b_ref;
      my $alt = $orig_b_alt;

      ($ref, $alt, $pos) = @{trim_sequences($ref, $alt, $pos, undef, 1, $direction)};

      if(my $orig_a_alt = $minimised_a_alleles{"$ref\_$alt\_$pos"}) {
        push @matches, {
          a_allele => $orig_a_alt,
          a_index  => $a_indexes{$orig_a_alt},
          b_allele => $orig_b_alt,
          b_index  => $i
        };
        last DIRECTION;
      }
    }

    $i++;
  }

  return \@matches;
}

# accessory method for get_matched_variant_alleles()
sub _get_trim_directions {
  my ($ref, $alt) = @_;
  if(length($ref) > 1 || length($alt) > 1) {
    return [0, 1];
  }
  else {
    return [0];
  }
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


=head2 raw_freqs_from_gts

  Arg [1]    : arrayref of Bio::EnsEMBL::Variation::SampleGenotype $gts
  Arg [2]    : arrayref of Bio::EnsEMBL::Variation::Population $pops
  Arg [3]    : hashref of sample IDs per population ID $pop_hash, as generated by
               Bio::EnsEMBL::Variation::Population::DBSQL::PopulationAdaptor::_get_sample_population_hash
               {
                 pop_id1 => {
                   sample_id1 => 1,
                   sample_id2 => 1
                 },
                 pop_id2 => {
                   sample_id3 => 1,
                   sample_id4 => 1
                 }
               }
  Arg [4]    : (optional) hashref of missing allele counts by sample ID
               e.g. for alleles found in VCF but not in VariationFeature
  Example    : $raw_obj_hash = raw_freqs_from_gts($gts, $pops, $pop_hash, $missing)
  Description: Computes proto-objects representing Alleles and PopulationGenotypes
               with aggregated counts and frequencies from a list of genotype objects,
               limited to the given populations, with population/sample mapping
               described in $pop_hash
  Returntype : hashref with two keys:
               {
                 Allele => $arrayref_of_alleles,
                 PopulationGenotype => $arrayref_of_populationgenotypes
               }
  Exceptions : none
  Caller     : Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor
  Status     : stable

=cut


sub raw_freqs_from_gts {
  my ($genotypes, $pops, $pop_sample_hash, $missing) = @_;

  # allow for missing being undef
  $missing ||= {};

  # do initial grouping by subsnp
  my %by_ss = ();
  push @{$by_ss{$_->subsnp || ''}}, $_ for @$genotypes;
  
  my (@alleles, @pop_genotypes, $missing_by_pop);
    
  # go through given populations
  foreach my $pop(@$pops) {
    my $pop_id = $pop->dbID;
    
    foreach my $ss(keys %by_ss) {

      my (%a_counts, $a_total, %g_counts, $g_total);

      # get all genotypes for samples in this population
      foreach my $genotype(grep {$pop_sample_hash->{$pop_id}->{$_->{_sample_id} ||= $_->sample->dbID}} @{$by_ss{$ss}}) {

        # count genotypes
        $g_total++;
        $g_counts{$genotype->genotype_string(1)}++;

        # count alleles
        foreach my $allele(@{$genotype->genotype}) {
          $a_total++;
          $a_counts{$allele}++;
        }
      }
      
      # add missing allele data
      foreach my $sample_id(keys %{$pop_sample_hash->{$pop_id}}) {
        if(my $missing_by_sample = $missing->{$sample_id}) {
          foreach my $allele(keys %$missing_by_sample) {
            $missing_by_pop->{$pop_id}->{$allele} = 1;
            $a_total += $missing_by_sample->{$allele};
          }
        }
      }

      # add missing pop data
      $g_total++ for grep {$missing->{$_}} keys %{$pop_sample_hash->{$pop_id}};
      
      next unless (%a_counts or %g_counts) && ($a_total or $g_total);

      # calculate freqs and create pseudo-objects
      push @alleles, map {
        {
          allele     => $_,
          count      => $a_counts{$_} || 0,
          frequency  => $a_counts{$_} ? $a_counts{$_} / $a_total : 0,
          population => $pop,
          subsnp     => $ss eq '' ? undef : $ss,
          _missing_alleles => $missing_by_pop->{$pop_id} || {},
        }
      } keys %a_counts;

      push @pop_genotypes, map {
        {
          genotype   => [split /\|/, $_],
          count      => $g_counts{$_} || 0,
          frequency  => $g_counts{$_} ? $g_counts{$_} / $g_total : 0,
          population => $pop,
          subsnp     => $ss eq '' ? undef : $ss
        }
      } keys %g_counts;
    }
  }

  return {
    Allele             => \@alleles,
    PopulationGenotype => \@pop_genotypes,
  };
}


1;

#
# Ensembl module for Bio::EnsEMBL::Variation::TranscriptVariation
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::TranscriptVariation

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::TranscriptVariation;

  $tr_var = Bio::EnsEMBL::Variation::TranscriptVariation->new
    (-variation_feature => $feature,
     -transcript        => $transcript,
     -pep_allele_string => 'N/K',
     -cdna_start        => 1127,
     -cdna_end          => 1127,
     -translation_start => 318,
     -translation_end   => 318,
     -type              => 'NON_SYNONYMOUS_CODING');


  print "variation: ", $tr_var->variation_feature()->variation_name(), "\n";
  print "transcript: ", $tr_var->transcript->stable_id(), "\n";
  print "type: ", $tr_var->type(), "\n";
  print "cdna coords: ", $tr_var->cdna_start(), '-', $tr_var->cdna_end(), "\n";
  print "pep coords: ", $tr_var->translation_start(), '-',
        $tr_var->translation_end(), "\n";
  print "amino acid change: ", $tr_var->pep_allele_string(), "\n";


=head1 DESCRIPTION

A TranscriptVariation object represents a variation feature which is in close
proximity to an Ensembl transcript.  A TranscriptVariation object has several
attributes which define the relationship of the variation to the transcript.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::TranscriptVariation;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Storable;

our @ISA = ('Bio::EnsEMBL::Storable');

our %VALID_TYPES = ('INTRONIC' => 1,
                    'UPSTREAM' => 1,
                    'DOWNSTREAM' => 1,
                    'SYNONYMOUS_CODING' => 1,
                    'NON_SYNONYMOUS_CODING', => 1,
                    'FRAMESHIFT_CODING' => 1,
                    '5PRIME_UTR' => 1,
                    '3PRIME_UTR' => 1);

=head2 new

  Arg [-VARIATION_FEATURE] :
    Bio::EnsEMBL::Variation::VariationFeature - The variation feature that was
    used to derive this TranscriptVariation.

  Arg [-TRANSCRIPT] :
    Bio::EnsEMBL::Transcript - The transcript affected by this
    TranscriptVariation

  Arg [-PEP_ALLELE_STRING] :
    string - A '/' delimited string representing amino acids altered by this
    variation

  Arg [-CDNA_START] :
    The start of this variation on the associated transcript in cdna
    coordinates

  Arg [-CDNA_END] :
    The end of this variation on the associated transcript in cdna coordinates

  Arg [-TRANSLATION_START] :
    The start of this variation on the translation of the associated transcript
    in peptide coordinates

  Arg [-TRANSLATION_END] :
    The end of this variation on the translation of the associated transcript
    in peptide coordinates

  Arg [-TYPE] :
    The type of this TranscriptVariation.  Must be one of:
    'INTRONIC', 'UPSTREAM', 'DOWNSTREAM', 'SYNONYMOUS_CODING',
    'NON_SYNONYMOUS_CODING', 'FRAMESHIFT_CODING', '5PRIME_UTR', '3PRIME_UTR'

  Example    : 
    $tr_var = Bio::EnsEMBL::Variation::TranscriptVariation->new
      (-variation_feature => $feature,
       -transcript        => $transcript,
       -pep_allele_string => 'N/K',
       -cdna_start        => 1127,
       -cdna_end          => 1127,
       -translation_start => 318,
       -translation_end   => 318,
       -type              => 'NON_SYNONYMOUS_CODING');

  Description: Constructor. Instantiates a
               Bio::EnsEMBL::Variation::TranscriptVariation object
  Returntype : Bio::EnsEMBL::Variation::TranscriptVariation
  Exceptions : throw on bad argument
  Caller     : general, TranscriptVariationAdaptor

=cut

sub new {
  my $class = shift;

  my ($vf, $tr, $pep_allele, $cdna_start,$cdna_end, $tl_start,$tl_end, $type) =
    rearrange([qw(VARIATION_FEATURE TRANSCRIPT PEP_ALLELE_STRING CDNA_START
                  CDNA_END TRANSLATION_START TRANSLATION_END TYPE)], @_);

  if(defined($vf) &&
     (!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature'))) {
    throw('VariationFeature argument expected');
  }

  if(defined($tr) && (!ref($tr) || !$tr->isa('Bio::EnsEMBL::Transcript'))) {
    throw('Transcript argument expected');
  }

  if(defined($type)) {
    $type = uc($type);
    if(!$VALID_TYPES{$type}) {
      my $valid = join(',',map({"'$_'"} keys(%VALID_TYPES)));
      throw("Type argument must be one of: $valid");
    }
  }

  if(defined($cdna_start) && ($cdna_start !~ /^\d+$/ || $cdna_start < 1)) {
    throw('CDNA start must be greater than or equal to 1');
  }

  if(defined($cdna_end) && ($cdna_end !~ /^\d+$/ || $cdna_start < 0)) {
    throw('CDNA end must be greater than or equal to 0');
  }

  if(defined($tl_start) && ($tl_start !~ /^\d+$/ || $tl_start < 1)) {
    throw('Translation start must be greater than or equal to 1');
  }

  if(defined($tl_end) && ($tl_end !~ /^\d+$/ || $tl_start < 0)) {
    throw('Translation end must be greater than or equal to 0');
  }

  return bless {'variation_feature' => $vf,
                'transcript'        => $tr,
                'pep_allele_string' => $pep_allele,
                'cdna_start'        => $cdna_start,
                'cdna_end'          => $cdna_end,
                'translation_start' => $tl_start,
                'translation_end'   => $tl_end,
                'type'              => $type}, $class;
}



=head2 transcript

  Arg [1]    : (optional) Bio::EnsEMBL::Transcript $transcript
  Example    : print $trvar->transcript()->stable_id(), "\n";
  Description: Getter/Setter for the Transcript that is affected by this
               TranscriptVariation.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub transcript {
  my $self = shift;

  if(@_) {
    my $tr = shift;
    if(defined($tr) && (!ref($tr) || !$tr->isa('Bio::EnsEMBL::Transcript'))) {
      throw('Bio::EnsEMBL::Transcript argument expected');
    }
    $self->{'transcript'} = $tr;
  }

  return $self->{'transcript'};
}




=head2 variation_feature

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    : print $trvar->variation_feature()->variation_name(), "\n";
  Description: Getter/Setter for the VariationFeature associated with this
               transcript variation.
  Returntype : Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub variation_feature {
  my $self = shift;

  if(@_) {
    my $vf = shift;
    if(defined($vf) && 
       (!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature'))) {
      throw('Bio::EnsEMBL::Variation::VariationFeature argument expected');
    }
    $self->{'variation_feature'} = $vf;
  }

  return $self->{'variation_feature'};
}




=head2 pep_allele_string

  Arg [1]    : string $allele (optional) 
               The new value to set the pep_allele_string attribute to
  Example    : $pep_allele_string = $obj->pep_allele_string()
  Description: Getter/Setter for the pep_allele_string attribute.  The
               pep allele string is a '/' delimited string of amino acid
               codes representing the change to the peptide made by this
               variation.  A '-' represents one half of a insertion/deletion.
               The reference allele (the one found in the ensembl peptide)
               should be first.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub pep_allele_string{
  my $self = shift;
  return $self->{'pep_allele_string'} = shift if(@_);
  return $self->{'pep_allele_string'};
}



=head2 cdna_start

  Arg [1]    : (optional) int $start
  Example    : $cdna_start = $trvar->cdna_start();
  Description: Getter/Setter for the start position of this variation on the
               transcript in CDNA coordinates.
  Returntype : int
  Exceptions : throw if $start is not an int
               throw if $start < 1
  Caller     : general

=cut

sub cdna_start {
  my $self = shift;

  if(@_) {
    my $cdna_start = shift;
    if(defined($cdna_start) && ($cdna_start !~ /^\d+$/ || $cdna_start < 1)) {
      throw('cdna start must be an integer greater than 0');
    }
    $self->{'cdna_start'} = $cdna_start;
  }
  return $self->{'cdna_start'};
}



=head2 cdna_end

  Arg [1]    : (optional) int $end
  Example    : $cdna_end = $trvar->cdna_end();
  Description: Getter/Setter for the end position of this variation on the
               transcript in cdna coordinates.
  Returntype : int
  Exceptions : throw if $end is not an int
               throw if $end < 0
  Caller     : general

=cut

sub cdna_end {
  my $self = shift;

  if(@_) {
    my $cdna_end = shift;
    if(defined($cdna_end) && ($cdna_end !~ /^\d+$/ || $cdna_end < 0)) {
      throw('cdna end must be an integer greater than or equal to 0');
    }
    $self->{'cdna_end'} = $cdna_end;
  }

  return $self->{'cdna_end'};
}



=head2 translation_start

  Arg [1]    : (optional) int $tl_start
  Example    : $tl_start = $trvar->translation_start();
  Description: Getter/Setter for the start position of this variation on the
               translation of the associated transcript in peptide coordinates.
  Returntype : int
  Exceptions : throw if $start is not an int
               throw if $start < 0
  Caller     : general

=cut

sub translation_start {
  my $self = shift;

  if(@_) {
    my $tl_start = shift;
    if(defined($tl_start) && ($tl_start !~ /^\d+$/ || $tl_start < 1)) {
      throw('translation start must be an integer greater than or equal to 1');
    }
    $self->{'translation_start'} = $tl_start;
  }

  return $self->{'translation_start'};
}


=head2 translation_end

  Arg [1]    : (optional) int $tl_end
  Example    : $tl_end = $trvar->translation_end();
  Description: Getter/Setter for the end position of this variation on the
               translation of the associated transcript in peptide coordinates.
  Returntype : int
  Exceptions : throw if $end is not an int
               throw if $end < 0
  Caller     : general

=cut

sub translation_end {
  my $self = shift;

  if(@_) {
    my $tl_end = shift;
    if(defined($tl_end) && ($tl_end !~ /^\d+$/ || $tl_end < 0)) {
      throw('translation end must be an integer greater than or equal to 0');
    }
    $self->{'translation_end'} = $tl_end;
  }

  return $self->{'translation_end'};
}


=head2 type

  Arg [1]    : (optional) string $type
  Example    : if($tr_var->type() eq 'INTRONIC') { do_something(); }
  Description: Getter/Setter for the type of this transcript variation.
               Allowed values are:  'INTRONIC','UPSTREAM','DOWNSTREAM',
               'SYNONYMOUS_CODING','NON_SYNONYMOUS_CODING','FRAMESHIFT_CODING',
               '5PRIME_UTR','3PRIME_UTR'
  Returntype : string
  Exceptions : throw if provided argument is not one of the valid strings
  Caller     : general

=cut

sub type {
  my $self = shift;

  if(@_) {
    my $type = shift;

    if(defined($type)) {
      $type = uc($type);
      if(!$VALID_TYPES{$type}) {
        my $valid = join(',',map({"'$_'"} keys(%VALID_TYPES)));
        throw("Type argument must be one of: $valid");
      }
    }
    $self->{'type'} = $type;
  }

  return $self->{'type'};
}


1;

# EnsEMBL module for ConsequenceType
# Copyright EMBL-EBI/Sanger center 2005
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Variation::ConsequenceType

=head1 SYNOPSIS


=head1 DESCRIPTION

Represents the effect of a Variation in a Transcript

=cut

package Bio::EnsEMBL::Variation::ConsequenceType;

use strict;
use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&type_variation);

use Bio::EnsEMBL::Utils::Exception qw(warning);
use Bio::EnsEMBL::Variation::VariationFeature; #to get the consequence_types priorities

=head2 new

  Arg [1]    : (optional) int $transcript_id
  Arg [2]    : (optional) int $variation_feature_id
  Arg [2]    : (optional) int $start
  Arg [3]    : (optional) int $end
  Arg [4]    : (optional) int $strand
  Arg [5]    : (optional) refarray $alleles
  Example    : $synonym = Bio::EnsEMBL::Variation::ConsequenceType->new($transcript_id,$variation_feature_id,$start,$end,$strand,['A','C']);
  Description: Creates a new ConsequenceType
  Returntype : Bio::EnsEMBL::Variation::ConsequenceType
  Exceptions : none
  Caller     : general

=cut

sub new {
  my ($caller, $transcript_id, $variation_feature_id, $start, $end, $strand, $alleles) = @_;

  my $class = ref($caller) || $caller;

  return bless( {'transcript_id'   => $transcript_id,
		 'variation_feature_id' => $variation_feature_id,
		 'alleles'  => $alleles,
		 'start' => $start,
		 'end'   => $end,
		 'strand' => $strand}, $class );
}


=head2 transcript_id

  Arg [1]    : (optional) int $transcript_id
  Example    : $transcript_id = $consequence_type->transcript_id;
  Description: Getter/Setter for the internal id of the transcript_id calculated
               the effect of the Variation
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub transcript_id {
  my $self = shift;

  if(@_) {
    $self->{'transcript_id'} = shift;
  }

  return $self->{'transcript_id'};
}


=head2 variation_feature_id

  Arg [1]    : (optional) int $variation_feature_id
  Example    : $variation_feature_id = $consequence_type->variation_feature_id;
  Description: Getter/Setter for the variation_feature affecting the transcript
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub variation_feature_id {
  my $self = shift;

  if(@_) {
    $self->{'variation_feature_id'} = shift;
  }

  return $self->{'variation_feature_id'};
}

=head2 alleles

  Arg [1]    : (optional) array ref $alleles
  Example    : @alleles = @{$consequence_type->alleles};
  Description: Getter/Setter for the alleles for the variation
  Returntype : reference to array
  Exceptions : none
  Caller     : general

=cut

sub alleles {
  my $self = shift;

  if(@_) {
    $self->{'alleles'} = shift;
  }

  return $self->{'alleles'};
}

=head2 start

  Arg [1]    : (optional) int $start
  Example    : $start = $consequence_type->start
  Description: Getter/Setter for the start of the variation in the sequence
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub start {
  my $self = shift;

  if(@_) {
    $self->{'start'} = shift;
  }
  
  return $self->{'start'}
}


=head2 end

  Arg [1]    : (optional) int $end
  Example    : $end = $consequence_type->end
  Description: Getter/Setter for the end of the variation in the sequence
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub end {
  my $self = shift;

  if(@_) {
    $self->{'end'} = shift;
  }
  
  return $self->{'end'}
}


=head2 strand

  Arg [1]    : (optional) int $strand
  Example    : $strand = $consequence_type->strand
  Description: Getter/Setter for the strand of the variation in the sequence
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub strand {
  my $self = shift;

  if(@_) {
    $self->{'strand'} = shift;
  }
  
  return $self->{'strand'}
}


=head2 aa_start

  Arg [1]    : (optional) int $aa_start
  Example    : $aa_start = $consequence_type->aa_start
  Description: Getter/Setter for the start of the aa in peptide coordinates
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub aa_start {
  my $self = shift;

  if(@_) {
    $self->{'aa_start'} = shift;
  }
  
  return $self->{'aa_start'}
}

=head2 aa_end

  Arg [1]    : (optional) int $aa_end
  Example    : $aa_end = $consequence_type->aa_end
  Description: Getter/Setter for the end of the aa in peptide coordinates
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub aa_end {
  my $self = shift;

  if(@_) {
    $self->{'aa_end'} = shift;
  }
  
  return $self->{'aa_end'}
}

=head2 cdna_start

  Arg [1]    : (optional) int $cdna_start
  Example    : $cdna_start = $consequence_type->cdna_start
  Description: Getter/Setter for the start of the variation in the cdna coordinates
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub cdna_start {
  my $self = shift;

  if(@_) {
    $self->{'cdna_start'} = shift;
  }
  
  return $self->{'cdna_start'}
}

=head2 cdna_end

  Arg [1]    : (optional) int $cdna_end
  Example    : $cdna_end = $consequence_type->cdna_end
  Description: Getter/Setter for the end of the variation in cdna coordinates
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub cdna_end {
  my $self = shift;

  if(@_) {
    $self->{'cdna_end'} = shift;
  }
  
  return $self->{'cdna_end'}
}


=head2 type

  Arg [1]    : (optional) string $type 
               (possible values "INTRONIC", "UPSTREAM", "DOWNSTREAM", "SYNONYMOUS_CODING",
                 "NON_SYNONYMOUS_CODING", "FRAMESHIFT_CODING", "5PRIME_UTR", "3PRIME_UTR",
		 "STOP_GAINED", "STOP_LOST", "UTR", "INTERGENIC")
  Example    : $consequence_type = $consequence_type->type
  Description: Getter/Setter for the consequence type of the variation in the transcript
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub type {
  my $self = shift;

  my %consequence_types = %Bio::EnsEMBL::Variation::VariationFeature::CONSEQUENCE_TYPES;
  if(@_) {
      my $type = shift;
      if (defined $consequence_types{$type}){
	  $self->{'type'} = $type;
      }
      else{
	  warning("Trying to set the consequence type to a not valid value. Possible values: ",keys %consequence_types,"\n");
      }
  }
  return $self->{'type'}
}

=head2 aa_alleles

  Arg [1]    : (optional) string $aa_alleles
  Example    : $aa_alleles = $consequence_type->aa_alleles
  Description: Getter/Setter for the aa that changes in the transcript
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub aa_alleles {
  my $self = shift;

  if(@_) {
    $self->{'aa_alleles'} = shift;
  }
  
  return $self->{'aa_alleles'}
}

1;


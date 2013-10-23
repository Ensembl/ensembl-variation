=head1 LICENSE

Copyright 2013 Ensembl

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
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

# EnsEMBL module for ConsequenceType
#
#
# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Variation::ConsequenceType

=head1 SYNOPSIS


=head1 DESCRIPTION

Represents the effect of a Variation in a Transcript

=cut

package Bio::EnsEMBL::Variation::ConsequenceType;

use strict;

use Bio::EnsEMBL::Utils::Exception qw(warning);

#conversion of consequence type to bit value

#there is a special type, SARA, that only applies to the effect of the Alleles, not Variations, and is equivalent
#to Same As Reference Allele, meaning that the Allele is the same as in reference sequence, so has no effect 
#but it is not stored anywhere in the database and need no conversion at all
#when creating the VariationFeature object, thus the absence in the hash
our %CONSEQUENCE_TYPES = (
  'ESSENTIAL_SPLICE_SITE' => 1,
  'STOP_GAINED' => 2,
  'STOP_LOST' => 4,
  'COMPLEX_INDEL' => 8,
  'FRAMESHIFT_CODING' => 16,
  'NON_SYNONYMOUS_CODING' => 32,
  'SPLICE_SITE' => 64,
  'PARTIAL_CODON' => 128,
  'SYNONYMOUS_CODING' => 256,
  'REGULATORY_REGION' => 512,
  'WITHIN_MATURE_miRNA' => 1024,
  '5PRIME_UTR' => 2048,
  '3PRIME_UTR' => 2094,
  'UTR'        => 4096,
  'INTRONIC' => 8192,
  'NMD_TRANSCRIPT' => 16384,
  'WITHIN_NON_CODING_GENE' => 32768,
  'UPSTREAM' => 65536,
  'DOWNSTREAM' => 131072,
  'HGMD_MUTATION' => 262144,
  'NO_CONSEQUENCE' => 524288,
  'INTERGENIC' => 1048576,
  '_'          => 2097152,
);


our %CONSEQUENCE_DESCRIPTIONS = (
  'ESSENTIAL_SPLICE_SITE'  => 'In the first 2 or the last 2 basepairs of an intron',
  'STOP_GAINED'            => 'In coding sequence, resulting in the gain of a stop codon',
  'STOP_LOST'              => 'In coding sequence, resulting in the loss of a stop codon',
  'COMPLEX_INDEL'          => 'Insertion or deletion that spans an exon/intron or coding sequence/UTR border',
  'FRAMESHIFT_CODING'      => 'In coding sequence, resulting in a frameshift',
  'NON_SYNONYMOUS_CODING'  => 'In coding sequence and results in an amino acid change in the encoded peptide sequence',
  'SPLICE_SITE'            => '1-3 bps into an exon or 3-8 bps into an intron',
  'PARTIAL_CODON'          => 'Located within the final, incomplete codon of a transcript whose end coordinate is unknown',
  'SYNONYMOUS_CODING'      => 'In coding sequence, not resulting in an amino acid change (silent mutation)',
  'REGULATORY_REGION'      => 'In regulatory region annotated by Ensembl',
  'WITHIN_MATURE_miRNA'    => 'Located within a microRNA',
  '5PRIME_UTR'             => 'In 5 prime untranslated region',
  '3PRIME_UTR'             => 'In 3 prime untranslated region',
  'INTRONIC'               => 'In intron',
  'NMD_TRANSCRIPT'         => 'Located within a transcript predicted to undergo nonsense-mediated decay',
  'WITHIN_NON_CODING_GENE' => 'Located within a gene that does not code for a protein',
  'UPSTREAM'               => 'Within 5 kb upstream of the 5 prime end of a transcript',
  'DOWNSTREAM'             => 'Within 5 kb downstream of the 3 prime end of a transcript',
  'HGMD_MUTATION'          => 'Mutation from the HGMD database - consequence unknown',
  'INTERGENIC'             => 'More than 5 kb either upstream or downstream of a transcript',
);


our %CONSEQUENCE_LABELS = (
  'ESSENTIAL_SPLICE_SITE'  => 'Essential splice site',
  'STOP_GAINED'            => 'Stop gained',
  'STOP_LOST'              => 'Stop lost',
  'COMPLEX_INDEL'          => 'Complex in/del',
  'FRAMESHIFT_CODING'      => 'Frameshift coding',
  'NON_SYNONYMOUS_CODING'  => 'Non-synonymous coding',
  'SPLICE_SITE'            => 'Splice site',
  'PARTIAL_CODON'          => 'Partial codon',
  'SYNONYMOUS_CODING'      => 'Synonymous coding',
  'REGULATORY_REGION'      => 'Regulatory region',
  'WITHIN_MATURE_miRNA'    => 'Within mature miRNA',
  '5PRIME_UTR'             => '5 prime UTR',
  '3PRIME_UTR'             => '3 prime UTR',
  'INTRONIC'               => 'Intronic',
  'NMD_TRANSCRIPT'         => 'NMD transcript',
  'WITHIN_NON_CODING_GENE' => 'Within non-coding gene',
  'UPSTREAM'               => 'Upstream',
  'DOWNSTREAM'             => 'Downstream',
  'HGMD_MUTATION'          => 'HGMD mutation',
  'INTERGENIC'             => 'Intergenic',
);

# hash storing whether consequence affects peptide sequence
our %AFFECTS_PEPTIDE = (
  'STOP_GAINED'				=> 1,
  'STOP_LOST'				=> 1,
  'COMPLEX_INDEL'			=> 1,
  'FRAMESHIFT_CODING'		=> 1,
  'NON_SYNONYMOUS_CODING'	=> 1,
  'PARTIAL_CODON'			=> 1,
);


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
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

=cut

sub cdna_end {
  my $self = shift;

  if(@_) {
    $self->{'cdna_end'} = shift;
  }
  
  return $self->{'cdna_end'}
}


=head2 type

  Arg [1]    : string $type 
               (possible types 'FRAMESHIFT_CODING','STOP_GAINED','STOP_LOST','NON_SYNONYMOUS_CODING',
		'SYNONYMOUS_CODING','REGULATORY_REGION','WITHIN_MATURE_miRNA','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM','WITHIN_NON_CODING_GENE','INTERGENIC', 'SARA')
  Example    : $consequence_type = $consequence_type->type
  Description: Getter/Setter for consequence type of the variation in the transcript
  Returntype : none
  Exceptions : warning if the consequence type is not recognised
  Caller     : general
  Status     : At Risk

=cut

sub type {
  my $self = shift;

  if(@_) {
      my $type = shift;
      #there is a special type, SARA, that only applies to the effect of the Alleles, and is equivalent
      #to Same As Reference Allele, which is not stored anywhere in the database and need no conversion at all
      #when creating the VariationFeature object, thus the absence in the hash
      if (defined $CONSEQUENCE_TYPES{$type} || $type eq 'SARA'){
	  push @{$self->{'type'}}, $type;
      }
      else{
	  warning("Trying to set the consequence type to a not valid value. Possible values: ",keys %CONSEQUENCE_TYPES,"\n");
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
  Status     : At Risk

=cut

sub aa_alleles {
  my $self = shift;

  if(@_) {
    $self->{'aa_alleles'} = shift;
  }
  
  return $self->{'aa_alleles'}
}


=head2 codon

  Arg [1]    : (optional) string $codon
  Example    : $codon = $consequence_type->codon
  Description: Getter/Setter for the codon affected by that Allele in the transcript
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub codon {
  my $self = shift;

  if(@_) {
    $self->{'codon'} = shift;
  }
  
  return $self->{'codon'}
}


=head2 codons

  Arg [1]    : (optional) string $codons
  Example    : $codons = $consequence_type->codons
  Description: Getter/Setter for the possible codons
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub codons {
  my $self = shift;

  if(@_) {
    $self->{'codons'} = shift;
  }
  
  return $self->{'codons'}
}

=head2 cds_start

  Arg [1]    : (optional) int $cds_start
  Example    : $cds_start = $consequence_type->cds_start
  Description: Getter/Setter for the start of the variation in the coding sequence
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub cds_start {
  my $self = shift;

  if(@_) {
    $self->{'cds_start'} = shift;
  }
  
  return $self->{'cds_start'}
}

=head2 cds_end

  Arg [1]    : (optional) int $cds_end
  Example    : $cds_end = $consequence_type->cds_end
  Description: Getter/Setter for the end of the variation in the coding sequence
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub cds_end {
  my $self = shift;

  if(@_) {
    $self->{'cds_end'} = shift;
  }
  
  return $self->{'cds_end'}
}

=head2 display_consequence

  Arg [1]    : (optional) string $consequence_type
  Example    : $display_consequence = $ct->display_consequence();
  Description: Getter for the consequence type to display,
               when more than one
  Returntype : string
  Exceptions : throw on incorrect argument
  Caller     : webteam
  Status     : At Risk

=cut

sub display_consequence{
    my $self = shift;
 
    my $highest_priority;
    #get the value to display from the consequence_type attribute
    $highest_priority = 'INTERGENIC';
    foreach my $ct (@{$self->type}){
	if ($CONSEQUENCE_TYPES{$ct} < $CONSEQUENCE_TYPES{$highest_priority}){
	    $highest_priority = $ct;
	}
    }

    return $highest_priority;
}

sub empty_type{
    my $self = shift;

    $self->{'type'} = ();
    return $self->type;
}

1;


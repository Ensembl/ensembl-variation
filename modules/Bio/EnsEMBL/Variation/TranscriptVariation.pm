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

=head1 NAME

Bio::EnsEMBL::Variation::TranscriptVariation

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::TranscriptVariation;

    my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
        -transcript        => $transcript,
        -variation_feature => $var_feat
    );

    print "consequence type: ", (join ",", @{$tv->consequence_type}), "\n";
    print "cdna coords: ", $tv->cdna_start, '-', $tv->cdna_end, "\n";
    print "cds coords: ", $tv->cds_start, '-', $tv->cds_end, "\n";
    print "pep coords: ", $tv->translation_start, '-',$tv->translation_end, "\n";
    print "amino acid change: ", $tv->pep_allele_string, "\n";
    print "codon change: ", $tv->codons, "\n";
    print "allele sequences: ", (join ",", map { $_->variation_feature_seq } 
        @{ $tv->get_all_TranscriptVariationAlleles }), "\n";

=head1 DESCRIPTION

A TranscriptVariation object represents a variation feature which is in close
proximity to an Ensembl transcript. A TranscriptVariation object has several
attributes which define the relationship of the variation to the transcript.

=cut

package Bio::EnsEMBL::Variation::TranscriptVariation;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref check_ref);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap within_cds);
use Bio::EnsEMBL::Variation::BaseTranscriptVariation;

use base qw(Bio::EnsEMBL::Variation::BaseTranscriptVariation Bio::EnsEMBL::Variation::VariationFeatureOverlap);

=head2 new

  Arg [-TRANSCRIPT] : 
    The Bio::EnsEMBL::Transcript associated with the given VariationFeature

  Arg [-VARIATION_FEATURE] :
    The Bio::EnsEMBL::VariationFeature associated with the given Transcript

  Arg [-ADAPTOR] :
    A Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor

  Arg [-DISAMBIGUATE_SINGLE_NUCLEOTIDE_ALLELES] :
    A flag indiciating if ambiguous single nucleotide alleles should be disambiguated
    when constructing the TranscriptVariationAllele objects, e.g. a Variationfeature
    with an allele string like 'T/M' would be treated as if it were 'T/A/C'. We limit
    ourselves to single nucleotide alleles to avoid the combinatorial explosion if we
    allowed longer alleles with potentially many ambiguous bases.

  Example : 
    my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
        -transcript        => $transcript,
        -variation_feature => $var_feat
    );

  Description: Constructs a new TranscriptVariation instance given a VariationFeature
               and a Transcript, most of the work is done in the VariationFeatureOverlap
               superclass - see there for more documentation.
  Returntype : A new Bio::EnsEMBL::Variation::TranscriptVariation instance 
  Exceptions : throws unless both VARIATION_FEATURE and TRANSCRIPT are supplied
  Status     : Stable

=cut 

sub new {
    my $class = shift;

    my %args = @_;

    # swap a '-transcript' argument for a '-feature' one for the superclass
    unless($args{'-feature'} ||= delete $args{'-transcript'}) {
      for my $arg (keys %args) {
        if (lc($arg) eq '-transcript') {
          $args{'-feature'} = delete $args{$arg};
        }
      }
    }

    # call the superclass constructor
    my $self = $class->SUPER::new(%args) || return undef;

    # rebless the alleles from vfoas to tvas
    map { bless $_, 'Bio::EnsEMBL::Variation::TranscriptVariationAllele' } 
        @{ $self->get_all_BaseVariationFeatureOverlapAlleles };
    
    return $self;
}

sub get_TranscriptVariationAllele_for_allele_seq {
    my ($self, $allele_seq) = @_;
    return $self->SUPER::get_VariationFeatureOverlapAllele_for_allele_seq($allele_seq);
}

=head2 add_TranscriptVariationAllele

  Arg [1]    : A Bio::EnsEMBL::Variation::TranscriptVariationAllele instance
  Description: Add an allele to this TranscriptVariation
  Returntype : none
  Exceptions : throws if the argument is not the expected type
  Status     : Stable

=cut

sub add_TranscriptVariationAllele {
    my ($self, $tva) = @_;
    assert_ref($tva, 'Bio::EnsEMBL::Variation::TranscriptVariationAllele') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;
    return $self->SUPER::add_VariationFeatureOverlapAllele($tva);
}

=head2 get_reference_TranscriptVariationAllele

  Description: Get the object representing the reference allele of this TranscriptVariation
  Returntype : Bio::EnsEMBL::Variation::TranscriptVariationAllele instance
  Exceptions : none
  Status     : Stable

=cut

sub get_reference_TranscriptVariationAllele {
    my $self = shift;
    return $self->SUPER::get_reference_VariationFeatureOverlapAllele(@_);
}

=head2 get_all_alternate_TranscriptVariationAlleles

  Description: Get a list of the alternate alleles of this TranscriptVariation
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariationAllele objects
  Exceptions : none
  Status     : Stable

=cut

sub get_all_alternate_TranscriptVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_alternate_VariationFeatureOverlapAlleles(@_);
}

=head2 get_all_TranscriptVariationAlleles

  Description: Get a list of the all the alleles, both reference and alternate, of 
               this TranscriptVariation
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariationAllele objects
  Exceptions : none
  Status     : Stable

=cut

sub get_all_TranscriptVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_VariationFeatureOverlapAlleles(@_);
}


=head2 cdna_allele_string

  Description: Return a '/' delimited string of the alleles of this variation with respect
               to the associated transcript
  Returntype : string
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cdna_allele_string {
    my $self = shift;
    
    unless ($self->{_cdna_allele_string}) {
        $self->{_cdna_allele_string} = join '/', map { $_->feature_seq } @{ $self->get_all_TranscriptVariationAlleles };
    }
    
    return $self->{_cdna_allele_string};
}

=head2 pep_allele_string

  Description: Return a '/' delimited string of amino acid codes representing
               all the possible changes made to the peptide by this variation
  Returntype : string
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub pep_allele_string {
    my $self = shift;
    
    unless ($self->{_pep_allele_string}) {
        
        my @peptides = grep { defined } map { $_->peptide } @{ $self->get_all_TranscriptVariationAlleles };
        
        $self->{_pep_allele_string} = join '/', @peptides;
    }
    
    return $self->{_pep_allele_string};
}

=head2 codons

  Description: Return a '/' delimited string of all possible codon sequences.  
               The variant sequence within the codon will be capitalised while 
               the rest of the codon sequence will be in lower case
  Returntype : string
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub codons {
    my $self = shift;
    
    unless ($self->{_display_codon_allele_string}) {
  
        my @display_codons = grep { defined } map { $_->display_codon } @{ $self->get_all_TranscriptVariationAlleles }; 

        $self->{_display_codon_allele_string} = join '/', @display_codons;
    }

    return $self->{_display_codon_allele_string};
}

=head2 codon_position

  Description: For variations that fall in the CDS, returns the base in the
               codon that this variation falls in
  Returntype : int - 1, 2 or 3
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub codon_position {
    my $self = shift;
    
    unless ($self->{_codon_position}) {

        my $cdna_start = $self->cdna_start;
        
        my $tran_cdna_start = $self->transcript->cdna_coding_start;
       
        # we need to take into account the exon phase
        my $exon_phase = $self->transcript->start_Exon->phase;

        my $phase_offset = $exon_phase > 0 ? $exon_phase : 0;

        if (defined $cdna_start && defined $tran_cdna_start) {
            $self->{_codon_position} = (($cdna_start - $tran_cdna_start + $phase_offset) % 3) + 1;
        }
    }
    
    return $self->{_codon_position};
}

=head2 affects_peptide

  Description: Check if any of this TranscriptVariation's alleles change the
               resultant peptide sequence
  Returntype : boolean
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub affects_peptide {
    my $self = shift;
    return scalar grep { $_->SO_term =~ /stop|missense|frameshift|inframe|initiator/ } map {@{$_->get_all_OverlapConsequences}} @{ $self->get_all_alternate_TranscriptVariationAlleles }; 
}


sub _protein_function_predictions {
    
    my ($self, $analysis) = @_;

    my $tran = $self->transcript;

    my $matrix = $tran->{_variation_effect_feature_cache}->{protein_function_predictions}->{$analysis};

    unless ($matrix || exists($tran->{_variation_effect_feature_cache}->{protein_function_predictions}->{$analysis})) {
        return undef unless $self->{adaptor} && $self->{adaptor}->db;
        
        my $pfpma = $self->{adaptor}->db->get_ProteinFunctionPredictionMatrixAdaptor;
           
        $matrix = $pfpma->fetch_by_analysis_translation_md5($analysis, $self->_translation_md5);

        $tran->{_variation_effect_feature_cache}->{protein_function_predictions}->{$analysis} = $matrix;
    }

    return $matrix; 
}

=head2 hgvs_genomic

  Description: Return the strings representing the genomic-level effect of each of the alleles 
               of this variation in HGVS format
  Returntype : hashref where the key is the allele sequence and then value is the HGVS string 
  Exceptions : none
  Status     : At Risk

=cut

sub hgvs_genomic {
    return _hgvs_generic(@_,'genomic');
}


=head2 hgvs_transcript

  Description: Return the strings representing the CDS-level effect of each of the alleles 
               of this variation in HGVS format
  Returntype : hashref where the key is the allele sequence and then value is the HGVS string 
  Exceptions : none
  Status     : At Risk

=cut

sub hgvs_transcript {
    return _hgvs_generic(@_,'transcript');
}

=head2 hgvs_protein

  Description: Return the strings representing the protein-level effect of each of the alleles 
               of this variation in HGVS format
  Returntype : hashref where the key is the allele sequence and then value is the HGVS string 
  Exceptions : none
  Status     : At Risk

=cut

sub hgvs_protein {
    return _hgvs_generic(@_,'protein');
}

=head 

# We haven't implemented support for these methods yet

sub hgvs_rna {
    return _hgvs_generic(@_,'rna');
}

sub hgvs_mitochondrial {
    return _hgvs_generic(@_,'mitochondrial');
}
=cut

sub _hgvs_generic {
    my $self = shift;
    my $reference = pop;
    my $hgvs = shift;
    
    #The rna and mitochondrial modes have not yet been implemented, so return undef in case we get a call to these
    return undef if ($reference =~ m/rna|mitochondrial/);
    
    # The HGVS subroutine
    my $sub = qq{hgvs_$reference};
    
    # Loop over the TranscriptVariationAllele objects associated with this TranscriptVariation
    foreach my $tv_allele (@{ $self->get_all_alternate_TranscriptVariationAlleles }) {
        
        #If an HGVS hash was supplied and the allele exists as key, set the HGVS notation for this allele
        if (defined($hgvs)) {
            my $notation = $hgvs->{$tv_allele->variation_feature_seq()};
            $tv_allele->$sub($notation) if defined $notation;
        }
        # Else, add the HGVS notation for this allele to the HGVS hash
        else {
	    $hgvs->{$tv_allele->variation_feature_seq()} = $tv_allele->$sub();
        }
    }
    
    return $hgvs;
}

sub _prefetch_for_vep {
    my $self = shift;
    
    $self->cdna_coords;
    $self->cds_coords;
    $self->translation_coords;
    $self->pep_allele_string;
}

sub _get_ref_allele {
  my $self = shift;
  my ($ref_feature, $vf, $no_ref_check, $use_feature_ref) = @_;

  # this method is explicitly for retrieving the ref allele from the transcript sequence
  # which may differ from the underlying slice
  # if we don't want this, default to using the SUPER method
  return $self->SUPER::_get_ref_allele(@_) unless $use_feature_ref;

  # insertion
  return '-' if $vf->{start} > $vf->{end};

  # we only want to do this if the variant maps completely to the spliced cDNA
  # we need to check the full mapping set retrieved from the TranscriptMapper
  # in case the variant spans a whole intron, for example
  my $cdna_coords = $self->cdna_coords();
  my ($cdna_start, $cdna_end) = ($self->cdna_start, $self->cdna_end);

  return $self->SUPER::_get_ref_allele(@_) unless scalar @$cdna_coords <= 2 && $cdna_start && $cdna_end;

  my $tr = $self->transcript;
  my $tr_spliced_seq = $tr->{_variation_effect_feature_cache}->{spliced_seq} ||= $tr->spliced_seq;
  my $ref_allele = substr($tr_spliced_seq, $cdna_start - 1, ($cdna_end - $cdna_start) + 1);

  reverse_comp(\$ref_allele) unless $vf->strand eq $tr->strand;

  return $ref_allele;
}


1;

=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

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
        @{ $tv->get_all_TranscriptVariationAlleles }, "\n";

=head1 DESCRIPTION

A TranscriptVariation object represents a variation feature which is in close
proximity to an Ensembl transcript. A TranscriptVariation object has several
attributes which define the relationship of the variation to the transcript.

=cut

package Bio::EnsEMBL::Variation::TranscriptVariation;

use strict;
use warnings;

use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref check_ref);
use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap within_cds);

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlap);

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
  Status     : At Risk

=cut 

sub new {
    my $class = shift;

    my %args = @_;

    # swap a '-transcript' argument for a '-feature' one for the superclass

    for my $arg (keys %args) {
        if (lc($arg) eq '-transcript') {
            $args{'-feature'} = delete $args{$arg};
        }
    }

    # call the superclass constructor
    my $self = $class->SUPER::new(%args) || return undef;

    # rebless the alleles from vfoas to tvas
    map { bless $_, 'Bio::EnsEMBL::Variation::TranscriptVariationAllele' } 
        @{ $self->get_all_TranscriptVariationAlleles };
    
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
  Status     : At Risk

=cut

sub add_TranscriptVariationAllele {
    my ($self, $tva) = @_;
    assert_ref($tva, 'Bio::EnsEMBL::Variation::TranscriptVariationAllele');
    return $self->SUPER::add_VariationFeatureOverlapAllele($tva);
}

=head2 get_reference_TranscriptVariationAllele

  Description: Get the object representing the reference allele of this TranscriptVariation
  Returntype : Bio::EnsEMBL::Variation::TranscriptVariationAllele instance
  Exceptions : none
  Status     : At Risk

=cut

sub get_reference_TranscriptVariationAllele {
    my $self = shift;
    return $self->SUPER::get_reference_VariationFeatureOverlapAllele(@_);
}

=head2 get_all_alternate_TranscriptVariationAlleles

  Description: Get a list of the alternate alleles of this TranscriptVariation
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptVariationAllele objects
  Exceptions : none
  Status     : At Risk

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
  Status     : At Risk

=cut

sub get_all_TranscriptVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_VariationFeatureOverlapAlleles(@_);
}

=head2 transcript_stable_id

  Description: Returns the stable_id of the associated Transcript
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub transcript_stable_id {
    my $self = shift;
    return $self->SUPER::_feature_stable_id(@_);
}

=head2 transcript

  Arg [1]    : (optional) Bio::EnsEMBL::Transcript
  Description: Get/set the associated Transcript
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : throws if argument is wrong type
  Status     : At Risk

=cut

sub transcript {
    my ($self, $transcript) = @_;
    assert_ref($transcript, 'Bio::EnsEMBL::Transcript') if $transcript;
    return $self->SUPER::feature($transcript, 'Transcript');
}

=head2 feature

  Arg [1]    : (optional) Bio::EnsEMBL::Transcript
  Description: Get/set the associated Transcript (overriding the superclass feature method)
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : throws if argument is wrong type
  Status     : At Risk

=cut

sub feature {
    my $self = shift;
    return $self->transcript(@_);
}

=head2 cdna_start

  Arg [1]    : (optional) int $start
  Example    : $cdna_start = $tv->cdna_start;
  Description: Getter/Setter for the start position of this variation on the
               transcript in cDNA coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cdna_start {
    my ($self, $cdna_start) = @_;
    
    $self->{cdna_start} = $cdna_start if defined $cdna_start;
    
    unless (exists $self->{cdna_start}) {
        my $cdna_coords = $self->cdna_coords;
        
        if (@$cdna_coords != 1 || $cdna_coords->[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
            $self->{cdna_start} = undef;
            $self->{cdna_end}   = undef;
        }
        else {
            $self->{cdna_start} = $cdna_coords->[0]->start;
            $self->{cdna_end}   = $cdna_coords->[0]->end;
        }
    }
    
    return $self->{cdna_start};
}

=head2 cdna_end

  Arg [1]    : (optional) int $end
  Example    : $cdna_end = $tv->cdna_end;
  Description: Getter/Setter for the end position of this variation on the
               transcript in cDNA coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cdna_end {
    my ($self, $cdna_end) = @_;
    
    $self->{cdna_end} = $cdna_end if defined $cdna_end;
    
    # call cdna_start to calculate the start and end
    $self->cdna_start unless exists $self->{cdna_end};
    
    return $self->{cdna_end};
}

=head2 cds_start

  Arg [1]    : (optional) int $start
  Example    : $cds_start = $tv->cds_start;
  Description: Getter/Setter for the start position of this variation on the
               transcript in CDS coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cds_start {
    my ($self, $cds_start) = @_;
    
    $self->{cds_start} = $cds_start if defined $cds_start;
    
    unless (exists $self->{cds_start}) {
        my $cds_coords = $self->cds_coords;
        
        if (@$cds_coords != 1 || $cds_coords->[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
            $self->{cds_start}  = undef;
            $self->{cds_end}    = undef;
        }
        else {
            my $exon_phase = $self->transcript->start_Exon->phase;
            
            $self->{cds_start} = $cds_coords->[0]->start + ($exon_phase > 0 ? $exon_phase : 0);
            $self->{cds_end}   = $cds_coords->[0]->end   + ($exon_phase > 0 ? $exon_phase : 0);
        }
    }
    
    return $self->{cds_start};
}

=head2 cds_end

  Arg [1]    : (optional) int $end
  Example    : $cds_end = $tv->cds_end;
  Description: Getter/Setter for the end position of this variation on the
               transcript in CDS coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cds_end {
    my ($self, $cds_end) = @_;
    
    $self->{cds_end} = $cds_end if defined $cds_end;
    
    # call cds_start to calculate the start and end
    $self->cds_start unless exists $self->{cds_end};
    
    return $self->{cds_end};
}

=head2 translation_start

  Arg [1]    : (optional) int $start
  Example    : $translation_start = $tv->translation_start;
  Description: Getter/Setter for the start position of this variation on the
               transcript in peptide coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub translation_start {
    my ($self, $translation_start) = @_;
    
    $self->{translation_start} = $translation_start if defined $translation_start;
    
    unless (exists $self->{translation_start}) {
        my $translation_coords = $self->translation_coords;
        
        if (@$translation_coords != 1 || $translation_coords->[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
            $self->{translation_start} = undef;
            $self->{translation_end}   = undef;
        }
        else {
            $self->{translation_start} = $translation_coords->[0]->start;
            $self->{translation_end}   = $translation_coords->[0]->end;
        }
    }

    return $self->{translation_start};
}


=head2 translation_end

  Arg [1]    : (optional) int $end
  Example    : $transaltion_end = $tv->translation_end;
  Description: Getter/Setter for the end position of this variation on the
               transcript in peptide coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub translation_end {
    my ($self, $translation_end) = @_;
    
    $self->{translation_end} = $translation_end if defined $translation_end;
    
    # call translation_start to calculate the start and end
    $self->translation_start unless exists $self->{translation_end};
    
    return $self->{translation_end};
}

=head2 cdna_coords

  Description: Use the TranscriptMapper to calculate the cDNA 
               coordinates of this variation
  Returntype : listref of Bio::EnsEMBL::Coordinate and Bio::EnsEMBL::Gap objects
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cdna_coords {
    my $self = shift;
    
    unless ($self->{_cdna_coords}) {
        my $vf   = $self->variation_feature;
        my $tran = $self->transcript; 
        $self->{_cdna_coords} = [ $self->_mapper->genomic2cdna($vf->start, $vf->end, $tran->strand) ];
    }
    
    return $self->{_cdna_coords};
}

=head2 cds_coords

  Description: Use the TranscriptMapper to calculate the CDS 
               coordinates of this variation
  Returntype : listref of Bio::EnsEMBL::Coordinate and Bio::EnsEMBL::Gap objects
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cds_coords {
    my $self = shift;
    
    unless ($self->{_cds_coords}) {
        my $vf   = $self->variation_feature;
        my $tran = $self->transcript; 
        $self->{_cds_coords} = [ $self->_mapper->genomic2cds($vf->start, $vf->end, $tran->strand) ];
    }
    
    return $self->{_cds_coords};
}

=head2 translation_coords

  Description: Use the TranscriptMapper to calculate the peptide
               coordinates of this variation
  Returntype : listref of Bio::EnsEMBL::Coordinate and Bio::EnsEMBL::Gap objects
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub translation_coords {
    my $self = shift;
    
    unless ($self->{_translation_coords}) {
        my $vf   = $self->variation_feature;
        my $tran = $self->transcript; 
        $self->{_translation_coords} = [ $self->_mapper->genomic2pep($vf->start, $vf->end, $tran->strand) ];
    }
    
    return $self->{_translation_coords};
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
  Status     : At Risk

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

=head2 affects_cds

  Description: Check if any of this TranscriptVariation's alleles lie within
               the CDS of the Transcript
  Returntype : boolean
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub affects_cds {
    my $self = shift;
    return scalar grep { within_cds($_) } @{ $self->get_all_alternate_TranscriptVariationAlleles }; 
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
    return scalar grep { $_->SO_term =~ /stop|non_syn|frameshift|inframe|initiator/ } map {@{$_->get_all_OverlapConsequences}} @{ $self->get_all_alternate_TranscriptVariationAlleles }; 
}

=head2 get_overlapping_ProteinFeatures

  Description: Find any ProteinFeatures (e.g. pfam or interpro domains etc.) that
               the associated variation feature lies in
  Returntype : listref of Bio::EnsEMBL::ProteinFeatures (possibly empty)
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub get_overlapping_ProteinFeatures {
    my $self = shift;

    unless (exists $self->{_protein_features}) {

        $self->{_protein_features } = [];

        my $tl = $self->transcript->translation;

        if (defined $tl) {
            
            my $tl_start = $self->translation_start;
            my $tl_end   = $self->translation_end;

            if (defined $tl_start && defined $tl_end) {
                for my $feat (@{ $tl->get_all_ProteinFeatures }) {
                    if (overlap($feat->start, $feat->end, $tl_start, $tl_end)) { 
                        push @{ $self->{_protein_features} }, $feat;
                    }
                }
            }
        }
    }

    return $self->{_protein_features};
}

=head2 exon_number

  Description: Identify which exon this variant falls in   
  Returntype : '/'-separated string containing the exon number and the total 
               number of exons in this transcript, or undef if this variant 
               does not fall in an exon
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub exon_number {
    my $self = shift;
    $self->_exon_intron_number unless exists $self->{exon_number};
    return $self->{exon_number};
}

=head2 intron_number
  
  Description: Identify which intron this variant falls in   
  Returntype : '/'-separated string containing the intron number and the total 
               number of introns in this transcript, or undef if this variant 
               does not fall in an intron
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub intron_number {
    my $self = shift;
    $self->_exon_intron_number unless exists $self->{intron_number};
    return $self->{intron_number};
}

sub _exon_intron_number {
    my $self = shift;

    # work out which exon or intron this variant falls in

    # ensure the keys exist so even if we don't fall in an exon 
    # or intron we'll only call this method once

    $self->{exon_number} = $self->{intron_number} = undef;

    my $vf = $self->variation_feature;    
    
    my $vf_start = $vf->start;
    my $vf_end   = $vf->end;

    my $strand = $self->transcript->strand;

    my $exons = $self->_exons;

    my $tot_exons = scalar(@$exons);

    my $exon_count = 0;

    my $prev_exon;

    for my $exon (@$exons) {

        $exon_count++;
        
        if (overlap($vf_start, $vf_end, $exon->start, $exon->end)) {
            $self->{exon_number} = sprintf "%d/%d", $exon_count, $tot_exons;
            last;
        }

        if ($prev_exon) {
            my $intron_start = $strand == 1 ? $prev_exon->end + 1 : $exon->end + 1;
            my $intron_end   = $strand == 1 ? $exon->start - 1 : $prev_exon->start - 1;

            if ($prev_exon && overlap($vf_start, $vf_end, $intron_start, $intron_end)) {
                $self->{intron_number} = sprintf "%d/%d", $exon_count - 1, $tot_exons - 1;
                last;
            }
        }

        $prev_exon = $exon;
    }
}

sub _intron_effects {
    my $self = shift;

    # internal method used by Bio::EnsEMBL::Variation::Utils::VariationEffect
    # when calculating various consequence types
    
    # this method is a major bottle neck in the effect calculation code so 
    # we cache results and use local variables instead of method calls where
    # possible to speed things up - caveat bug-fixer!
    
    unless ($self->{_intron_effects}) {
        
        my $vf = $self->variation_feature;
        
        my $intron_effects = {};
        
        my $found_effect = 0;
        
        my $vf_start = $vf->start;
        my $vf_end   = $vf->end;

        my $insertion = $vf_start == $vf_end+1;

        for my $intron (@{ $self->_introns }) {

            my $intron_start = $intron->start;
            my $intron_end   = $intron->end;
            
            # under various circumstances the genebuild process can introduce
            # artificial short (<= 12 nucleotide) introns into transcripts
            # (e.g. to deal with errors in the reference sequence etc.), we
            # don't want to categorise variations that fall in these introns
            # as intronic, or as any kind of splice variant
            
            my $frameshift_intron = ( abs($intron_end - $intron_start) <= 12 );
            
            if ($frameshift_intron) {
                if (overlap($vf_start, $vf_end, $intron_start, $intron_end)) {
                    $intron_effects->{within_frameshift_intron} = 1;
                    next;
                }
            }

            if (overlap($vf_start, $vf_end, $intron_start, $intron_start+1)) {
                $intron_effects->{start_splice_site} = 1;
            }
            
            if (overlap($vf_start, $vf_end, $intron_end-1, $intron_end)) {
                $intron_effects->{end_splice_site} = 1;
            }
            
            # we need to special case insertions between the donor and acceptor sites

            if (overlap($vf_start, $vf_end, $intron_start+2, $intron_end-2) or 
                ($insertion && ($vf_start == $intron_start+2 || $vf_end == $intron_end-2)) ) {
                $intron_effects->{intronic} = 1;
            }
            
            # the definition of splice_region (SO:0001630) is "within 1-3 bases 
            # of the exon or 3-8 bases of the intron", the intron start is the 
            # first base of the intron so we only need to add or subtract 7 from 
            # it to get the correct coordinate. We also need to special case 
            # insertions between the edge of an exon and a donor or acceptor site
            # and between a donor or acceptor site and the intron
            
            if ( overlap($vf_start, $vf_end, $intron_start-3, $intron_start-1) or
                 overlap($vf_start, $vf_end, $intron_start+2, $intron_start+7) or
                 overlap($vf_start, $vf_end, $intron_end-7,   $intron_end-2  ) or
                 overlap($vf_start, $vf_end, $intron_end+1,   $intron_end+3  ) or
                 ($insertion && ( 
                     $vf_start == $intron_start || 
                     $vf_end == $intron_end ||
                     $vf_start == $intron_start+2 ||
                     $vf_end == $intron_end-2
                    ) )) { 
                   
                $intron_effects->{splice_region} = 1;
            }
        }
        
        $self->{_intron_effects} = $intron_effects;       
    }

    return $self->{_intron_effects};
}

# NB: the methods below all cache their data in the associated transcript itself, this
# gives a significant speed up when you are calculating the effect of all variations
# on a transcript, and means that the cache will be freed when the transcript itself
# is garbage collected rather than us having to maintain a transcript feature cache 
# ourselves

sub _introns {
    my $self = shift;
    
    my $tran = $self->transcript;
    
    my $introns = $tran->{_variation_effect_feature_cache}->{introns} ||= $tran->get_all_Introns;
    
    return $introns;
}

sub _exons {
    my $self = shift;

    my $tran = $self->transcript;

    my $exons = $tran->{_variation_effect_feature_cache}->{exons} ||= $tran->get_all_Exons;

    return $exons;
}

sub _translateable_seq {
    my $self = shift;
    
    my $tran = $self->transcript;
    
    my $tran_seq = $tran->{_variation_effect_feature_cache}->{translateable_seq} ||= $tran->translateable_seq;
    
    return $tran_seq;
}

sub _mapper {
    my $self = shift;
    
    my $tran = $self->transcript;
    
    my $mapper = $tran->{_variation_effect_feature_cache}->{mapper} ||= $tran->get_TranscriptMapper;
    
    return $mapper;
}
  
sub _peptide {
    my $self = shift;
    
    my $tran = $self->transcript;
    
    my $peptide = $tran->{_variation_effect_feature_cache}->{peptide};
    
    unless ($peptide) {
        my $translation = $tran->translate;
        $peptide = $translation ? $translation->seq : undef;
        $tran->{_variation_effect_feature_cache}->{peptide} = $peptide;
    }
    
    return $peptide;
}

sub _translation_md5 {
    my $self = shift;

    my $tran = $self->transcript;
    
    unless (exists $tran->{_variation_effect_feature_cache}->{translation_md5}) {
        $tran->{_variation_effect_feature_cache}->{translation_md5} = 
            $self->_peptide ? md5_hex($self->_peptide) : undef;
    }

    return $tran->{_variation_effect_feature_cache}->{translation_md5};
}

sub _codon_table {
    my $self = shift;
    
    my $tran = $self->transcript;
    
    my $codon_table = $tran->{_variation_effect_feature_cache}->{codon_table};
    
    unless ($codon_table) {
        # for mithocondrial dna we need to to use a different codon table
        my $attrib = $tran->slice->get_all_Attributes('codon_table')->[0]; 
        
        # default to the vertebrate codon table which is denoted as 1
        $codon_table = $attrib ? $attrib->value : 1;
        
        $tran->{_variation_effect_feature_cache}->{codon_table} = $codon_table
    }
    
    return $codon_table;
}

sub _protein_function_predictions {
    
    my ($self, $analysis) = @_;

    my $tran = $self->transcript;

    my $matrix = $tran->{_variation_effect_feature_cache}->{protein_function_predictions}->{$analysis};

    unless ($matrix || exists($tran->{_variation_effect_feature_cache}->{protein_function_predictions}->{$analysis})) {
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

=head2 hgvs_coding

  Description: Return the strings representing the CDS-level effect of each of the alleles 
               of this variation in HGVS format
  Returntype : hashref where the key is the allele sequence and then value is the HGVS string 
  Exceptions : none
  Status     : At Risk

=cut

sub hgvs_coding {
    deprecate('HGVS coding support has been moved to hgvs_transcript. This method will be removed in the next release.');
    return _hgvs_generic(@_,'transcript');
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
    
    #ÊThe rna and mitochondrial modes have not yet been implemented, so return undef in case we get a call to these
    return undef if ($reference =~ m/rna|mitochondrial/);
    
    # The HGVS subroutine
    my $sub = qq{hgvs_$reference};
    
    # Loop over the TranscriptVariationAllele objects associated with this TranscriptVariation
    foreach my $tv_allele (@{ $self->get_all_alternate_TranscriptVariationAlleles }) {
        
        #ÊIf an HGVS hash was supplied and the allele exists as key, set the HGVS notation for this allele
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


1;

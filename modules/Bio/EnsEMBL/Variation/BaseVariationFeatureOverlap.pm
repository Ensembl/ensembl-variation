package Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap within_cds);

sub new {
    my $class = shift;

    my (
        $base_variation_feature,
        $feature, 
        $no_transfer
    ) = rearrange([qw(
            BASE_VARIATION_FEATURE
            FEATURE
            NO_TRANSFER
        )], @_);
}

sub new_fast {
    my ($class, $hashref) = @_;
    return bless $hashref, $class;
}

sub feature {
    my ($self, $feature) = @_;
    $self->{feature} = $feature if $feature;
    return $self->{feature};
}

sub base_variation_feature {
    my ($self, $bvf) = @_;
    $self->{base_variation_feature} = $bvf if $bvf;
    return $self->{base_variation_feature};
}

=head2 add_BaseVariationFeatureOverlapAllele

  Arg [1]    : A Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele instance
  Description: Add an allele to this BaseVariationFeatureOverlap
  Returntype : none
  Exceptions : throws if the argument is not the expected type
  Status     : At Risk

=cut

sub add_BaseVariationFeatureOverlapAllele {
    my ($self, $vfoa) = @_;

    assert_ref($vfoa, 'Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele');

    if ($vfoa->is_reference) {
        $self->{reference_allele} = $vfoa;
    }
    else {
        my $alt_alleles = $self->{alt_alleles} ||= [];
        push @$alt_alleles, $vfoa;
    }
}

=head2 get_reference_BaseVariationFeatureOverlapAllele

  Description: Get the object representing the reference allele of this BaseVariationFeatureOverlapAllele
  Returntype : Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele instance
  Exceptions : none
  Status     : At Risk

=cut

sub get_reference_BaseVariationFeatureOverlapAllele {
    my $self = shift;
    return $self->{reference_allele};
}

=head2 get_all_alternate_BaseVariationFeatureOverlapAlleles

  Description: Get a list of the alternate alleles of this BaseVariationFeatureOverlapAllele
  Returntype : listref of Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele objects
  Exceptions : none
  Status     : At Risk

=cut

sub get_all_alternate_BaseVariationFeatureOverlapAlleles {
    my $self = shift;

    $self->{alt_alleles} ||= [];
    
    return $self->{alt_alleles};
}

=head2 get_all_BaseVariationFeatureOverlapAlleles

  Description: Get a list of the all the alleles, both reference and alternate, of this
               BaseVariationFeatureOverlap
  Returntype : listref of Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele objects
  Exceptions : none
  Status     : At Risk

=cut

sub get_all_BaseVariationFeatureOverlapAlleles {
    my $self = shift;
	my @alts = @{ $self->get_all_alternate_BaseVariationFeatureOverlapAlleles };
	my $ref = $self->get_reference_BaseVariationFeatureOverlapAllele;
	unshift @alts, $ref if defined($ref);
	
    return \@alts;
}

=head2 consequence_type

  Arg [1]    : (optional) String $term_type
  Description: Get a list of all the unique consequence terms of the alleles of this 
               VariationFeatureOverlap. By default returns Ensembl display terms
               (e.g. 'NON_SYNONYMOUS_CODING'). $term_type can also be 'label'
               (e.g. 'Non-synonymous coding'), 'SO' (Sequence Ontology, e.g.
               'non_synonymous_codon') or 'NCBI' (e.g. 'missense')
  Returntype : listref of strings
  Exceptions : none
  Status     : At Risk

=cut

sub consequence_type {
    my $self = shift;
	my $term_type = shift;
	
	my $method_name;
	
    # delete cached term
    if(defined($term_type)) {
        delete $self->{_consequence_type};
		$method_name = $term_type.($term_type eq 'label' ? '' : '_term');
		$method_name = 'display_term' unless defined $self->most_severe_OverlapConsequence && $self->most_severe_OverlapConsequence->can($method_name);
    }
	
	$method_name ||= 'display_term';
    
    unless ($self->{_consequence_type}) {
        
        # use a hash to ensure we don't include redundant terms (because more than one
        # allele may have the same consequence display_term)

        my %cons_types;

        for my $allele (@{ $self->get_all_alternate_BaseVariationFeatureOverlapAlleles }) {
            for my $cons (@{ $allele->get_all_OverlapConsequences }) {
                $cons_types{$cons->$method_name} = $cons->rank;
            }
        }

        # sort the consequence types by rank such that the more severe terms are earlier in the list

        $self->{_consequence_type} = [ sort { $cons_types{$a} <=> $cons_types{$b} } keys %cons_types ];
    }
    
    return $self->{_consequence_type};
}

=head2 most_severe_OverlapConsequence

  Description: Get the OverlapConsequence considered (by Ensembl) to be the most severe 
               consequence of all the alleles of this VariationFeatureOverlap 
  Returntype : Bio::EnsEMBL::Variation::OverlapConsequence
  Exceptions : none
  Status     : At Risk

=cut

sub most_severe_OverlapConsequence {
    my $self = shift;
    
    unless ($self->{_most_severe_consequence}) {
        
        my $highest;
        
        for my $allele (@{ $self->get_all_alternate_BaseVariationFeatureOverlapAlleles }) {
            for my $cons (@{ $allele->get_all_OverlapConsequences }) {
                $highest ||= $cons;
                if ($cons->rank < $highest->rank) {
                    $highest = $cons;
                }
            }
        }
        
        $self->{_most_severe_consequence} = $highest;
    }
    
    return $self->{_most_severe_consequence};
}

=head2 display_consequence

  Arg [1]    : (optional) String $term_type
  Description: Get the term for the most severe OverlapConsequence of this 
               VariationFeatureOverlap. By default returns Ensembl display terms
               (e.g. 'NON_SYNONYMOUS_CODING'). $term_type can also be 'label'
               (e.g. 'Non-synonymous coding'), 'SO' (Sequence Ontology, e.g.
               'non_synonymous_codon') or 'NCBI' (e.g. 'missense')
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub display_consequence {
    my $self = shift;
	my $term_type = shift;
	
	my $method_name;
	
    # delete cached term
    if(defined($term_type)) {
		$method_name = $term_type.($term_type eq 'label' ? '' : '_term');
		$method_name = 'display_term' unless @{$self->get_all_OverlapConsequences} && $self->get_all_OverlapConsequences->[0]->can($method_name);
    }
	
	$method_name ||= 'display_term';
    
    my $worst_conseq = $self->most_severe_OverlapConsequence;

    return $worst_conseq ? $worst_conseq->$method_name : '';
}

sub _intron_effects {
    my $self = shift;

    # internal method used by Bio::EnsEMBL::Variation::Utils::VariationEffect
    # when calculating various consequence types
    
    # this method is a major bottle neck in the effect calculation code so 
    # we cache results and use local variables instead of method calls where
    # possible to speed things up - caveat bug-fixer!
    
    unless ($self->{_intron_effects}) {
        
        my $vf = $self->base_variation_feature;
        
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
    
    my $feat = $self->feature;

    if ($feat->isa('Bio::EnsEMBL::Transcript')) {

        my $introns = $feat->{_variation_effect_feature_cache}->{introns} ||= $feat->get_all_Introns;

        return $introns;
    
    }

    return [];

}



1;


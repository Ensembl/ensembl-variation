=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Variation::VariationFeatureOverlap

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::VariationFeatureOverlap;

    my $vfo = Bio::EnsEMBL::Variation::VariationFeatureOverlap->new(
        -feature            => $feature,
        -variation_feature  => $var_feat
    );

    print "consequence type: ", (join ",", @{ $vfo->consequence_type }), "\n";
    print "most severe consequence: ", $vfo->display_consequence, "\n";

=head1 DESCRIPTION

A VariationFeatureOverlap represents a VariationFeature which is in close
proximity to another Ensembl Feature. It is the superclass of feature-specific
objects such as TranscriptVariation and RegulatoryFeatureVariation, and has
methods common to all such objects. You will not normally instantiate this
class directly, instead instantiating one of the feature-specific subclasses.

=cut

package Bio::EnsEMBL::Variation::VariationFeatureOverlap;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(expand);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

=head2 new

  Arg [-FEATURE] : 
    The Bio::EnsEMBL::Feature associated with the given VariationFeature

  Arg [-VARIATION_FEATURE] :
    The Bio::EnsEMBL::VariationFeature associated with the given Feature

  Arg [-ADAPTOR] :
    A Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor

  Arg [-DISAMBIGUATE_SINGLE_NUCLEOTIDE_ALLELES] :
    A flag indiciating if ambiguous single nucleotide alleles should be disambiguated
    when constructing the VariationFeatureOverlapAllele objects, e.g. a Variationfeature
    with an allele string like 'T/M' would be treated as if it were 'T/A/C'. We limit
    ourselves to single nucleotide alleles to avoid the combinatorial explosion if we
    allowed longer alleles with potentially many ambiguous bases.

  Example : 
    my $vfo = Bio::EnsEMBL::Variation::VariationFeatureOverlap->new(
        -feature           => $feature,
        -variation_feature => $var_feat
    );

  Description: Constructs a new VariationFeatureOverlap instance given a VariationFeature
               and a Feature
  Returntype : A new Bio::EnsEMBL::Variation::VariationFeatureOverlap instance 
  Exceptions : throws unless both VARIATION_FEATURE and FEATURE are supplied, or if the
               supplied ADAPTOR is the wrong class
  Status     : At Risk

=cut 

sub new {

    my $class = shift;

    my (
        $variation_feature, 
        $feature, 
        $adaptor, 
        $ref_feature, 
        $disambiguate_sn_alleles,
        $no_ref_check,
        $no_transfer
    ) = rearrange([qw(
            VARIATION_FEATURE 
            FEATURE 
            ADAPTOR 
            REF_FEATURE 
            DISAMBIGUATE_SINGLE_NUCLEOTIDE_ALLELES
            NO_REF_CHECK
            NO_TRANSFER
        )], @_);

    assert_ref($variation_feature, 'Bio::EnsEMBL::Variation::VariationFeature');
    assert_ref($feature, 'Bio::EnsEMBL::Feature');
    assert_ref($adaptor, 'Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor') if $adaptor;

    $ref_feature ||= $variation_feature->slice;

    # we need to ensure the Feature and the VariationFeature live on the same slice
    # so we explicitly transfer the Feature here
    unless(defined $no_transfer && $no_transfer == 1) {
        $feature = $feature->transfer($variation_feature->slice) 
            or throw("Unable to transfer the supplied feature to the same slice as the variation feature");
    }

    my $self = bless {
        variation_feature   => $variation_feature,
        feature             => $feature,
        adaptor             => $adaptor,
        ref_feature         => $ref_feature,
    }, $class;
    
    
    my $ref_allele;
    
    # we take the reference allele sequence from the reference sequence, not from the allele string
    unless($no_ref_check) {
        $ref_allele = $ref_feature->subseq(
            $variation_feature->start, 
            $variation_feature->end, 
            $variation_feature->strand
        );
    }

    # get the variation feature allele string, expand it, and split it into separate alleles
    
    my $allele_string = $variation_feature->allele_string;
    
    expand(\$allele_string);

    my @alleles = split /\//, $allele_string;
    
    $ref_allele = $alleles[0] if $no_ref_check;
    $ref_allele = '-' unless $ref_allele;
  
    if ($disambiguate_sn_alleles) {
        
        # if this flag is set, disambiguate any ambiguous single nucleotide alleles, so  
        # e.g. an allele string like T/M would be equivalent to an allele string of T/A/C
        # we only do this for single nucleotide alleles to avoid the combinatorial explosion
        # of long allele strings with potentially many ambiguous bases (because ensembl 
        # genomes want this functionality)

        my @possible_alleles;

        for my $allele (@alleles) {
            
            if ($allele !~ /^[ACGT-]+$/ && length($allele) == 1) {
                for my $possible ( split //, unambiguity_code($allele) ) {
                    push @possible_alleles, $possible;
                }
            }
            else {
                # the allele is either unambiguous or longer than 1 nucleotide, so add it unaltered
                push @possible_alleles, $allele;
            }
        }

        @alleles = @possible_alleles;
    }

    # make sure the alleles are unique
    
    # we also want to deal with alleles like (T)0 which expand into 
    # an empty string and we want to treat this as a deletion, so 
    # we replace
    # any empty strings with '-'
    
    @alleles = keys %{ { map { ($_ || '-') => 1 } @alleles } };

    # create an object representing the reference allele
    
    my $ref_vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new(
        -variation_feature_overlap   => $self,
        -variation_feature_seq       => $ref_allele,
        -is_reference                => 1,
    );
    
    $self->add_VariationFeatureOverlapAllele($ref_vfoa);

    # create objects representing the alternate alleles
    
    for my $allele (@alleles) {
        
        next if $allele eq $ref_allele;
        
        my $vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new(
            -variation_feature_overlap  => $self,
            -variation_feature_seq      => $allele,
            -is_reference               => 0,
        );
       
        $self->add_VariationFeatureOverlapAllele($vfoa);
    }
    
    return $self;
}

sub dbID {
    my $self = shift;
    
    unless ($self->{dbID}) {
        # we don't really have a dbID, so concatenate all the dbIDs of our alleles

        $self->{dbID} = join '_', map { $_->dbID } @{ $self->get_all_alternate_VariationFeatureOverlapAlleles };
    }

    return $self->{dbID};
}

sub new_fast {
    my ($class, $hashref) = @_;
    return bless $hashref, $class;
}

=head2 variation_feature

  Arg [1]    : (optional) A Bio::EnsEMBL::Variation::VariationFeature
  Description: Get/set the associated VariationFeature, lazy-loading it if required
  Returntype : Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throws if the argument is the wrong type
  Status     : At Risk

=cut

sub variation_feature {
    my ($self, $variation_feature) = @_;

    if ($variation_feature) {
        assert_ref($variation_feature, 'Bio::EnsEMBL::Variation::VariationFeature');
        $self->{variation_feature} = $variation_feature;
    }
    
    if (my $vf_id = $self->{_variation_feature_id}) {
        
        # lazy-load the VariationFeature
        
        if (my $adap = $self->{adaptor}) {
            if (my $vfa = $adap->db->get_VariationFeatureAdaptor) {
                if (my $vf = $vfa->fetch_by_dbID($vf_id)) {
                    $self->{variation_feature} = $vf;
                    delete $self->{_variation_feature_id};
                }
            }
        }
    }
    
    return $self->{variation_feature};
}

sub _variation_feature_id {

    # get the dbID of the variation feature, using the VariationFeature object 
    # if we have one, or the internal hash value if we don't
    
    my $self = shift;
    
    if (my $vf = $self->{variation_feature}) {
        return $vf->dbID;
    }
    elsif (my $id = $self->{_variation_feature_id}) {
        return $id;
    }
    else {
        return undef;
    }
}

=head2 feature

  Arg [1]    : (optional) A Bio::EnsEMBL::Feature
  Description: Get/set the associated Feature, lazy-loading it if required
  Returntype : Bio::EnsEMBL::Feature
  Exceptions : throws isf the argument is the wrong type
  Status     : At Risk

=cut

sub feature {
    my ($self, $feature, $type) = @_;
    
    if ($feature) {
        assert_ref($feature, 'Bio::EnsEMBL::Feature');
        $self->{feature} = $feature;
    }
 
    if ($type && !$self->{feature}) {
    
        # try to lazy load the feature
        
        if (my $adap = $self->{adaptor}) {
            
            my $get_method = 'get_'.$type.'Adaptor';
           
            # XXX: this can doesn't work because the method is AUTOLOADed, need to rethink this...
            #if ($adap->db->dnadb->can($get_method)) {
                if (my $fa = $adap->db->dnadb->$get_method) {
                    
                    # if we have a stable id for the feature use that
                    if (my $feature_stable_id = $self->{_feature_stable_id}) {
                        if (my $f = $fa->fetch_by_stable_id($feature_stable_id)) {
                            $self->{feature} = $f;
                            delete $self->{_feature_stable_id};
                        }
                    }
                    elsif (my $feature_label = $self->{_feature_label}) {
                        # get a slice covering the vf
                        
                        #for my $f ($fa->fetch_all_by_Slice_constraint)
                    }
                }
            #}
            else {
                warn "Cannot get an adaptor for type: $type";
            }
        }
    }
    
    return $self->{feature};
}

sub _fetch_feature_for_stable_id {
    
    # we shouldn't actually need this method as there will apparently
    # soon be core support for fetching any feature by its stable id, 
    # but I'm waiting for core to add this...

    my ($self, $feature_stable_id) = @_;
    
    my $type_lookup = {
        G   => { type => 'Gene',                 group => 'core' },
        T   => { type => 'Transcript',           group => 'core'  },
        R   => { type => 'RegulatoryFeature',    group => 'funcgen' },
    };
    
    if ($feature_stable_id =~ /^ENS[A-Z]*([G|R|T])\d+$/) {
        
        my $type  = $type_lookup->{$1}->{type};
        my $group = $type_lookup->{$1}->{group};
        
        if (my $adap = $self->{adaptor}) {
            
            my $get_method = 'get_'.$type.'Adaptor';
            
            if ($adap->db->dnadb->can($get_method)) {
                if (my $fa = $adap->db->dnadb->$get_method) {
                    
                    # if we have a stable id for the feature use that
                    if (my $feature_stable_id = $self->{_feature_stable_id}) {
                        if (my $f = $fa->fetch_by_stable_id($feature_stable_id)) {
                            $self->{feature} = $f;
                            delete $self->{_feature_stable_id};
                        }
                    }
                    elsif (my $feature_label = $self->{_feature_label}) {
                        # get a slice covering the vf
                        
                        
                        #for my $f ($fa->fetch_all_by_Slice_constraint)
                    }
                }
            }
            else {
                warn "Cannot get an adaptor for type: $type";
            }
    }
    }
}

sub _fetch_adaptor_for_group {
    my ($self, $group) = @_;
    
}

sub _feature_stable_id {
    my $self = shift;
    if ($self->feature && $self->feature->can('stable_id')) {
        return $self->feature->stable_id;
    }
    elsif (my $id = $self->{_feature_stable_id}) {
        return $id;
    }
    else {
        return undef;
    }
}

sub get_VariationFeatureOverlapAllele_for_allele_seq {
    my ($self, $allele_seq) = @_;
    return $self->{_alleles_by_seq}->{$allele_seq};
}

=head2 add_VariationFeatureOverlapAllele

  Arg [1]    : A Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele instance
  Description: Add an allele to this VariationFeatureOverlap
  Returntype : none
  Exceptions : throws if the argument is not the expected type
  Status     : At Risk

=cut

sub add_VariationFeatureOverlapAllele {
    my ($self, $vfoa) = @_;

    assert_ref($vfoa, 'Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele');

    if ($vfoa->is_reference) {
        $self->{reference_allele} = $vfoa;
    }
    else {
        my $alt_alleles = $self->{alt_alleles} ||= [];
        push @$alt_alleles, $vfoa;
    }

    $self->{_alleles_by_seq}->{ $vfoa->variation_feature_seq } = $vfoa;
}

=head2 get_reference_VariationFeatureOverlapAllele

  Description: Get the object representing the reference allele of this VariationFeatureOverlapAllele
  Returntype : Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele instance
  Exceptions : none
  Status     : At Risk

=cut

sub get_reference_VariationFeatureOverlapAllele {
    my $self = shift;
    return $self->{reference_allele};
}

=head2 get_all_alternate_VariationFeatureOverlapAlleles

  Description: Get a list of the alternate alleles of this VariationFeatureOverlapAllele
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele objects
  Exceptions : none
  Status     : At Risk

=cut

sub get_all_alternate_VariationFeatureOverlapAlleles {
    my $self = shift;

    $self->{alt_alleles} ||= [];
    
    return $self->{alt_alleles};
}

=head2 get_all_VariationFeatureOverlapAlleles

  Description: Get a list of the all the alleles, both reference and alternate, of this
               VariationFeatureOverlap
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele objects
  Exceptions : none
  Status     : At Risk

=cut

sub get_all_VariationFeatureOverlapAlleles {
    my $self = shift;
    return [ 
        $self->get_reference_VariationFeatureOverlapAllele, 
        @{ $self->get_all_alternate_VariationFeatureOverlapAlleles } 
    ];
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

        for my $allele (@{ $self->get_all_alternate_VariationFeatureOverlapAlleles }) {
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
        
        for my $allele (@{ $self->get_all_alternate_VariationFeatureOverlapAlleles }) {
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

sub _convert_to_sara {
    my $self = shift;
    
    my $ref_allele = $self->{reference_allele};
    $ref_allele->_convert_to_sara;
    
    $self->{alt_alleles} = [$ref_allele];
}

sub _rearrange_alleles {
    my $self = shift;
    my $keep_alleles = shift;
    
    # fix alt alleles
    my $alt_alleles = $self->{alt_alleles};
    my @new_alleles = grep {$keep_alleles->{$_->variation_feature_seq}} @$alt_alleles;
    $self->{alt_alleles} = scalar @new_alleles ? \@new_alleles : $alt_alleles;
    
    # copy to ref allele if homozygous non-ref
    $self->{reference_allele} = $self->{alt_alleles}->[0] if scalar keys %$keep_alleles == 1;
}

1;

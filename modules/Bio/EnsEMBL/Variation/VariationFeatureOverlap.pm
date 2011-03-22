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

package Bio::EnsEMBL::Variation::VariationFeatureOverlap;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(expand);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

sub new {

    my $class = shift;

    my (
        $variation_feature, 
        $feature, 
        $adaptor, 
        $ref_feature, 
        $disambiguate_sn_alleles
    ) = rearrange([qw(
            VARIATION_FEATURE 
            FEATURE 
            ADAPTOR 
            REF_FEATURE 
            DISAMBIGUATE_SINGLE_NUCLEOTIDE_ALLELES
        )], @_);

    throw("VariationFeature argument required") unless $variation_feature;
    throw("Feature argument required") unless $feature;

    $ref_feature ||= $variation_feature->slice;

    my $self = bless {
        variation_feature   => $variation_feature,
        feature             => $feature,
        adaptor             => $adaptor,
        ref_feature         => $ref_feature,
    }, $class;
    
    # we take the reference allele sequence from the reference sequence, not from the allele string 
    
    my $ref_allele = $ref_feature->subseq(
        $variation_feature->start, 
        $variation_feature->end, 
        $variation_feature->strand
    );
    
    $ref_allele = '-' unless $ref_allele;

    # get the variation feature allele string, expand it, and split it into separate alleles
    
    my $allele_string = $variation_feature->allele_string;
    
    expand(\$allele_string);

    my @alleles = split /\//, $allele_string;
  
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
    
    @alleles = keys %{ { map { $_ => 1 } @alleles } };

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

sub new_fast {
    my ($class, $hashref) = @_;
    return bless $hashref, $class;
}

sub dbID {
    my ($self, $dbID) = @_;
    $self->{dbID} = $dbID if defined $dbID;
    return $self->{dbID};
}

sub variation_feature {
    my ($self, $variation_feature) = @_;
    
    $self->{variation_feature} = $variation_feature if $variation_feature;
    
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

sub variation_feature_id {
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

sub feature {
    my ($self, $feature, $type) = @_;
    
    $self->{feature} = $feature if $feature;
 
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

sub feature_stable_id {
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

sub add_VariationFeatureOverlapAllele {
    my ($self, $vfoa) = @_;

    unless ($vfoa && ref $vfoa && $vfoa->isa('Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele')) {
        throw("Argument must be a VariationFeatureOverlapAllele");
    }

    if ($vfoa->is_reference) {
        $self->{reference_allele} = $vfoa;
    }
    else {
        my $alt_alleles = $self->{alt_alleles} ||= [];
        push @$alt_alleles, $vfoa;
    }
}

sub get_reference_VariationFeatureOverlapAllele {
    my $self = shift;
    return $self->{reference_allele};
}

sub get_all_alternate_VariationFeatureOverlapAlleles {
    my $self = shift;
    return $self->{alt_alleles};
}

sub get_all_VariationFeatureOverlapAlleles {
    my $self = shift;
    return [ 
        $self->get_reference_VariationFeatureOverlapAllele, 
        @{ $self->get_all_alternate_VariationFeatureOverlapAlleles } 
    ];
}

sub consequence_type {
    my $self = shift;
    
    unless ($self->{_consequence_type}) {
        
        # build up a unique list of all the consequence display terms
        # of this VariationFeatureOverlap's alleles

        my %cons_types;

        for my $allele (@{ $self->get_all_alternate_VariationFeatureOverlapAlleles }) {
            for my $cons (@{ $allele->consequence_types }) {
                $cons_types{$cons->display_term}++
            }
        }
        
        $self->{_consequence_type} = [ keys %cons_types ];
    }
    
    return $self->{_consequence_type};
}

sub most_severe_consequence {
    my $self = shift;
    
    unless ($self->{_most_severe_consequence}) {
        
        my $highest;
        
        for my $allele (@{ $self->get_all_alternate_VariationFeatureOverlapAlleles }) {
            for my $cons (@{ $allele->consequence_types }) {
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

sub display_consequence {
    my $self = shift;
    return $self->most_severe_consequence->display_term;
}

1;

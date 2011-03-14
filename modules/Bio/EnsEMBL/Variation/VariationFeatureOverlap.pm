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

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(expand);
use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

sub new {

    my $class = shift;

    my ($variation_feature, $feature, $adaptor, $ref_feature) = 
        rearrange([qw(VARIATION_FEATURE FEATURE ADAPTOR REF_FEATURE)], @_);

    $ref_feature ||= $variation_feature->slice;

    my $self = bless {
        variation_feature   => $variation_feature,
        feature             => $feature,
        adaptor             => $adaptor,
        ref_feature         => $ref_feature,
    }, $class;
    
    # now look at each allele of the VariationFeature in turn
    
    # get the allele string, expand it, and split it into separate alleles
    
    my $ref_allele = $ref_feature->subseq(
        $variation_feature->start, 
        $variation_feature->end, 
        $variation_feature->strand
    );
    
    $ref_allele = '-' unless $ref_allele;
        
    my $allele_string = $variation_feature->allele_string;
    
    expand(\$allele_string);

    my @alleles = split /\//, $allele_string;
  
    # create an object representing the reference allele
    
    my $ref_vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
        variation_feature_overlap   => $self,
        variation_feature_seq       => $ref_allele,
        is_reference                => 1,
    });
    
    $self->reference_allele($ref_vfoa);

    # create objects representing the alternate alleles
    
    my @alt_alleles;
    
    for my $allele (@alleles) {
        
        next if $allele eq $ref_allele;
        
        my $vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
            variation_feature_overlap   => $self,
            variation_feature_seq       => $allele,
        });
        
        push @alt_alleles, $vfoa,
    }
    
    $self->alt_alleles(@alt_alleles);
    
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
    
    return $self->{feature};
}

sub _fetch_feature_for_stable_id {
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
    if ($self->{feature} && $self->{feature}->can('stable_id')) {
        return $self->{feature}->stable_id;
    }
    elsif (my $id = $self->{_feature_stable_id}) {
        return $id;
    }
    else {
        return undef;
    }
}

sub reference_allele {
    my ($self, $reference_allele) = @_;
    $self->{reference_allele} = $reference_allele if $reference_allele;
    return $self->{reference_allele};
}

sub alt_alleles {
    my ($self, @new_alt_alleles) = @_;
    
    my $alt_alleles = $self->{alt_alleles} ||= [];

    push @$alt_alleles, @new_alt_alleles if @new_alt_alleles;

    return $alt_alleles;
}

sub alleles {
    my ($self) = @_;
    return [ $self->reference_allele, @{ $self->alt_alleles } ];
}

sub consequence_type {
    my $self = shift;
    
    unless ($self->{_consequence_type}) {
        
        my %cons_types;

        for my $allele (@{ $self->alt_alleles }) {
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
        
        for my $allele (@{ $self->alt_alleles }) {
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

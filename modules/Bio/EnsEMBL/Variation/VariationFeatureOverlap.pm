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

use Bio::EnsEMBL::Utils::Sequence qw(expand);
use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

sub new {
    my ($class, $hashref) = @_;
    
    my $self = bless $hashref, $class;
    
    # now look at each allele of the VariationFeature in turn
    
    # get the allele string, expand it, and split it into separate alleles
    
    my $vf          = $self->{variation_feature};
    my $tran        = $self->{transcript};
    my $ref_feature = $self->{ref_feature} || $vf->slice;
    
    my $ref_allele = $ref_feature->subseq($vf->start, $vf->end, $vf->strand);
    
    $ref_allele = '-' unless $ref_allele;
        
    my $allele_string = $vf->allele_string;
    
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

sub overlap_consequences {
    my ($self, $overlap_consequences) = @_;
    
    $self->{overlap_consequences} = $overlap_consequences if $overlap_consequences;
    
    unless ($self->{overlap_consequences}) {
        
        # try to load the consequence objects from the database
        
        my @cons;
        
        # get an adaptor either from us, or from the associated variation feature
        #if (my $adap = $self->{adaptor} || $self->variation_feature->{adaptor}) {
        if (my $adap = $self->{adaptor}) {
            
            # get the list of possible overlap consequences
            #for my $cons (@{ $adap->db->get_AttributeAdaptor->fetch_all_OverlapConsequences }) {
            for my $cons (values %{ $adap->_overlap_consequences }) {
                
                # check that this consequence type applies to this feature type
                my $ens_classes = $adap->ensembl_classes_for_SO_term($cons->feature_SO_term);
                
                my $feat_class = ref $self->feature;
                
                if (grep { $_ eq $feat_class } @{ $ens_classes }) {
                    
#                    # also check if the biotypes match (or if the biotype is not defined)
#                    my $biotype = $adap->ensembl_biotype_for_SO_term($cons->feature_SO_term);
#                    
#                    if (defined $biotype && $self->feature->can('biotype')) {
#                        #next unless $self->feature->biotype eq $biotype;
#                    }
                    
                    # OK, this consequence type applies to this feature
                    push @cons, $cons;
                }   
            }
        }
        else {
            warn "Can't load OverlapConsequence objects without an adaptor";
        }
        
        $self->{overlap_consequences} = \@cons;
    }
    
    return $self->{overlap_consequences};
}

sub consequence_type {
    my $self = shift;
    
    unless ($self->{_consequence_type}) {
        my @cons = map { 
            map { $_->ensembl_term } @{ $_->consequence_types } 
        } @{ $self->alt_alleles };
        
        $self->{_consequence_type} = \@cons;
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
    return $self->most_severe_consequence->ensembl_term;
}

1;

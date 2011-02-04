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

my $feature_SO_term_to_ensembl_type = {
    'mRNA'                  =>  {
        feature_type    => 'Bio::EnsEMBL::Transcript',
        variant_type    => 'Bio::EnsEMBL::Variation::TranscriptVariationNew',
    },
    'regulatory_region'     => {
        feature_type    => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
        variant_type    => 'Bio::EnsEMBL::Variation::RegualtoryFeatureVariation',
    },
    'binding_site'     => {
        feature_type    => 'Bio::EnsEMBL::Funcgen::MotifFeature',
        variant_type    => 'Bio::EnsEMBL::Variation::MotifFeatureVariation',
    },
};

sub new {
    my ($class, $hashref) = @_;
    
    my $self = bless $hashref, $class;
    
    # now look at each allele of the VariationFeature in turn
    
    # get the allele string, expand it, and split it into separate alleles
    
    my $vf   = $self->{variation_feature};
    my $tran = $self->{transcript};
    
    my $allele_string = $vf->allele_string;
    
    expand(\$allele_string);
    
    unless ($allele_string =~ /\//) {
        # for the HGMDs and CNV probes just set the alt allele 
        # to the reference allele
        $allele_string .= "/$allele_string";
    }
    
    my @alleles = split /\//, $allele_string;
  
    # create an object representing the reference allele
    
    my $ref_allele = shift @alleles;

    my $ref_vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
        variation_feature_overlap   => $self,
        variation_feature_seq       => $ref_allele,
        is_reference                => 1,
    });
    
    $self->reference_allele($ref_vfoa);

    # create objects representing the alternate alleles
    
    my @alt_alleles;
    
    for my $allele (@alleles) {
        
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

sub feature_type_id {
    my ($self, $feature_type_id) = @_;
    $self->{feature_type_id} = $feature_type_id if defined $feature_type_id;
    # XXX: find the correct feature type id
    return $self->{feature_type_id} || 1;
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
    if (my $tran = $self->{variation_feature}) {
        return $tran->stable_id;
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
 
    if (my $rf_id = $self->{_feature_stable_id} && $type) {
    
        # try to lazy load the feature
        
        if (my $adap = $self->{adaptor}) {
            
            my $get_method = 'get_'.$type.'Adaptor';
            
            if ($adap->db->dnadb->can($get_method)) {
                if (my $fa = $adap->db->dnadb->$get_method) {
                    if (my $f = $fa->fetch_by_stable_id($rf_id)) {
                        $self->{feature} = $f;
                        delete $self->{_feature_stable_id};
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

sub feature_stable_id {
    my $self = shift;
    if (my $f = $self->{feature}) {
        return $f->stable_id;
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
        if (my $adap = $self->{adaptor} || $self->variation_feature->{adaptor}) {
            
            # get the list of possible overlap consequences
            for my $cons (@{ $adap->db->get_OverlapConsequenceAdaptor->fetch_all }) {
                
                # check that this consequence type applies to this feature type
                my $ens_class = $adap->ensembl_class_for_SO_term($cons->feature_SO_term);
                
                if ($ens_class eq ref $self->feature) {
                    
                    # also check if the biotypes match (or if the biotype is not defined)
                    my $biotype = $adap->ensembl_biotype_for_SO_term($cons->feature_SO_term);
                    
                    if (defined $biotype && $self->feature->can('biotype')) {
                        next unless $self->feature->biotype eq $biotype;
                    }
                    
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

sub display_consequence {
    my $self = shift;
    
    my $highest;
    
    for my $allele (@{ $self->alt_alleles }) {
        for my $cons (@{ $allele->consequence_types }) {
            $highest ||= $cons;
            if ($cons->rank < $highest->rank) {
                $highest = $cons;
            }
        }
    }
    
    return $highest->ensembl_term;
}

1;

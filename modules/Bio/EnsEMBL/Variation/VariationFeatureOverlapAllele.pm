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

package Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Constants qw(@OVERLAP_CONSEQUENCES);
use Scalar::Util qw(weaken);

sub new {
    my $class = shift;

    my (
        $variation_feature_overlap, 
        $variation_feature_seq, 
        $is_reference, 
    ) = rearrange([qw(
            VARIATION_FEATURE_OVERLAP 
            VARIATION_FEATURE_SEQ 
            IS_REFERENCE
        )], @_);

    die "VariationFeatureOverlap argument required" unless $variation_feature_overlap;
    die "Allele sequence required" unless $variation_feature_seq;

    my $self = bless {
        variation_feature_overlap   => $variation_feature_overlap,
        variation_feature_seq       => $variation_feature_seq,
        is_reference                => $is_reference,
    }, $class;

    return $self;
}

sub new_fast {
    my ($class, $hashref) = @_;
    my $self = bless $hashref, $class;
    # avoid a memory leak, because the vfo also has a reference to us
    weaken $self->{variation_feature_overlap} if $self->{variation_feature_overlap};
    return $self;
}

sub variation_feature_overlap {
    my ($self, $variation_feature_overlap) = @_;
    
    if ($variation_feature_overlap) {
        $self->{variation_feature_overlap} = $variation_feature_overlap;
        # avoid a memory leak, because the vfo also has a reference to us
        weaken $self->{variation_feature_overlap};
    }
    
    return $self->{variation_feature_overlap};
}

sub variation_feature {
    my $self = shift;
    return $self->variation_feature_overlap->variation_feature;
}

sub feature {
    my $self = shift;
    return $self->variation_feature_overlap->feature;
}

sub feature_seq {
    # the sequence of this allele relative to the feature
    my ($self, $feature_seq) = @_;
    
    $self->{feature_seq} = $feature_seq if $feature_seq;
    
    unless ($self->{feature_seq}) {
        
        # check if we need to reverse complement the vf_seq
        
        my $vf = $self->variation_feature_overlap->variation_feature;
        my $feature = $self->variation_feature_overlap->feature;
        
        if ($vf->strand != $feature->strand) {
            my $vf_seq = $self->{variation_feature_seq};
            reverse_comp(\$vf_seq);
            $self->{feature_seq} = $vf_seq;
        }
        else {
            $self->{feature_seq} = $self->{variation_feature_seq};
        }
    }
    
    return $self->{feature_seq};
}

sub variation_feature_seq {
    # the sequence of this allele relative to the variation feature
    my ($self, $variation_feature_seq) = @_;
    $self->{variation_feature_seq} = $variation_feature_seq if $variation_feature_seq;
    return $self->{variation_feature_seq};
}

sub is_reference {
    my ($self, $is_reference) = @_;
    $self->{is_reference} = $is_reference if defined $is_reference;
    return $self->{is_reference};
}

sub dbID {
    my ($self, $dbID) = @_;
    $self->{dbID} = $dbID if defined $dbID;
    return $self->{dbID};
}

sub allele_string {
    my ($self) = @_;
    
    my $ref = $self->variation_feature_overlap->get_reference_VariationFeatureOverlapllele->variation_feature_seq;
    
    # for the HGMDs and CNV probes where the alleles are artificially set to be
    # the same, just return the reference sequence
    
    if ($ref eq $self->variation_feature_seq) {
        return $ref;
    }
    else {
        return $ref.'/'.$self->variation_feature_seq;
    }
}

sub get_all_OverlapConsequences {
    my ($self, @new_consequences) = @_;
    
    my $cons = $self->{consequence_types};
    
    if (@new_consequences) {
        $cons ||= [];
        push @$cons, @new_consequences;
    }
    
    unless (defined $cons) {
        
        # calculate consequences on the fly
        
        $cons = [];
        
        for my $oc (@OVERLAP_CONSEQUENCES) {
            if ($oc->feature_class eq ref $self->feature) {
                if ($oc->predicate->($self)) {
                    push @$cons, $oc;
                }
            }
        }
    }
    
    $self->{consequence_types} = $cons;
    
    return $cons;
}

sub _is_redundant {
    
    
}

sub SO_isa {
    my ($self, $query) = @_;
    
    if (my $adap = $self->variation_feature_overlap->{adaptor}) {
        if (my $ota = $adap->db->dnadb->get_OntologyTermAdaptor) {
            my $term = $ota->fetch_by_accession();
            my @parents = $ota->fetch_by_child_term($term);
        }
    }
    
    for my $cons (@{ $self->get_all_OverlapConsequences }) {
        if ($cons->SO_term eq $query) {
            return 1;
        }
    } 
}

1;

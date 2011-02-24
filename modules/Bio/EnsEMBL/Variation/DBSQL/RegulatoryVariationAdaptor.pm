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

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::RegulatoryVariationAdaptor;

use Bio::EnsEMBL::Variation::RegulatoryFeatureVariation;
use Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele;
use Bio::EnsEMBL::Variation::MotifFeatureVariation;
use Bio::EnsEMBL::Variation::MotifFeatureVariationAllele;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor);

sub store {
    my ($self, $vfo) = @_;
    
    my $dbh = $self->dbc->db_handle;
    
    my $sth = $dbh->prepare_cached(q{
        INSERT INTO regulatory_variation_allele (
            variation_feature_id,
            feature_stable_id,
            feature_type,
            allele_string,
            is_somatic,
            consequence_types
        ) VALUES (?,?,?,?,?,?)
    });
    
    my $feature_stable_id;
    
    if ($vfo->feature->can('stable_id')) {
        $feature_stable_id = $vfo->feature->stable_id;
    }
    elsif ($vfo->feature->can('_stable_id')) {
        $feature_stable_id = $vfo->feature->_stable_id;
    }
    else {
        die "No stable_id method for feature";
    }
    
    for my $allele (@{ $vfo->alt_alleles }) {
        
        $sth->execute(
            $vfo->variation_feature->dbID,
            $feature_stable_id,
            $self->SO_term_for_ensembl_feature($vfo->feature),
            $allele->allele_string,
            $vfo->variation_feature->is_somatic,
            (join ",", map { $_->SO_term } @{ $allele->consequence_types }),
        );
    }
}

sub fetch_all_by_RegulatoryFeatures {
    my $self = shift;
    return $self->fetch_all_by_Features(@_);
}

sub fetch_all_somatic_by_RegulatoryFeatures {
    my $self = shift;
    return $self->fetch_all_somatic_by_Features(@_);
}

sub fetch_all_by_RegulatoryFeatures_with_constraint {
    my $self = shift;
    return $self->fetch_all_by_Features_with_constraint(@_);
}

sub fetch_all_by_MotifFeatures {
    my $self = shift;
    return $self->fetch_all_by_Features(@_);
}

sub fetch_all_somatic_by_MotifFeatures {
    my $self = shift;
    return $self->fetch_all_somatic_by_Features(@_);
}

sub fetch_all_by_MotifFeatures_with_constraint {
    my $self = shift;
    return $self->fetch_all_by_Features_with_constraint(@_);
}

sub _objs_from_sth {
    my ($self, $sth) = @_;
    
    my (
        $regulatory_variation_allele_id,
        $variation_feature_id, 
        $feature_stable_id, 
        $feature_type,
        $allele_string,
        $consequence_types,
    );
    
    $sth->bind_columns(
        \$regulatory_variation_allele_id,
        \$variation_feature_id, 
        \$feature_stable_id, 
        \$feature_type,
        \$allele_string,
        \$consequence_types,
    );
    
    my %rvs;
    
    while ($sth->fetch) {
        
        my ($ref_allele, $alt_allele)   = split /\//, $allele_string;
       
        # for HGMD mutations etc. just set the alt allele to the ref allele
        $alt_allele ||= $ref_allele;
        
        my $key = $variation_feature_id.'_'.$feature_stable_id.'_'.$feature_type;
        
        my $class = $self->ensembl_variant_class_for_SO_term($feature_type);
            
        my $allele_class = $class.'Allele';
        
        my $rv = $rvs{$key};
        
        unless ($rv) {
            
            $rv = $class->new_fast({
                #dbID                    => $variation_feature_id.','.$feature_stable_id,
                _variation_feature_id   => $variation_feature_id,
                _feature_stable_id      => $feature_stable_id,
                adaptor                 => $self,
            });
            
            $rvs{$key} = $rv;
            
            my $ref_allele = $allele_class->new_fast({
                is_reference                => 1,
                variation_feature_allele    => $ref_allele,
            });
        }
        
        my @cons_types = map { $self->_overlap_consequences->{$_} } 
            split /,/, $consequence_types; # / comment exists to satisfy eclipse!
        
        my $allele = $allele_class->new_fast({
            is_reference                => 0,
            variation_feature_allele    => $alt_allele,
            consequence_types           => \@cons_types,
        });
        
        $rv->alt_alleles($allele);
    }
    
    return [values %rvs];
}

sub _tables {
    return (
        ['regulatory_variation_allele']
    );
}

sub _columns {
    return qw(
        regulatory_variation_allele_id 
        variation_feature_id 
        feature_stable_id 
        feature_type
        allele_string 
        consequence_types
    );
}

1;

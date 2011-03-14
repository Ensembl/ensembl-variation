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

package Bio::EnsEMBL::Variation::DBSQL::RegulationVariationAdaptor;

use Bio::EnsEMBL::Variation::RegulatoryFeatureVariation;
use Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele;
use Bio::EnsEMBL::Variation::MotifFeatureVariation;
use Bio::EnsEMBL::Variation::MotifFeatureVariationAllele;
use Bio::EnsEMBL::Variation::ExternalFeatureVariation;
use Bio::EnsEMBL::Variation::ExternalFeatureVariationAllele;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor);

sub store {
    my ($self, $vfo) = @_;
    
    my $dbh = $self->dbc->db_handle;
    
    my $sth = $dbh->prepare_cached(q{
        INSERT INTO regulation_variation (
            variation_feature_id,
            feature_stable_id,
            feature_interdb_stable_id,
            feature_type_attrib_id,
            allele_string,
            somatic,
            consequence_types,
            target_feature_stable_id
        ) VALUES (?,?,?,?,?,?,?,?)
    });
    
    for my $allele (@{ $vfo->alt_alleles }) {
        
        my ($feat_stable_id, $interdb_stable_id);
        
        if ($vfo->feature->can('stable_id')) {
            $feat_stable_id = $vfo->feature->stable_id;
        }
        elsif ($vfo->feature->can('interdb_stable_id')) {
            $interdb_stable_id = $vfo->feature->interdb_stable_id;
        }
        else {
            die "Can't store vfo, no identifier for feature!";
        }
        
        my $feat_class_attrib_id = $self->AttributeAdaptor->attrib_id_for_type_value(
            'ens_feature_class', 
            ref $allele->feature
        );
       
        die "No feature class for ".(ref $allele->feature) unless $feat_class_attrib_id;

        $sth->execute(
            $vfo->variation_feature_id,
            $feat_stable_id,
            $interdb_stable_id,
            $feat_class_attrib_id,
            $allele->allele_string,
            $vfo->variation_feature->is_somatic,
            (join ',', map { $_->SO_term } @{ $allele->consequence_types }),
            $vfo->target_feature_stable_id,
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
    return $self->fetch_all_by_Features_no_stable_id(@_);
}

sub fetch_all_somatic_by_MotifFeatures {
    my $self = shift;
    return $self->fetch_all_somatic_by_Features_no_stable_id(@_);
}

sub fetch_all_by_MotifFeatures_with_constraint {
    my $self = shift;
    return $self->fetch_all_by_Features_no_stable_id_with_constraint(@_);
}

sub fetch_all_by_ExternalFeatures {
    my $self = shift;
    return $self->fetch_all_by_Features_no_stable_id(@_);
}

sub fetch_all_somatic_by_ExternalFeatures {
    my $self = shift;
    return $self->fetch_all_somatic_by_Features_no_stable_id(@_);
}

sub fetch_all_by_ExternalFeatures_with_constraint {
    my $self = shift;
    return $self->fetch_all_by_Features_no_stable_id_with_constraint(@_);
}

sub fetch_all_by_Features_no_stable_id {
    my ($self, $features) = @_;
    return $self->fetch_all_by_Features_no_stable_id_with_constraint($features, 'somatic = 0');
}

sub fetch_all_somatic_by_Features_no_stable_id {
    my ($self, $features) = @_;
    return $self->fetch_all_by_Features_no_stable_id_with_constraint($features, 'somatic = 1');
}

sub fetch_all_by_Features_no_stable_id_with_constraint {
    my ($self, $features, $constraint) = @_;
    
    my $dbh = $self->dbc->db_handle;
   
    my %feats_by_id = map { $_->interdb_stable_id => $_ } @$features;
    
    my $id_str = join',', map {"'$_'"} keys %feats_by_id;
    
    my $full_constraint = "feature_interdb_stable_id in ( $id_str )";
    $full_constraint .= " AND $constraint" if $constraint;
    
    my $vfos = $self->generic_fetch($full_constraint);
    
    for my $vfo (@$vfos) {
        if ($vfo->{_feature_interdb_stable_id}) {
            my $feat_id = delete $vfo->{_feature_interdb_stable_id};
            $vfo->{feature} = $feats_by_id{$feat_id};
        }
    }
    
    return $vfos;
}

sub _objs_from_sth {
    my ($self, $sth) = @_;
    
    my (
        $funcgen_variation_allele_id,
        $variation_feature_id, 
        $feature_stable_id,
        $feature_interdb_stable_id,
        $feature_type_attrib_id,
        $allele_string,
        $consequence_types,
        $target_feature_stable_id,
    );
    
    $sth->bind_columns(
        \$funcgen_variation_allele_id,
        \$variation_feature_id, 
        \$feature_stable_id,
        \$feature_interdb_stable_id,
        \$feature_type_attrib_id,
        \$allele_string,
        \$consequence_types,
        \$target_feature_stable_id,
    );
    
    my %fgvs;
    
    while ($sth->fetch) {
        
        my ($ref_allele, $alt_allele)   = split /\//, $allele_string;
       
        # for HGMD mutations etc. just set the alt allele to the ref allele
        $alt_allele ||= $ref_allele;
        
        my $key = $variation_feature_id.'_'.$feature_type_attrib_id.'_';
        
        $key .= $feature_stable_id ? $feature_stable_id : $feature_interdb_stable_id;
        
        my $feature_class = $self->AttributeAdaptor->attrib_value_for_id($feature_type_attrib_id);
        
        my $var_class = $self->AttributeAdaptor->ensembl_variant_class_for_feature_class($feature_class);
            
        my $allele_class = $var_class.'Allele';
        
        my $fgv = $fgvs{$key};
        
        unless ($fgv) {
            
            $fgv = $var_class->new_fast({
                _variation_feature_id       => $variation_feature_id,
                _feature_stable_id          => $feature_stable_id,
                _feature_interdb_stable_id  => $feature_interdb_stable_id,
                _target_feature_stable_id   => $target_feature_stable_id,
                adaptor                     => $self,
            });
            
            $fgvs{$key} = $fgv;
            
            my $reference_allele = $allele_class->new_fast({
                is_reference                => 1,
                variation_feature_allele    => $ref_allele,
            });
        }
        
        my @cons_types = map { $self->AttributeAdaptor->OverlapConsequence_for_SO_term($_) } 
            split /,/, $consequence_types; # / comment exists to satisfy eclipse!
        
        my $allele = $allele_class->new_fast({
            is_reference                => 0,
            variation_feature_allele    => $alt_allele,
            consequence_types           => \@cons_types,
        });
        
        $fgv->alt_alleles($allele);
    }
    
    return [values %fgvs];
}

sub _tables {
    return (
        ['regulation_variation']
    );
}

sub _columns {
    return qw(
        regulation_variation_id 
        variation_feature_id 
        feature_stable_id 
        feature_interdb_stable_id
        feature_type_attrib_id
        allele_string 
        consequence_types
        target_feature_stable_id
    );
}

1;

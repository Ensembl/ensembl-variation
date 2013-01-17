=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Variation::DBSQL::RegulatoryFeatureVariationAdaptor

=head1 SYNOPSIS
    my $reg = 'Bio::EnsEMBL::Registry';
  
    $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
    my $rfa  = $reg->get_adaptor('human', 'funcgen', 'RegulatoryFeature');
    my $va   = $reg->get_adaptor('human', 'variation', 'Variation');
    my $vfa  = $reg->get_adaptor('human', 'variation', 'VariationFeature');
    my $rfva = $reg->get_adaptor('human','variation','RegulatoryFeatureVariation');

    # fetch RegulatoryFeature by dbID
    my $rf_stable_id = 'ENSR00000316845';
    my $rf           = $rfa->fetch_by_stable_id($rf_stable_id);

    # fetch all RegulatoryFeatureVariations falling in the MotifFeature
    my $rfvs = $rfva->fetch_all_by_RegulatoryFeatures([$rf]);

    # fetch by VariationFeatures
    my $variation_name = 'rs191666497';
    my $v   = $va->fetch_by_name($variation_name);
    my $vfs = $vfa->fetch_all_by_Variation($v);
    
    # fetch all RegulatoryFeatureVariations for the list of VariationFeatures
    for my $vf (@$vfs) {
        print $vf->variation_name, "\n";
    }
    $rfvs = $rfva->fetch_all_by_VariationFeatures($vfs);
    for my $rfv (@$rfvs) {
        print $rfv->regulatory_feature_stable_id, "\n";
    }
 
=head1 DESCRIPTION

This adaptor allows you to fetch RegulatoryFeatureVariation objects either by the RegulatoryFeature
the associated VariationFeature falls in, or by VariationFeature directly. Storing
RegulatoryFeatureVariation objects in a variation schema database is also supported. In the
database there will a separate row for each alternative allele of a RegulatoryFeatureVariation, 
but the methods here will fetch all alleles associated with the RegulatoryFeatureVariation
at once.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::RegulatoryFeatureVariationAdaptor;

use Bio::EnsEMBL::Variation::RegulatoryFeatureVariation;
use Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);

use base qw(Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor);

=head2 store

  Arg [1]    : Bio::EnsEMBL::Variation::RegulatoryFeatureVariation $rfv
  Description: Store the RegulatoryFeatureVariation in the database
  Status     : At risk

=cut

sub store {
    my ($self, $rfv) = @_;
    my $dbh = $self->dbc->db_handle;
    my $sth = $dbh->prepare_cached(q{
        INSERT DELAYED INTO regulatory_feature_variation (
            variation_feature_id,
            feature_stable_id,
            feature_type,
            allele_string,
            somatic,
            consequence_types
        ) VALUES (?,?,?,?,?,?)
    });

    for my $allele (@{ $rfv->get_all_alternate_RegulatoryFeatureVariationAlleles }) {
        $sth->execute(
            $rfv->variation_feature->dbID,
            $rfv->feature->stable_id,
            $rfv->regulatory_feature->feature_type->name,
            $allele->allele_string,
            $rfv->variation_feature->is_somatic,
            (join ',', map { $_->SO_term } @{ $allele->get_all_OverlapConsequences }),
        );
    }
}

=head2 fetch_all_by_RegulatoryFeatures

  Arg [1]    : listref of Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Description: Fetch all germline RegulatoryFeatureVariations associated with the
               given list of RegulatoryFeatures
  Returntype : listref of Bio::EnsEMBL::Variation::RegulatoryFeatureVariation
  Status     : Stable

=cut

sub fetch_all_by_RegulatoryFeatures {
    my ($self, $regulatory_features) = @_;
    return $self->fetch_all_by_RegulatoryFeatures_with_constraint($regulatory_features, 'somatic = 0');
}

=head2 fetch_all_somatic_by_RegulatoryFeatures

  Arg [1]    : listref of Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Description: Fetch all somatic RegulatoryFeatureVariations associated with the
               given list of RegulatoryFeatures
  Returntype : listref of Bio::EnsEMBL::Variation::RegulatoryFeatureVariation
  Status     : Stable

=cut

sub fetch_all_somatic_by_RegulatoryFeatures {
    my ($self, $regulatory_features) = @_;
    return $self->fetch_all_by_RegulatoryFeatures_with_constraint($regulatory_features, 'somatic = 1');
}

=head2 fetch_all_by_RegulatoryFeatures_SO_terms

  Arg [1]    : listref of Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Arg [2]    : listref of SO terms
  Description: Fetch all germline RegulatoryFeatureVariations associated with the
               given list of RegulatoryFeatures and with consequences from given list of SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::RegulatoryFeatureVariation
  Status     : At risk

=cut

sub fetch_all_by_RegulatoryFeatures_SO_terms {
    my ($self, $regulatory_features, $terms) = @_;
    my $constraint = $self->_get_consequence_constraint($terms);
    return $self->fetch_all_by_RegulatoryFeatures_with_constraint($regulatory_features, $constraint.' AND somatic = 0');
}

=head2 fetch_all_somatic_by_RegulatoryFeatures_SO_terms

  Arg [1]    : listref of Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Arg [2]    : listref of SO terms
  Description: Fetch all somatic RegulatoryFeatureVariations associated with the
               given list of RegulatoryFeatures and with consequences from given list of SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::RegulatoryFeatureVariation
  Status     : At risk

=cut

sub fetch_all_somatic_by_RegulatoryFeatures_SO_terms {
    my ($self, $regulatory_features, $terms) = @_;
    my $constraint = $self->_get_consequence_constraint($terms);
    return $self->fetch_all_by_RegulatoryFeatures_with_constraint($regulatory_features, $constraint.' AND somatic = 1');
}

=head2 fetch_all_by_VariationFeatures_SO_terms

  Arg [1]    : listref of Bio::EnsEMBL::Variation::VariationFeature
  Arg [2]    : listref of SO terms
  Description: Fetch all germline RegulatoryFeatureVariations associated with the
               given list of VariationFeatures and with consequences from given list of SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::RegulatoryFeatureVariation
  Status     : At risk

=cut

sub fetch_all_by_VariationFeatures_SO_terms {
    my ($self, $vfs, $regulatory_features, $terms, $without_children, $included_so) = @_;
    my $constraint = $self->_get_consequence_constraint($terms, $without_children, $included_so);
    if (!$constraint) {
      return [];
    }
    return $self->SUPER::fetch_all_by_VariationFeatures_with_constraint($vfs, $regulatory_features, $constraint);
}

=head2 count_all_by_VariationFeatures_SO_terms

  Arg [1]    : listref of Bio::EnsEMBL::Variation::VariationFeatures
  Arg [2]    : listref of SO terms
  Description: Count RegulatoryFeatureVariations associated with given
               VariationFeatures and with consequences from given list of SO terms
  Returntype : int
  Status     : At risk

=cut

sub count_all_by_VariationFeatures_SO_terms {
    my ($self, $vfs, $regulatory_features, $terms, $included_so) = @_;
    my $constraint = $self->_get_consequence_constraint($terms, 1, $included_so);
    if (!$constraint) {
      return 0;
    }
    return $self->SUPER::count_all_by_VariationFeatures_with_constraint($vfs, $regulatory_features, $constraint);
}

sub fetch_all_by_RegulatoryFeatures_with_constraint {
    my ($self, $regulatory_features, $constraint) = @_;
    return $self->SUPER::fetch_all_by_Features_with_constraint($regulatory_features, $constraint);
}

sub _objs_from_sth {
    my ($self, $sth) = @_;
    #warn $sth->sql;
    my (
        $regulatory_feature_variation_id,
        $variation_feature_id, 
        $feature_stable_id, 
        $feature_type,
        $allele_string,
        $somatic,
        $consequence_types,
    );
    
    $sth->bind_columns(
        \$regulatory_feature_variation_id,
        \$variation_feature_id, 
        \$feature_stable_id, 
        \$feature_type,
        \$allele_string,
        \$somatic,
        \$consequence_types,
    );
    
    my %rfvs;
    
    while ($sth->fetch) {
        my ($ref_allele, $alt_allele)   = split /\//, $allele_string;
        my $key = $variation_feature_id.'_'.$feature_stable_id;
        my $rfv = $rfvs{$key};
        unless ($rfv) {
            $rfv = Bio::EnsEMBL::Variation::RegulatoryFeatureVariation->new_fast({
                _variation_feature_id => $variation_feature_id,
                _feature_stable_id    => $feature_stable_id,
                feature_type          => $feature_type,
                adaptor               => $self,
            });
            
            $rfvs{$key} = $rfv;
            
            my $ref_allele = Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele->new_fast({
                is_reference                 => 1,
                variation_feature_seq        => $ref_allele,
                regulatory_feature_variation => $rfv,
                dbID                         => $regulatory_feature_variation_id,
            });

            $rfv->add_RegulatoryFeatureVariationAllele($ref_allele);
        }

        my $overlap_consequences = [ map { $OVERLAP_CONSEQUENCES{$_} } split /,/, $consequence_types ];
        my $allele = Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele->new_fast({
            is_reference                 => 0,
            variation_feature_seq        => $alt_allele,
            regulatory_feature_variation => $rfv, 
            overlap_consequences         => $overlap_consequences,
            dbID                         => $regulatory_feature_variation_id,
        });
        $rfv->add_RegulatoryFeatureVariationAllele($allele);
    }
    return [values %rfvs];
}

sub _tables {
    return (
        ['regulatory_feature_variation', 'rfv']
    );
}

sub _columns {
    return 
    qw(
        regulatory_feature_variation_id 
        variation_feature_id 
        feature_stable_id 
        feature_type
        allele_string
        somatic 
        consequence_types 
    );
}

1;

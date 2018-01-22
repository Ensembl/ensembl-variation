=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

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
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

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
  $sth->finish();
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

sub _fetch_all_by_VariationFeatures_no_DB {
  my ($self, $vfs, $features) = @_;
  
  # get features?
  if(!$features || !@$features) {
    my $slices = $self->_get_ranged_slices_from_VariationFeatures($vfs);

    if(
      my $fg_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
        -species  => $self->db->species, 
        -type     => 'RegulatoryFeature',
      )
    ) {
      @$features = map {@{$fg_adaptor->fetch_all_by_Slice($_)}} @$slices;
    }
  }
  
  my @return;
  
  foreach my $f(@{$features || []}) {

    my $f_slice = $f->slice;
    
    foreach my $vf(grep {overlap($f->{start}, $f->{end}, $_->{start}, $_->{end})} @$vfs) {
      my $vfo = Bio::EnsEMBL::Variation::RegulatoryFeatureVariation->new(
        -variation_feature => $vf,
        -feature           => $f,
        -adaptor           => $self,
        -no_ref_check      => 1,
        -no_transfer       => ($vf->slice + 0) == ($f_slice + 0)
      );
      
      $vf->add_RegulatoryFeatureVariation($vfo);
      push @return, $vfo;
    }
  }
  
  return \@return;
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

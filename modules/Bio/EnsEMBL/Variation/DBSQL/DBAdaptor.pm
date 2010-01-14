#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::DBAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::DBAdaptor

=head1 SYNOPSIS



=head1 DESCRIPTION

This module provides a connection to an Ensembl variation database and
provides a means to obtain ObjectAdaptors.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;


package Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;


use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our @ISA = ('Bio::EnsEMBL::DBSQL::DBAdaptor');

sub get_available_adaptors{
    my %pairs = (
		 'Population' => 'Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor',
		 'Individual' => 'Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor',
		 'Variation' => 'Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor',
		 'VariationFeature' => 'Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor',
		 'StructuralVariation' => 'Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor',
		 'VariationAnnotation' => 'Bio::EnsEMBL::Variation::DBSQL::VariationAnnotationAdaptor',
		 'AlleleFeature' => 'Bio::EnsEMBL::Variation::DBSQL::AlleleFeatureAdaptor',
		 'LDFeatureContainer' => 'Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor',
#    		 'IndividualGenotype' => 'Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor',
		 'IndividualGenotype' => 'Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor',
		 'PopulationGenotype' => 'Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor',
		 'TranscriptVariation' => 'Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor',
		 'VariationGroup' => 'Bio::EnsEMBL::Variation::DBSQL::VariationGroupAdaptor',
		 'AlleleGroup' => 'Bio::EnsEMBL::Variation::DBSQL::AlleleGroupAdaptor',
		 'VariationGroupFeature' => 'Bio::EnsEMBL::Variation::DBSQL::VariationGroupFeatureAdaptor',
		 'MetaCoordContainer' => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
		 'MetaContainer' => 'Bio::EnsEMBL::Variation::DBSQL::MetaContainer',
		 'ReadCoverage' => 'Bio::EnsEMBL::Variation::DBSQL::ReadCoverageAdaptor',
#		 'CompressedGenotype' => 'Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor',
		 #'VariationFeatureCollection' => 'Bio::EnsEMBL::Variation::Collection::VariationFeature',
		 'VariationAnnotation' => 'Bio::EnsEMBL::Variation::DBSQL::VariationAnnotationAdaptor',
		 'ReadCoverageCollection' => 'Bio::EnsEMBL::Variation::DBSQL::ReadCoverageCollectionAdaptor',
                 'VariationSet' => 'Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor',
                 );
    return (\%pairs);
}

1;

=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

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

=head1 METHODS

=cut

use strict;
use warnings;


package Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;


use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our @ISA = ('Bio::EnsEMBL::DBSQL::DBAdaptor');

our $DEFAULT_INCLUDE_FAILED_VARIATIONS = 0;
our $DEFAULT_INCLUDE_NON_SIGNIFICANT_PHENOTYPES = 0;

sub get_available_adaptors{
    my %pairs = (
		'Population'                      => 'Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor',
		'Individual'                      => 'Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor',
		'Variation'                       => 'Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor',
		'VariationFeature'                => 'Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor',
		'StructuralVariation'             => 'Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor',
		'SupportingStructuralVariation'   => 'Bio::EnsEMBL::Variation::DBSQL::SupportingStructuralVariationAdaptor',
		'StructuralVariationFeature'      => 'Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor',
		'Study'                           => 'Bio::EnsEMBL::Variation::DBSQL::StudyAdaptor',
		'AlleleFeature'                   => 'Bio::EnsEMBL::Variation::DBSQL::AlleleFeatureAdaptor',
		'LDFeatureContainer'              => 'Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor',
		'IndividualGenotype'              => 'Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor',
		'IndividualGenotypeFeature'       => 'Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeFeatureAdaptor',
		'PopulationGenotype'              => 'Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor',
		'TranscriptVariation'             => 'Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor',
		'MetaCoordContainer'              => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
		'MetaContainer'                   => 'Bio::EnsEMBL::Variation::DBSQL::MetaContainer',
		'ReadCoverage'                    => 'Bio::EnsEMBL::Variation::DBSQL::ReadCoverageAdaptor',
		'GenotypeCode'                    => 'Bio::EnsEMBL::Variation::DBSQL::GenotypeCodeAdaptor',
		'VariationAnnotation'             => 'Bio::EnsEMBL::Variation::DBSQL::VariationAnnotationAdaptor',
		'ReadCoverageCollection'          => 'Bio::EnsEMBL::Variation::DBSQL::ReadCoverageCollectionAdaptor',
		'VariationSet'                    => 'Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor',
		'OverlapConsequence'              => 'Bio::EnsEMBL::Variation::DBSQL::OverlapConsequenceAdaptor',
		'Attribute'                       => 'Bio::EnsEMBL::Variation::DBSQL::AttributeAdaptor',
		'Allele'                          => 'Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor',
		'ProteinFunctionPredictionMatrix' => 'Bio::EnsEMBL::Variation::DBSQL::ProteinFunctionPredictionMatrixAdaptor',
    'Phenotype'                       => 'Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor',
    'RegulatoryFeatureVariation'      => 'Bio::EnsEMBL::Variation::DBSQL::RegulatoryFeatureVariationAdaptor',
    'MotifFeatureVariation'           => 'Bio::EnsEMBL::Variation::DBSQL::MotifFeatureVariationAdaptor',
    'PhenotypeFeature'                => 'Bio::EnsEMBL::Variation::DBSQL::PhenotypeFeatureAdaptor',
    'Publication'                     => 'Bio::EnsEMBL::Variation::DBSQL::PublicationAdaptor',
    'StructuralVariationSample'       => 'Bio::EnsEMBL::Variation::DBSQL::StructuralVariationSampleAdaptor',
    );
	
    return (\%pairs);
}


=head2 include_failed_variations

  Arg [1]    : int $newval (optional)
  Example    :
		#ÊGet a DBAdaptor for the human variation database
		my $dba = $registry->get_DBAdaptor('human','variation');
		
		#ÊConfigure the DBAdaptor to return failed variations when using
		#Êfetch methods in the various object adaptors
		$dba->include_failed_variations(1);
		
		#ÊGet a variation set adaptor
		my $vs_adaptor = $dba->get_VariationSetAdaptor();
		
		#ÊGet a variation set for the 1000 genomes high coverage Yoruba trio data
		my $vs = $vs_adaptor->fetch_by_name('1000 genomes - High coverage - Trios - YRI');
		
		# Get the iterator for the variations belonging to this variation set.
		#ÊThis will now include variations that has been flagged as being failed.
		#ÊThe default behaviour is not to return these.
		my $it = $vs->get_Variation_Iterator();
		
		# Iterate over the variations
		while ($it->has_next()) {
		
		    # Get the next variation object in the iterator
		    my $v = $it->next();
		    
		    # Check if the variation is flagged as failed
		    if ($v->is_failed()) {
			# Do something...
		    }
		    # If not, check if any of its subsnps have been flagged as failed
		    elsif ($v->has_failed_subsnps()) {
			#ÊDo something else...
		    }
		    else {
			#ÊDo something else...
		    }
		}
		
  Description: Getter/Setter for the behaviour of the adaptors connected through this
	       DBAdaptor when it comes to variations and alleles that have been flagged as failed.
	       The default behaviour is not to return these variations or alleles in e.g. the
	       'fetch_all_by...'-type methods. If this flag is set, those methods will
	       instead also return failed variations and alleles. Note that a variation is considered
	       failed when the variation itself is failed. If only some alleles belonging
	       to the variation are failed, the entire variation will not be considered
	       to be failed.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub include_failed_variations {
    my $self = shift;
    my $include = shift;
    
    #ÊIf the flag should be modified, do that
    if (defined($include)) {$self->{'include_failed_variations'} = $include;}
    
    #ÊIn case the flag has not been set at all, set it to the default value
    unless (exists($self->{'include_failed_variations'})) {$self->{'include_failed_variations'} = $DEFAULT_INCLUDE_FAILED_VARIATIONS;}
    
    return $self->{'include_failed_variations'};
}


# API-internal method for getting the constraint to filter out failed variations. Assumes that the
# failed_variation table has been (left) joined to the query and that the table alias is either supplied
# or equals 'fv'
sub _exclude_failed_variations_constraint {
    my $self = shift;
    my $table_alias = shift;
    
    # If not specified, assume that the failed_variation table alias is 'fv'
    $table_alias ||= 'fv';
    
    return $self->_exclude_failed_constraint('variation_id',$table_alias);
}

# API-internal method for getting the constraint to filter out failed structural variations. Assumes that the
# failed_structural_variation table has been (left) joined to the query and that the table alias is either supplied
# or equals 'fsv'
sub _exclude_failed_structural_variations_constraint {
    my $self = shift;
    my $table_alias = shift;
    
    # If not specified, assume that the failed_structural_variation table alias is 'fsv'
    $table_alias ||= 'fsv';
    
    return $self->_exclude_failed_constraint('structural_variation_id',$table_alias);
}

# API-internal method for getting the constraint to filter out failed alleles. Assumes that the
# failed_allele table has been (left) joined to the query and that the table alias is either supplied
# or equals 'fa'
sub _exclude_failed_alleles_constraint {
    my $self = shift;
    my $table_alias = shift;
    
    # If not specified, assume that the failed_variation table alias is 'fv'
    $table_alias ||= 'fa';
    
    return $self->_exclude_failed_constraint('allele_id',$table_alias);
}

sub _exclude_failed_constraint {
    my $self = shift;
    my $key_column = shift;
    my $table_alias = shift;
    
    #ÊIf we should include failed objects, no extra condition is needed
    return qq{ 1 } if ($self->include_failed_variations());
    
    # Otherwise, add a constraint on the alias table to have the key_column NULL
    my $constraint = qq{
	(
	    $table_alias.$key_column IS NULL
	)
    };
    
    return $constraint;
}


=head2 include_non_significant_phenotype_associations

  Arg [1]    : int $newval (optional)
  Example    :
    # Get a DBAdaptor for the human variation database
    my $dba = $registry->get_DBAdaptor('human','phenotypefeature');
    
    # Configure the DBAdaptor to return non significant phenotype associations when using
    # fetch methods in the various object adaptors
    $dba->include_non_significant_phenotype_associations(1);
    
  Description: Getter/Setter for the behaviour of the adaptors connected through this
         DBAdaptor when it comes to phenotype feature.
         The default behaviour is to return the phenotype features with significance results only e.g. the
         'fetch_all_by...'-type methods. If this flag is set, those methods will
         instead also return phenotype features with non significant results. 
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub include_non_significant_phenotype_associations {
    my $self = shift;
    my $include = shift;
    
    # If the flag should be modified, do that
    if (defined($include)) {$self->{'include_non_significant_phenotypes'} = $include;}
    
    # In case the flag has not been set at all, set it to the default value
    unless (exists($self->{'include_non_significant_phenotypes'})) {$self->{'include_non_significant_phenotypes'} = $DEFAULT_INCLUDE_NON_SIGNIFICANT_PHENOTYPES;}
    
    return $self->{'include_non_significant_phenotypes'};
}

1;

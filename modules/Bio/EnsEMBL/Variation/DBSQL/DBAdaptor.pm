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

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::DBAdaptor
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
our $DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME  = 1;

sub get_available_adaptors{
    my %pairs = (
      'Allele'                                 => 'Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor',
      'AlleleFeature'                          => 'Bio::EnsEMBL::Variation::DBSQL::AlleleFeatureAdaptor',
      'Attribute'                              => 'Bio::EnsEMBL::Variation::DBSQL::AttributeAdaptor',
      'GenotypeCode'                           => 'Bio::EnsEMBL::Variation::DBSQL::GenotypeCodeAdaptor',
      'Individual'                             => 'Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor',
      'IndividualGenotype'                     => 'Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor',
      'IndividualGenotypeFeature'              => 'Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeFeatureAdaptor',
      'LDFeatureContainer'                     => 'Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor',
      'MetaCoordContainer'                     => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
      'MetaContainer'                          => 'Bio::EnsEMBL::Variation::DBSQL::MetaContainer',
      'MotifFeatureVariation'                  => 'Bio::EnsEMBL::Variation::DBSQL::MotifFeatureVariationAdaptor',
      'Phenotype'                              => 'Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor',
      'PhenotypeFeature'                       => 'Bio::EnsEMBL::Variation::DBSQL::PhenotypeFeatureAdaptor',
      'Population'                             => 'Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor',
      'PopulationGenotype'                     => 'Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor',
      'ProteinFunctionPredictionMatrix'        => 'Bio::EnsEMBL::Variation::DBSQL::ProteinFunctionPredictionMatrixAdaptor',
      'Publication'                            => 'Bio::EnsEMBL::Variation::DBSQL::PublicationAdaptor',
      'ReadCoverage'                           => 'Bio::EnsEMBL::Variation::DBSQL::ReadCoverageAdaptor',
      'RegulatoryFeatureVariation'             => 'Bio::EnsEMBL::Variation::DBSQL::RegulatoryFeatureVariationAdaptor',
      'Sample'                                 => 'Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor',
      'SampleGenotype'                         => 'Bio::EnsEMBL::Variation::DBSQL::SampleGenotypeAdaptor',
      'SampleGenotypeFeature'                  => 'Bio::EnsEMBL::Variation::DBSQL::SampleGenotypeFeatureAdaptor',
      'Source'                                 => 'Bio::EnsEMBL::Variation::DBSQL::SourceAdaptor',
      'StrainSlice'                            => 'Bio::EnsEMBL::Variation::DBSQL::StrainSliceAdaptor',
      'StructuralVariation'                    => 'Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor',
      'StructuralVariationFeature'             => 'Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor',
      'StructuralVariationPopulationFrequency' => 'Bio::EnsEMBL::Variation::DBSQL::StructuralVariationPopulationFrequencyAdaptor',
      'StructuralVariationSample'              => 'Bio::EnsEMBL::Variation::DBSQL::StructuralVariationSampleAdaptor',
      'Study'                                  => 'Bio::EnsEMBL::Variation::DBSQL::StudyAdaptor',
      'SupportingStructuralVariation'          => 'Bio::EnsEMBL::Variation::DBSQL::SupportingStructuralVariationAdaptor',
      'TranscriptVariation'                    => 'Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor',
      'Variation'                              => 'Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor',
      'VariationFeature'                       => 'Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor',
      'VariationSet'                           => 'Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor',
      'VCFCollection'                          => 'Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor',
      'TranscriptHaplotype'                    => 'Bio::EnsEMBL::Variation::DBSQL::TranscriptHaplotypeAdaptor',
    );
	
    return (\%pairs);
}


=head2 include_failed_variations

  Arg [1]    : int $newval (optional)
  Example    :
		#Get a DBAdaptor for the human variation database
		my $dba = $registry->get_DBAdaptor('human','variation');
		
		#Configure the DBAdaptor to return failed variations when using
		#fetch methods in the various object adaptors
		$dba->include_failed_variations(1);
		
		#Get a variation set adaptor
		my $vs_adaptor = $dba->get_VariationSetAdaptor();
		
		#Get a variation set for the 1000 genomes high coverage Yoruba trio data
		my $vs = $vs_adaptor->fetch_by_name('1000 genomes - High coverage - Trios - YRI');
		
		# Get the iterator for the variations belonging to this variation set.
		#This will now include variations that has been flagged as being failed.
		#The default behaviour is not to return these.
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
		    elsif ($v->has_failed_alleles()) {
			#Do something else...
		    }
		    else {
			#Do something else...
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
    
    #If the flag should be modified, do that
    if (defined($include)) {$self->{'include_failed_variations'} = $include;}
    
    #In case the flag has not been set at all, set it to the default value
    unless (exists($self->{'include_failed_variations'})) {$self->{'include_failed_variations'} = $DEFAULT_INCLUDE_FAILED_VARIATIONS;}
    
    return $self->{'include_failed_variations'};
}


# API-internal method for getting the constraint to filter out variations not shown by default 
# (those failing standard QC which are not cited). 
# Assumes that the variation table has been (left) joined to the query and that the table 
# alias is either supplied or equals 'v'
sub _exclude_failed_variations_constraint {
    my $self = shift;
    my $table_alias = shift;
    
    # If not specified, assume that the variation table alias is 'v'
    $table_alias ||= 'v';
    
    return $table_alias . ".display =1";
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
    
    #If we should include failed objects, no extra condition is needed
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


=head2 shift_hgvs_variants_3prime

  Arg [1]    : int $newval (optional)
  Example    :
		#Get a DBAdaptor for the human variation database
		my $dba = $registry->get_DBAdaptor('human','variation');
		
		#Configure the DBAdaptor to return failed variations when using
		#fetch methods in the various object adaptors
		$dba->shift_hgvs_variants_3prime(1);
		
		#Proceed to extract HGVS annotation as normal

		
  Description: Getter/Setter for the behaviour of the adaptors connected through this
	       DBAdaptor when it comes to HGVS transcript and protein level annotation.
	       The default behaviour is to keep the positions as input. Formal HGVS annotation
               requires any variant to be described at the most 3 prime location possible which 
               changes the variation consequence in a very small number of cases.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub shift_hgvs_variants_3prime{
    my $self = shift;
    my $move_3prime = shift;

    #If the flag should be modified, do that
    if (defined($move_3prime)) {$self->{'shift_hgvs_variants_3prime'} = $move_3prime;}
    
    #In case the flag has not been set at all, set it to the default value
    unless (exists($self->{'shift_hgvs_variants_3prime'})) {$self->{'shift_hgvs_variants_3prime'} = $DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME;}
    
    return $self->{'shift_hgvs_variants_3prime'};
}

=head2 use_vcf

  Arg [1]    : int $newval (optional)
  Example    :
		# Get a LD feature container adaptor
		my $ldfca = $registry->get_DBAdaptor('human','variation','ldfeaturecontainer');
		
		# tell the adaptor to include genotypes from VCF files
		$ldfca->db->use_vcf(1);
		
		# fetch LD container as normal
    my $ldfc = $ldfca->fetch_by_Slice($slice, $population);

  Description: Getter/Setter the use_vcf property. Switching this to a
               non-zero value instructs the API to look in configured
               VCF files for genotype-based data. Setting to 1 tells the
               API to use both the DB and VCF; setting to 2 tells the API
               to use only VCF.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub use_vcf {
  my $self = shift;
  
  my ($new_use, $vcf_config_file) = @_;
  
  # default to 0
  $self->{use_vcf} = 0 if !exists($self->{use_vcf});
  
  # allow user to switch
  $self->{use_vcf} = $new_use if defined($new_use);
  
  $self->vcf_config_file($vcf_config_file) if defined($vcf_config_file);
  
  return $self->{use_vcf}
}

sub vcf_config {
  my $self = shift;

  if(!exists($self->{vcf_config})) {
    if(@_) {
      $self->{vcf_config} = shift;
    }
    else {
      $self->{vcf_config} = {};
    }
  }
  
  return $self->{vcf_config};
}

sub vcf_config_file {
  my $self = shift;
  
  $self->{vcf_config_file} = shift if @_;
  
  return $self->{vcf_config_file};
}

sub vcf_root_dir {
  my $self = shift;
  
  $self->{vcf_root_dir} = shift if @_;
  
  return $self->{vcf_root_dir};
}

sub vcf_tmp_dir {
  my $self = shift;
  
  $self->{vcf_tmp_dir} = shift if @_;
  
  return $self->{vcf_tmp_dir};
}

1;

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor

=head1 SYNOPSIS

  $vdb = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);
  $db  = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  $sa  = $db->get_SliceAdaptor();
  $lda = $vdb->get_LDFeatureContainerAdaptor();
  $vfa =  $vdb->get_VariationFeatureAdaptor();

  # Get a LDFeatureContainer in a region
  $slice = $sa->fetch_by_region('chromosome', 'X', 1e6, 2e6);

  $ldContainer = $lda->fetch_by_Slice($slice);

  print "Name of the ldContainer is: ", $ldContainer->name();

  # fetch ld featureContainer for a particular variation feature

  $vf = $vfa->fetch_by_dbID(145);

  $ldContainer = $lda->fetch_by_VariationFeature($vf);

  print "Name of the ldContainer: ", $ldContainer->name();


=head1 DESCRIPTION

This adaptor provides database connectivity for LDFeature objects.
LD Features may be retrieved from the Ensembl variation database by
several means using this module.

=head1 AUTHOR - Daniel Rios

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Variation::LDFeatureContainer;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');

=head2 fetch_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch genes on. Assuming it is always correct (in the top level)
  Arg [2]    : (optional) int $population_id. Population where we want to select the LD information
  Example    : $ldFeatureContainer = $ldfeaturecontainer_adaptor->fetch_by_Slice($slice);
  Description: Overwrites superclass method to add the name of the slice to the LDFeatureContainer.
  Returntype : Bio::EnsEMBL::Variation::LDFeatureContainer
  Exceptions : thrown on bad argument
  Caller     : general

=cut

sub fetch_by_Slice{
    my $self = shift;
    my $slice = shift;
    my $population_id = shift;
    my $ldFeatureContainer;
    if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
	throw('Bio::EnsEMBL::Slice arg expected');
    }

    #if a population is passed as an argument, select the LD in the region with the population
    if ($population_id){
	$ldFeatureContainer = $self->generic_fetch("pl.seq_region_id = " . $slice->get_seq_region_id() . " AND pl.seq_region_start >= " . $slice->start() . " AND pl.seq_region_end <= " . $slice->end() .  " AND pl.seq_region_start <= " . $slice->end() . " AND pl.population_id = " . $population_id);
    }
    else{
	$ldFeatureContainer = $self->generic_fetch("pl.seq_region_id = " . $slice->get_seq_region_id() . " AND pl.seq_region_start >= " . $slice->start() . " AND pl.seq_region_end <= " . $slice->end() .  " AND pl.seq_region_start <= " . $slice->end());
    }
    #and store the name of the slice in the Container
    $ldFeatureContainer->name($slice->name());
    return $ldFeatureContainer;
}

=head2 fetch_by_VariationFeature

  Arg [1]    : Bio::EnsEMBL:Variation::VariationFeature $vf
  Example    : my $ldFeatureContainer = $ldFetureContainerAdaptor->fetch_by_VariationFeature($vf);
  Description: Retrieves LDFeatureContainer for a given variation feature.  Most
               variations should only hit the genome once and only a return
               a single variation feature.
  Returntype : reference to Bio::EnsEMBL::Variation::LDFeatureContainer
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub fetch_by_VariationFeature {
  my $self = shift;
  my $vf  = shift;

  if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Bio::EnsEMBL::Variation::VariationFeature arg expected');
  }

  if(!defined($vf->dbID())) {
    throw("VariationFeature arg must have defined dbID");
  }

  my $ldFeatureContainer = $self->generic_fetch("pl.variation_feature_id_1 = ". $vf->dbID() . " OR pl.variation_feature_id_2 = " . $vf->dbID());

  #store the name of the variation feature that generates the object
  $ldFeatureContainer->name($vf->dbID());
  return $ldFeatureContainer;
}

=head2 get_global_population

    Example     : $global_population = $ldContainer->get_global_population();
    Description : Gets the global population used in the ld table to calculate the LD across all populations
    ReturnType  : int $pop_id
    Exceptions  : none
    Caller      : general

=cut

sub get_global_population{
    my $self = shift;
    my $global_population = '';

    my $sth = $self->prepare(qq{SELECT population_id FROM population WHERE name = 'Global population'});
    $sth->execute();
    $sth->bind_columns(\$global_population);
    $sth->fetch;
    $sth->finish;

    return $global_population;
    
}

# methods used by superclass to construct SQL
sub _tables { return ['pairwise_ld', 'pl']; }


sub _columns {
  return qw( pl.variation_feature_id_1 pl.variation_feature_id_2 pl.population_id
	     pl.seq_region_id pl.seq_region_start pl.seq_region_end 
	     pl.snp_distance_count pl.r2 pl.d_prime pl.sample_count);
}

#
# private method, creates ldfeatureContainer objects from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my %feature_container;
  
  my $vfa = $self->db()->get_VariationFeatureAdaptor;

  my %vf_objects; #hash containing the address of all the variation feature objects present in the ld table
  my ($variation_feature_id_1, $variation_feature_id_2, $population_id, $seq_region_id, $seq_region_start, $seq_region_end, $snp_distance_count,
      $r2, $d_prime,$sample_count);

  $sth->bind_columns(\$variation_feature_id_1, \$variation_feature_id_2, \$population_id, \$seq_region_id, \$seq_region_start, \$seq_region_end, \$snp_distance_count, \$r2, \$d_prime, \$sample_count);

  while($sth->fetch()) {
      my %ld_values;
      #get the id of the variations
      $ld_values{'d_prime'} = $d_prime;
      $ld_values{'r2'} = $r2;
      $ld_values{'snp_distance_count'} = $snp_distance_count;      
      $ld_values{'sample_count'} = $sample_count;
      if (!exists $vf_objects{$variation_feature_id_1}){
	  $vf_objects{$variation_feature_id_1} = $vfa->fetch_by_dbID($variation_feature_id_1);
      }      
      if (!exists $vf_objects{$variation_feature_id_2}){
	  $vf_objects{$variation_feature_id_2} = $vfa->fetch_by_dbID($variation_feature_id_2);
      }  
      
      $feature_container{$variation_feature_id_1 . '-' . $variation_feature_id_2}->{$population_id} =  \%ld_values;
  }
  $sth->finish();
  return Bio::EnsEMBL::Variation::LDFeatureContainer->new(
							  '-ldContainer'=> \%feature_container,
							  '-name' => '',
							  '-variationFeatures' => \%vf_objects
							  );
}

1;

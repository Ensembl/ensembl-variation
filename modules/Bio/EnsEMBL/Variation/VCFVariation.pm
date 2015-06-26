=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

# Ensembl module for Bio::EnsEMBL::Variation::VCFVariationFeature
#
#


=head1 NAME

Bio::EnsEMBL::Variation::BaseVariationFeature - Abstract base class for variation features

=head1 SYNOPSIS

None

=head1 DESCRIPTION

Abstract base class representing variation features. Should not be instantiated
directly.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::VCFVariation;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Source;

our @ISA = ('Bio::EnsEMBL::Variation::Variation');

sub new_from_VCFVariationFeature {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  # get and check VF
  my $vf = shift;
  assert_ref($vf, 'Bio::EnsEMBL::Variation::VCFVariationFeature');
  
  my $v = Bio::EnsEMBL::Variation::Variation->new_fast({
    name => $vf->variation_name,
    source => Bio::EnsEMBL::Variation::Source->new_fast({
      name => $vf->vcf_record->{metadata}->{source} || $vf->collection->id,
      description => $vf->vcf_record->{metadata}->{source} || $vf->collection->id,
    }),
  });
  
  bless $v, $class;
  
  $v->{variation_feature} = $vf;
  weaken($v->{variation_feature});
  
  return $v;
}

sub variation_feature {
  my $self = shift;
  $self->{variation_feature} = shift if @_;
  return $self->{variation_feature};
}

sub collection {
  my $self = shift;
  return $self->variation_feature->collection(@_);
}

sub vcf_record {
  my $self = shift;
  return $self->variation_feature->vcf_record(@_);
}

sub get_all_Alleles {
  my $self = shift;
  
  if(!exists($self->{alleles})) {
    
    # get genotypes as a hash keyed on individual dbID
    my %gt_hash = map {$_->individual->dbID() => $_} @{$self->get_all_IndividualGenotypes()};
    
    # get pop/ind hash
    # $hash->{$pop_dbID}->{$ind_dbID} => 1
    my $pop_hash = $self->collection->_get_Population_Individual_hash();
    
    my %pops_by_dbID = map {$_->dbID() => $_} @{$self->collection->get_all_Populations};
    
    my @alleles;
    
    foreach my $pop_id(keys %$pop_hash) {
      
      # get counts
      my %counts;
      $counts{$_}++ for map {@{$_->genotype}} map {$gt_hash{$_}} keys %{$pop_hash->{$pop_id}};
      
      # get total
      my $total = 0;
      $total += $_ for values %counts;
      
      # don't want to divide by zero
      if($total) {
        push @alleles, Bio::EnsEMBL::Variation::Allele->new_fast({
  				allele     => $_,
  				count      => $counts{$_},
  				frequency  => $counts{$_} / $total,
  				population => $pops_by_dbID{$pop_id},
  				variation  => $self,
        }) for keys %counts;
      }
    }
    
    # weaken ref back to self (variation)
    weaken $_->{variation} for @alleles;
    
    $self->{alleles} = \@alleles;
  }
  
  return $self->{alleles};
}

sub get_all_IndividualGenotypes {
  my $self = shift;
  
  if(!exists($self->{genotypes})) {
    
    # use method in VCFCollection to generate genotypes
    my $collection = $self->collection();
    $self->{genotypes} = $collection->_create_IndividualGenotypeFeatures($collection->get_all_Individuals, $self->vcf_record->get_individuals_genotypes, $self->variation_feature);
    
    # weaken ref back to variation feature
    weaken($_->{variation_feature}) for @{$self->{genotypes}};    
  }
   
  return $self->{genotypes};
}

1;
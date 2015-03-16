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

package Bio::EnsEMBL::Variation::VCFVariationFeature;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Variation::VCFVariation;

use Bio::EnsEMBL::Variation::VariationFeature;

our @ISA = ('Bio::EnsEMBL::Variation::VariationFeature');

sub _new_from_VCF_line {
    my $caller = shift;
    my $class = ref($caller) || $caller;
}

sub new_from_VariationFeature {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my ($vf, $collection, $vcf_record, $slice, $adaptor) = rearrange([qw(VARIATION_FEATURE COLLECTION VCF_RECORD SLICE ADAPTOR)], @_);
  
  # get and check VF
  assert_ref($vf, 'Bio::EnsEMBL::Variation::VariationFeature');
  
  bless $vf, $class;
  
  # assert_ref($collection, ('Bio::EnsEMBL::Variation::VCFCollection');
  $vf->{collection} = $collection;
  
  # assert_ref($vcf_record, 'Bio::EnsEMBL::IO::Parser::BaseVCF4');
  $vf->{vcf_record} = $vcf_record;

  # assert_ref($adaptor, 'Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor');
  $vf->{adaptor} = $adaptor;
  
  $vf->{slice} = $slice->seq_region_Slice;
  my $transferred = $vf->transfer($slice);
  
  return $transferred;
}

sub collection {
  my $self = shift;
  $self->{collection} = shift if @_;
  return $self->{collection};
}

sub vcf_record {
  my $self = shift;
  $self->{vcf_record} = shift if @_;
  return $self->{vcf_record};
}

sub variation {
  my $self = shift;
  
  $self->{variation} = shift if @ARGV;
  
  if(!exists($self->{variation})) {
    $self->{variation} = Bio::EnsEMBL::Variation::VCFVariation->new_from_VCFVariationFeature($self);
  }
  
  return $self->{variation};
}

sub get_all_IndividualGenotypes {
  return $_[0]->variation->get_all_IndividualGenotypes;
}

sub minor_allele {
  return $_[0]->variation->minor_allele;
}

sub minor_allele_count {
  return $_[0]->variation->minor_allele_count;
}

sub minor_allele_frequency {
  return $_[0]->variation->minor_allele_frequency;
}

sub source_object {
  return $_[0]->variation->source_object;
}

sub source_name {
  return $_[0]->variation->source_name;
}

sub source_version {
  return $_[0]->variation->source_version;
}

1;
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

# Ensembl module for Bio::EnsEMBL::Variation::VCFVariationFeature
#
#


=head1 NAME

Bio::EnsEMBL::Variation::VCFVariationFeature - A VariationFeature object derived from a VCF

=head1 SYNOPSIS

None

=head1 DESCRIPTION

A child class of Bio::EnsEMBL::Variation::VariationFeature representing an object
derived from a VCF file (via a VCFCollection). Overrides any methods from
parent class that may access data not present in this representation.

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

# sub _new_from_VCF_line {
#     my $caller = shift;
#     my $class = ref($caller) || $caller;
# }


=head2 new_from_VCFVariationFeature
  
  Arg [-VARIATION_FEATURE] :
    Bio::EnsEMBL::Variation::VariationFeature
    Parent class object used to derive this object

  Arg [-COLLECTION] :
    Bio::EnsEMBL::Variation::VCFCollection
    VCFCollection from which this object is generated

  Arg [-VCF_RECORD] :
    arrayref of strings
    Reference to split string as read from VCF file

  Arg [-SLICE] : 
    Bio::EnsEMBL::Slice
    Slice to which this variant maps

  Arg [-ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::VariationFeature
    Adaptor providing database connectivity

  Example    : my $vcf_vf = Bio::EnsEMBL::Variation::VCFVariationFeature->new_from_VariationFeature(
    -VARIATION_FEATURE => $vf,
    -COLLECTION        => $coll,
    -VCF_RECORD        => $record,
    -SLICE             => $slice,
    -ADAPTOR           => $vf_adaptor
  );

  Description: Create a VCFVariationFeature object from a VariationFeature
  Returntype : Bio::EnsEMBL::Variation::VCFVariationFeature
  Exceptions : throws if given VariationFeature is incorrect class
  Caller     : Bio::EnsEMBL::Variation::VCFCollection
  Status     : stable

=cut

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

  # add source object
  $vf->{source} = $collection->source;
    
  # remap to seq region slice
  $vf->{slice} = $slice->seq_region_Slice;
  my $transferred = $vf->transfer($slice);

  my $li = $transferred->location_identifier;
  $transferred->{dbID} = $li;

  # generate name if missing
  $transferred->{variation_name} ||= $li;
  
  return $transferred;
}


=head2 collection

  Example    : my $coll = $v->collection
  Description: Get the VCFCollection object used to generate this object
  Returntype : Bio::EnsEMBL::Variation::VCFCollection
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub collection {
  my $self = shift;
  $self->{collection} = shift if @_;
  return $self->{collection};
}


=head2 vcf_record

  Example    : my $record = $v->vcf_record
  Description: Get the VCF record used to generate this object as an
               array of strings, one per VCF column
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub vcf_record {
  my $self = shift;
  $self->{vcf_record} = shift if @_;
  return $self->{vcf_record};
}


=head2 variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::VCFVariation $v
  Example    : my $vf = $v->variation_feature()
  Description: Get/set the associated Variation object
  Returntype : Bio::EnsEMBL::Variation::VCFVariation
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub variation {
  my $self = shift;
  
  $self->{variation} = shift if @ARGV;
  
  if(!exists($self->{variation})) {
    $self->{variation} = Bio::EnsEMBL::Variation::VCFVariation->new_from_VCFVariationFeature($self);
  }
  
  return $self->{variation};
}


=head2 get_all_SampleGenotypes

  Args       : none
  Example    : $sample_genotypes = $var->get_all_SampleGenotypes()
  Description: Getter for SampleGenotypes for this Variation, returns empty list if 
               there are none 
  Returntype : Listref of SampleGenotypes
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_SampleGenotypes {
  return $_[0]->variation->get_all_SampleGenotypes;
}


=head2 minor_allele

  Example    : $ma = $obj->minor_allele()
  Description: Get the minor allele of this variation
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub minor_allele {
  return $_[0]->variation->minor_allele;
}


=head2 minor_allele_count

  Example    : $maf_count = $obj->minor_allele_count()
  Description: Get the sample count of the minor allele of this variation
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub minor_allele_count {
  return $_[0]->variation->minor_allele_count;
}


=head2 minor_allele_frequency

  Example    : $maf = $obj->minor_allele_frequency()
  Description: Get the frequency of the minor allele of this variation
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub minor_allele_frequency {
  return $_[0]->variation->minor_allele_frequency;
}

1;
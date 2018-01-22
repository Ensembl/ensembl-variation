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

# Ensembl module for Bio::EnsEMBL::Variation::VCFVariation
#
#


=head1 NAME

Bio::EnsEMBL::Variation::VCFVariation - A Variation object derived from a VCF

=head1 SYNOPSIS

None

=head1 DESCRIPTION

A child class of Bio::EnsEMBL::Variation::Variation representing an object
derived from a VCF file (via a VCFCollection). Overrides any methods from
parent class that may access data not present in this representation.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::VCFVariation;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Allele;
use Bio::EnsEMBL::Variation::PopulationGenotype;
use Bio::EnsEMBL::Variation::Source;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(raw_freqs_from_gts);

our @ISA = ('Bio::EnsEMBL::Variation::Variation');


=head2 new_from_VCFVariationFeature
  
  Arg [1]    : Bio::EnsEMBL::Variation::VCFVariationFeature $vf
  Example    : my $v = Bio::EnsEMBL::Variation::VCFVariation->new_from_VCFVariationFeature($vf)
  Description: Create a VCFVariation object from a VCFVariationFeature
  Returntype : Bio::EnsEMBL::Variation::VCFVariation
  Exceptions : throws if $vf is incorrect class
  Caller     : Bio::EnsEMBL::Variation::VCFVariationFeature
  Status     : stable

=cut

sub new_from_VCFVariationFeature {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  # get and check VF
  my $vf = shift;
  assert_ref($vf, 'Bio::EnsEMBL::Variation::VCFVariationFeature');
  
  my $v = Bio::EnsEMBL::Variation::Variation->new_fast({
    name          => $vf->variation_name,
    class_SO_term => $vf->class_SO_term,
    source        => $vf->source,
  });
  
  bless $v, $class;
  
  $v->{variation_feature} = $vf;
  weaken($v->{variation_feature});
  
  return $v;
}


=head2 adaptor
  
  Example    : my $ad = $v->adaptor()
  Description: Gets a database adaptor appropriate for the parent class
  Returntype : Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub adaptor {
  my $self = shift;
  
  if(!exists($self->{adaptor})) {
    $self->{adaptor} = $self->variation_feature->adaptor && $self->variation_feature->adaptor->db ? $self->variation_feature->adaptor->db->get_VariationAdaptor() : undef;
  }
  
  return $self->{adaptor};
}


=head2 variation_feature

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::VCFVariationFeature $vf  
  Example    : my $vf = $v->variation_feature()
  Description: Get/set the associated VariationFeature object
  Returntype : Bio::EnsEMBL::Variation::VCFVariationFeature
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub variation_feature {
  my $self = shift;
  $self->{variation_feature} = shift if @_;
  return $self->{variation_feature};
}


=head2 get_all_VariationFeatures

  Example    : my @vfs = @{$v->get_all_VariationFeatures()}
  Description: Get associated VariationFeature as a listref;
               implemented for compatibility with parent class
  Returntype : arrayref of Bio::EnsEMBL::Variation::VCFVariationFeature
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub get_all_VariationFeatures {
  my $self = shift;
  return [$self->variation_feature(@_)];
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
  return $self->variation_feature->collection(@_);
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
  return $self->variation_feature->vcf_record(@_);
}

=head2 get_all_Alleles

  Example    : @alleles = @{$v->get_all_Alleles()};
  Description: Retrieves all Alleles associated with this variation.
               Alleles are derived from INFO field keys/values or genotypes
               defined in the VCF record
  Returntype : Listref of Bio::EnsEMBL::Variation::Allele objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Alleles {
  my $self = shift;
  
  if(!exists($self->{alleles})) {

    # try to fetch from INFO first
    $self->{alleles} = [
      map {weaken($_->{variation}); $_}
      @{
        $self->collection->get_all_Alleles_by_VariationFeature(
          $self->variation_feature,
          undef, #$_[0] || undef,
          $self->vcf_record
        )
      }
    ];

    # revert to genotypes if that fetches nothing
    if(!@{$self->{alleles}}) {
      $self->{alleles} = [
        map {weaken($_->{variation}); $_}
        map {
          Bio::EnsEMBL::Variation::Allele->new_fast({
            %$_,
            variation => $self,
            frequency_subsnp_handle => ''
          })
        }
        @{$self->_raw_freq_objs_from_genotypes->{Allele}}
      ];
    }
  }
  
  return $self->{alleles};
}


=head2 get_all_PopulationGenotypes

  Example    : $pop_genotypes = $var->get_all_PopulationGenotypes()
  Description: Getter for PopulationGenotypes for this Variation, returns empty list if 
               there are none. 
  Returntype : Listref of PopulationGenotypes
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_PopulationGenotypes {
  my $self = shift;
  
  if(!exists($self->{population_genotypes})) {

    # population genotypes ATM can only be retrieved from genotypes
    # some VCFs have genotype counts in INFO, perhaps investigate in future?
    $self->{population_genotypes} = [
      map {weaken($_->{variation}); $_}
      map {
        Bio::EnsEMBL::Variation::PopulationGenotype->new_fast({
          %$_,
          variation => $self
        })
      }
      @{$self->_raw_freq_objs_from_genotypes->{PopulationGenotype}}
    ];
  }
  
  return $self->{population_genotypes};
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
  my $self = shift;
  
  if(!exists($self->{genotypes})) {

    $self->{genotypes} = [
      map {weaken($_->{variation}); $_}
      @{
        $self->collection->get_all_SampleGenotypeFeatures_by_VariationFeature(
          $self->variation_feature,
          undef, #$_[0] || undef,
          $self->vcf_record
        )
      }
    ];
  }
   
  return $self->{genotypes};
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
  my $self = shift;
  
  if(!exists($self->{minor_allele})) {
    my $counts = $self->_allele_counts;
    $self->{minor_allele} = scalar keys %$counts ? (sort {$counts->{$b} <=> $counts->{$a}} keys %{$counts})[1] : undef;    
  }
  
  return $self->{minor_allele};
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
  my $self = shift;
  
  if(!exists($self->{minor_allele_count})) {
    my $allele = $self->minor_allele;
    $self->{minor_allele_count} = $allele ? $self->_allele_counts->{$allele} : undef;
  }
  
  return $self->{minor_allele_count};
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
  my $self = shift;
  
  if(!exists($self->{minor_allele_frequency})) {
    my $total = 0;
    $total += $_ for values %{$self->_allele_counts};
    
    my $minor_allele_count = $self->minor_allele_count;
    
    $self->{minor_allele_frequency} = ($total && defined($minor_allele_count)) ? $self->_format_frequency($minor_allele_count / $total) : undef;
  }
  
  return $self->{minor_allele_frequency};
}


## stub method for compatibility with parent class
sub get_all_attributes {
  return {};
}


## get allele count data from VCF INFO field
sub _allele_counts {
  my $self = shift;
  
  if(!exists($self->{_allele_counts})) {
    
    my %counts = ();
    
    ## first try and use INFO fields, should be faster
    my $vcf_record = $self->vcf_record();
    my $info = $vcf_record->get_info();
    my @alts = @{$vcf_record->get_alternatives};
    my $ref = $vcf_record->get_reference;

    # convert REF/ALT to Ensembl style
    my %first_bases = map {substr($_, 0, 1) => 1} ($ref, @alts);
    if(scalar keys %first_bases == 1) {
      ($ref, @alts) = map {$_ = substr($_, 1); $_ ||= '-'; $_} ($ref, @alts);
    }
    
    # check we have required fields before attempting
    if(
      $info && ref($info) eq 'HASH' &&
      $info->{AN} &&
      (defined($info->{AC}) || defined($info->{AF})) &&
      $ref && scalar @alts
    ) {
      
      if($info->{AC}) {
        
        # may have one per ALT if multiple ALTS
        my @alt_counts = split(',', $info->{AC});
        
        # check for formatting error
        throw("ERROR: ALT count ".(scalar @alts)." doesn't match AC INFO field (".(scalar @alt_counts).")") unless scalar @alts == scalar @alt_counts;
        
        my $total = 0;
        
        for(my $i = 0; $i <= $#alt_counts; $i++) {
          $counts{$alts[$i]} = $alt_counts[$i];
          $total += $alt_counts[$i];
        }
        
        # get ref count
        $counts{$ref} = $info->{AN} - $total;
      }
      
      elsif(defined($info->{AF})) {
        
        # may have one per ALT if multiple ALTS
        my @alt_freqs = split(',', $info->{AF});
        
        # check for formatting error
        throw("ERROR: ALT count ".(scalar @alts)." doesn't match AF INFO field (".(scalar @alt_freqs).")") unless scalar @alts == scalar @alt_freqs;
        
        my $total = 0;
        
        for(my $i = 0; $i <= $#alt_freqs; $i++) {
          $counts{$alts[$i]} = sprintf("%.0f", $info->{AN} * $alt_freqs[$i]);
          $total += $alt_freqs[$i];
        }
        
        # get ref count
        $counts{$ref} = $info->{AN} - sprintf("%.0f", $info->{AN} * $total);
      }
    }
    
    ## otherwise fetch genotypes to get counts
    else {
      $counts{$_}++ for map {@{$_->genotype}} @{$self->get_all_SampleGenotypes};
    }
    
    $self->{_allele_counts} = \%counts;
  }
  
  return $self->{_allele_counts};
}


## gets intermediate hashrefs used to generate Allele and PopulationGenotype objects
## uses raw_freqs_from_gts util method
sub _raw_freq_objs_from_genotypes {
  my $self = shift;

  if(!exists($self->{_raw_freq_hashes})) {
    my $collection = $self->collection;

    $self->{_raw_freq_hashes} = raw_freqs_from_gts(
      $self->get_all_SampleGenotypes,
      $collection->get_all_Populations,
      $collection->_get_Population_Sample_hash
    );
  }

  return $self->{_raw_freq_hashes};
}


## formats floats to 4dp
sub _format_frequency {
  return sprintf("%.4g", $_[1] || 0);
}


1;
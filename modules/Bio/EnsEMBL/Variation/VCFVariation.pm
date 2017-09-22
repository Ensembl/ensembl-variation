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
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Allele;
use Bio::EnsEMBL::Variation::PopulationGenotype;
use Bio::EnsEMBL::Variation::Source;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(raw_freqs_from_gts);

our @ISA = ('Bio::EnsEMBL::Variation::Variation');

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

sub adaptor {
  my $self = shift;
  
  if(!exists($self->{adaptor})) {
    $self->{adaptor} = $self->variation_feature->adaptor && $self->variation_feature->adaptor->db ? $self->variation_feature->adaptor->db->get_VariationAdaptor() : undef;
  }
  
  return $self->{adaptor};
}

sub variation_feature {
  my $self = shift;
  $self->{variation_feature} = shift if @_;
  return $self->{variation_feature};
}

sub get_all_VariationFeatures {
  my $self = shift;
  return [$self->variation_feature(@_)];
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

sub minor_allele {
  my $self = shift;
  
  if(!exists($self->{minor_allele})) {
    my $counts = $self->_allele_counts;
    $self->{minor_allele} = scalar keys %$counts ? (sort {$counts->{$b} <=> $counts->{$a}} keys %{$counts})[1] : undef;    
  }
  
  return $self->{minor_allele};
}

sub minor_allele_count {
  my $self = shift;
  
  if(!exists($self->{minor_allele_count})) {
    my $allele = $self->minor_allele;
    $self->{minor_allele_count} = $allele ? $self->_allele_counts->{$allele} : undef;
  }
  
  return $self->{minor_allele_count};
}

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

sub _format_frequency {
  return sprintf("%.4g", $_[1] || 0);
}

sub get_all_attributes {
  return {};
}

1;
=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VCFAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::VCFAdaptor

=head1 SYNOPSIS



=head1 DESCRIPTION

This module represents a collection of VCF files, e.g. one for each chromosome
from the 1000 Genomes set. Each file is represented by a
Bio::EnsEMBL::IO::Parser::VCF4Tabix object

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::VCFCollection;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::IO::Parser::VCF4Tabix;
use Bio::EnsEMBL::Variation::IndividualGenotype;

our %TYPES = (
  'remote' => 1,
  'local'  => 1,
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my ($id, $type, $filename_template, $chromosomes, $individual_prefix, $populations, $adaptor) = rearrange([qw(ID TYPE FILENAME_TEMPLATE CHROMOSOMES INDIVIDUAL_PREFIX POPULATIONS ADAPTOR)], @_);
  
  throw("ERROR: No id defined for collection") unless $id;
  throw("ERROR: Collection type $type invalid") unless $type && defined($TYPES{$type});
  
  my %collection = (
    adaptor => $adaptor,
    id => $id,
    type => $type,
    individual_prefix => $individual_prefix,
    #_raw_populations => $populations,
  );
  
  bless(\%collection, $class);
  
  my $obj;
  
  foreach my $chr(@$chromosomes) {
    my $file = $filename_template;
    
    if($file =~ s/\#\#\#CHR\#\#\#/$chr/ || !$obj) {
      $obj = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($file);
    }
    
    $collection{files}->{$chr} = $obj;
  }
  
  return \%collection;
}

sub adaptor {
  my $self = shift;
  $self->{adaptor} = shift if @_;
  return $self->{adaptor};
}

sub id {
  my $self = shift;
  $self->{id} = shift if @_;
  return $self->{id};
}

sub type {
  my $self = shift;
  $self->{type} = shift if @_;
  return $self->{type};
}

sub individual_prefix {
  my $self = shift;
  $self->{individual_prefix} = shift if @_;
  return $self->{individual_prefix};
}

sub list_chromosomes {
  return [keys %{$_[0]->{files}}];
}

sub get_vcf_by_chr {
  return $_[0]->{files}->{$_[1]};
}

sub current {
  my $self = shift;
  $self->{current} = shift if @_;
  return $self->{current};
}

sub seek {
  my $self = shift;
  my ($c, $s, $e) = @_;
  
  my $vcf = $self->get_vcf_by_chr($c);
  return unless $vcf;
  
  # set current to the correct VCF
  $self->current($vcf);
  
  # now seek
  $vcf->seek($c, $s, $e);
  
  return $self->current;
}

sub seek_by_Slice {
  my $self = shift;
  my $slice = shift;
  
  my $vcf = $self->seek($slice->seq_region_name, $slice->start, $slice->end);
  $vcf->next();
  
  return defined($vcf->{record});
}

sub seek_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  
  my $vcf = $self->seek($vf->seq_region_name, $vf->seq_region_start - 1, $vf->seq_region_end + 1);
  
  # compare IDs
  my $count = 0;
  
  while($count++ < 10 && $vcf->next()) {
    last if !$vcf->{record};
    
    # if it has an ID, we can use that
    last if(grep {$vf->variation_name eq $_} @{$vcf->get_IDs});
    
    # otherwise compare coords
    last if $vcf->get_start == $vf->seq_region_start;
    
    # if we've gone too far, quit out
    return 0 if $vcf->get_start > $vf->seq_region_start + 1;
  }
  
  return defined($vcf->{record});
}

sub get_all_Individuals {
  my $self = shift;
  
  if(!defined($self->{individuals})) {
    
    # we should only need to get individuals from one chromosome's VCF
    my $chr = $self->list_chromosomes->[0];
    my $vcf = $self->get_vcf_by_chr($chr);
    
    my $prefix = $self->individual_prefix || '';
    my $ia = $self->adaptor->db->get_IndividualAdaptor();
    
    my $ind_names = $vcf->get_individuals();
    
    # do a fetch_all_by_name_list
    my %ind_objs = map {$_->name() => $_} @{$ia->fetch_all_by_name_list([map {$prefix.$_} @$ind_names])};
    
    # some may not be in DB
    foreach my $ind_name(@$ind_names) {
      
      # either use the DB one or create one
      my $ind = $ind_objs{$prefix.$ind_name} ||
        Bio::EnsEMBL::Variation::Individual->new(
          -name            => $prefix.$ind_name,
          -adaptor         => $ia,
          -type_individual => 'outbred',
          -display         => 'UNDISPLAYABLE',
          #-dbID            => --($self->{_individual_id}),
        );
      
      # store the raw name to easily match to data returned from other methods
      $ind->{_raw_name} = $ind_name;
      
      push @{$self->{individuals}}, $ind;
    }
  }
  
  return $self->{individuals};
}

sub get_all_Populations {
  my $self = shift;
  
  if(!defined($self->{populations})) {
    my $hash = $self->_get_Population_Individual_hash;
    my $pa = $self->adaptor->db->get_PopulationAdaptor();
    $self->{populations} = $pa->fetch_all_by_dbID_list([keys %$hash]);
  }
  
  return $self->{populations};
}

sub _get_Population_Individual_hash {
  my $self = shift;
  
  if(!defined($self->{_population_hash})) {
    
    # populations defined in config?
    if(defined($self->{_raw_populations})) {
      1;
    }
    
    # otherwise we'll have to fetch from the individuals
    else {
      my $inds = $self->get_all_Individuals();
      
      my @dbIDs = grep {defined($_)} map {$_->dbID || undef} @$inds;
      
      my $pa = $self->adaptor->db->get_PopulationAdaptor();
      my $hash = $pa->_get_individual_population_hash(\@dbIDs);
      $self->{_population_hash} = $hash;
    }
  }
  
  return $self->{_population_hash};
}

sub get_all_IndividualGenotypes_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  my $sample = shift;
  
  # seek to record for VariationFeature
  return [] unless $self->seek_by_VariationFeature($vf);
  
  my $vcf = $self->current();
  
  return $self->_create_IndividualGenotypes(
    $self->_limit_Individuals($self->get_all_Individuals, $sample),
    $vcf->get_individuals_genotypes,
    $vf
  );
}

sub get_all_IndividualGenotypes_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $sample = shift;
  
  return [] unless $self->seek_by_Slice($slice);
  
  my $vcf = $self->current();
  
  # get VariationFeatures
  my %vfs_by_pos;
  foreach my $vf(@{$slice->get_all_VariationFeatures()}) {
    push @{$vfs_by_pos{$vf->seq_region_start}}, $vf;
  }
  
  my @genotypes;
  my $individuals = $self->_limit_Individuals($self->get_all_Individuals, $sample);
  
  return [] unless scalar @$individuals;
  
  while($vcf->{record} && $vcf->get_start <= $slice->end) {
    my $start = $vcf->get_raw_start;
    
    # try to match this VCF record to VariationFeature at this position
    my $vf;
    foreach my $tmp_vf(@{$vfs_by_pos{$start} || []}) {
      $vf = $tmp_vf if(grep {$tmp_vf->variation_name eq $_} @{$vcf->get_IDs});
    }
    
    if($vf) {
      push @genotypes, @{$self->_create_IndividualGenotypes($individuals, $vcf->get_individuals_genotypes, $vf)};
    }
    
    $vcf->next();
  }
  
  return \@genotypes;
}

sub _limit_Individuals {
  my $self = shift;
  my $individuals = shift;
  my $sample = shift;
  
  # limit by individual or population?
  if(defined($sample)) {
    if($sample->isa('Bio::EnsEMBL::Variation::Individual')) {
      @$individuals = grep {$_->name eq $sample->name} @$individuals;
    }
    elsif($sample->isa('Bio::EnsEMBL::Variation::Population')) {
      my %limit = map {$_->name => 1} @{$sample->get_all_Individuals};
      @$individuals = grep {defined($limit{$_->name})} @$individuals;
    }
    else {
      throw("Argument $sample is not a Bio::EnsEMBL::Variation::Individual or ::Popluation");
    }
  }
  
  return $individuals;
}

sub _create_IndividualGenotypes {
  my $self = shift;
  my $individuals = shift;
  my $raw_gts = shift;
  my $vf = shift;
  
  my @genotypes;
  
  $self->{_gta} ||= $self->adaptor->db->get_IndividualGenotypeAdaptor;
  
  foreach my $ind(@$individuals) {
    next unless defined($raw_gts->{$ind->{_raw_name}});
    
    my $raw_gt = $raw_gts->{$ind->{_raw_name}};
    my $phased = ($raw_gt =~ /\|/ ? 1 : 0);
    
    my @bits = split(/\||\/|\\/, $raw_gt);
    
    # reverse complement alleles if VF is on -ve strand
    if($vf->seq_region_strand < 0) {
      reverse_comp($$_) for @bits
    }
    
    # create genotype objects
    push @genotypes, Bio::EnsEMBL::Variation::IndividualGenotype->new_fast({
      _variation_id => $vf->{_variation_id},
      variation_feature => $vf,
      individual => $ind,
      genotype => \@bits,
      phased => $phased,
      adaptor => $self->{_gta},
    });
  }
  
  return \@genotypes;
}

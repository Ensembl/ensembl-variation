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
use Bio::EnsEMBL::Utils::Scalar qw(check_ref);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VEP qw(parse_line);

use Bio::EnsEMBL::IO::Parser::VCF4Tabix;
use Bio::EnsEMBL::Variation::IndividualGenotype;
use Bio::EnsEMBL::Variation::Population;

our %TYPES = (
  'remote' => 1,
  'local'  => 1,
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my ($id, $type, $filename_template, $chromosomes, $ind_prefix, $pop_prefix, $ind_pops, $adaptor) = rearrange([qw(ID TYPE FILENAME_TEMPLATE CHROMOSOMES INDIVIDUAL_PREFIX POPULATION_PREFIX INDIVIDUAL_POPULATIONS ADAPTOR)], @_);
  
  throw("ERROR: No id defined for collection") unless $id;
  throw("ERROR: Collection type $type invalid") unless $type && defined($TYPES{$type});
  
  my %collection = (
    adaptor => $adaptor,
    id => $id,
    type => $type,
    individual_prefix => $ind_prefix,
    population_prefix => $pop_prefix,
    chromosomes => $chromosomes,
    filename_template => $filename_template,
    _use_db => 1,
    _raw_populations => $ind_pops,
  );
  
  bless(\%collection, $class);
  
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
  return $self->{individual_prefix} || '';
}

sub population_prefix {
  my $self = shift;
  $self->{population_prefix} = shift if @_;
  return $self->{population_prefix} || '';
}

sub filename_template {
  my $self = shift;
  $self->{filename_template} = shift if @_;
  return $self->{filename_template};
}

sub list_chromosomes {
  return $_[0]->{chromosomes};
}

sub use_db {
  my $self = shift;
  $self->{_use_db} = shift if @_;
  return $self->{_use_db};
}

sub get_vcf_by_chr {
  my $self = shift;
  my $chr = shift;
  
  if(!exists($self->{files}) || !exists($self->{files}->{$chr})) {
    my $obj;
    
    my $file = $self->filename_template;
    
    $file =~ s/\#\#\#CHR\#\#\#/$chr/;
    $obj = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($file);
    
    $self->{files}->{$chr} = $obj;
  }
  
  return $self->{files}->{$chr};
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
    
    my $prefix = $self->individual_prefix;
    my $ia = $self->use_db ? $self->adaptor->db->get_IndividualAdaptor() : undef;
    
    my $ind_names = $vcf->get_individuals();
    
    # do a fetch_all_by_name_list
    my %ind_objs;
    %ind_objs = map {$_->name() => $_} @{$ia->fetch_all_by_name_list([map {$prefix.$_} @$ind_names])} if $self->use_db;
    
    # some may not be in DB
    foreach my $ind_name(@$ind_names) {
      
      # either use the DB one or create one
      my $ind = $ind_objs{$prefix.$ind_name} ||
        Bio::EnsEMBL::Variation::Individual->new(
          -name            => $prefix.$ind_name,
          -adaptor         => $ia,
          -type_individual => 'outbred',
          -display         => 'UNDISPLAYABLE',
          -dbID            => --($self->{_individual_id}),
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

sub get_all_population_names {
  my $self = shift;
  
  if(!exists($self->{_population_names})) {
    
    my @names;
    
    if(defined($self->{_raw_populations})) {
      my $prefix = $self->population_prefix;
      @names =
        map {$prefix.$_}
        map {@{$self->{_raw_populations}->{$_}}}
        keys %{$self->{_raw_populations}};
    }
    
    else {
      @names = map {$_->name} @{$self->get_all_Populations};
    }
    
    $self->{_population_names} = \@names;
  }
  
  return $self->{_population_names};
}

sub has_Population {
  my $self = shift;
  my $pop = shift;
  
  my $name = ref($pop) eq 'SCALAR' ? $pop : $pop->name;
  
  return grep {$name eq $_} @{$self->get_all_population_names};
}

sub _get_Population_Individual_hash {
  my $self = shift;
  
  if(!exists($self->{_population_hash})) {
    my $hash;
    
    # populations defined in config?
    if(defined($self->{_raw_populations})) {
      my $pops;
      my $pa;
      my $prefix = $self->population_prefix;
      
      foreach my $ind(@{$self->get_all_Individuals}) {
        foreach my $pop(@{$self->{_raw_populations}->{$ind->name} || $self->{_raw_populations}->{$ind->{_raw_name}} || []}) {
          
          # try and fetch from DB
          if(!defined($pops->{$pop})) {
            if($self->use_db) {
              $pa ||= $self->adaptor->db->get_PopulationAdaptor();
              $pops->{$pop} = $pa->fetch_by_name($prefix.$pop);
            }
          }
          
          $pops->{$pop} ||= Bio::EnsEMBL::Variation::Population->new_fast({
            name => $prefix.$pop,
            dbID => --($self->{_population_id}),
            _raw_name => $pop,
          });
          
          $hash->{$pops->{$pop}->dbID}->{$ind->dbID} = 1;
          push @{$ind->{populations}}, $pops->{$pop};
        }
        
        $self->{populations} = [values %$pops];
      }
    }
    
    # otherwise we'll have to fetch from the individuals
    else {
      my $inds = $self->get_all_Individuals();
      
      my @dbIDs = grep {defined($_)} map {$_->dbID || undef} @$inds;
      
      my $pa = $self->adaptor->db->get_PopulationAdaptor();
      $hash = $pa->_get_individual_population_hash(\@dbIDs);
    }

    $self->{_population_hash} = $hash;
  }
  
  return $self->{_population_hash};
}

sub get_all_IndividualGenotypeFeatures_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  my $sample = shift;
  
  # seek to record for VariationFeature
  return [] unless $self->seek_by_VariationFeature($vf);
  
  my $vcf = $self->current();
  
  my $individuals = $self->_limit_Individuals($self->get_all_Individuals, $sample);
  my @individual_names = map {$_->{_raw_name}} @$individuals;
  
  return $self->_create_IndividualGenotypeFeatures(
    $individuals,
    $vcf->get_individuals_genotypes(\@individual_names),
    $vf
  );
}

sub get_all_IndividualGenotypeFeatures_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $sample = shift;
  
  return [] unless $self->seek_by_Slice($slice);
  
  my $vcf = $self->current();
  
  # get VariationFeatures if using DB
  my %vfs_by_pos;
  my $use_db = $self->use_db;
  
  if($use_db) {
    foreach my $vf(@{$slice->get_all_VariationFeatures()}) {
      push @{$vfs_by_pos{$vf->seq_region_start}}, $vf;
    }
  }
  
  my @genotypes;
  my $individuals = $self->_limit_Individuals($self->get_all_Individuals, $sample);
  my @individual_names = map {$_->{_raw_name}} @$individuals;
  
  return [] unless scalar @$individuals;
  
  while($vcf->{record} && $vcf->get_start <= $slice->end) {
    my $start = $vcf->get_raw_start;
    
    my $vf;
    
    # try to match this VCF record to VariationFeature at this position
    if($use_db) {
      foreach my $tmp_vf(@{$vfs_by_pos{$start} || []}) {
        $vf = $tmp_vf if(grep {$tmp_vf->variation_name eq $_} @{$vcf->get_IDs});
      }
    }
    
    # otherwise create a VariationFeature object
    else {
      $vf = parse_line({format => 'vcf'}, join("\t", @{$vcf->{record}}))->[0];
    }
    
    if($vf) {
      push @genotypes, @{$self->_create_IndividualGenotypeFeatures($individuals, $vcf->get_individuals_genotypes(\@individual_names), $vf)};
    }
    
    $vcf->next();
  }
  
  return \@genotypes;
}

sub get_all_LD_genotypes_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $sample = shift;
  
  return {} unless $self->seek_by_Slice($slice);
  
  my $vcf = $self->current();
  
  my @genotypes;
  my $individuals = $self->_limit_Individuals($self->get_all_Individuals, $sample);
  my @individual_names = map {$_->{_raw_name}} @$individuals;
  
  my %gts;
  
  while($vcf->{record} && $vcf->get_start <= $slice->end) {
    $gts{$vcf->get_raw_start} = $vcf->get_individuals_genotypes(\@individual_names);
    $vcf->next();
  }
  
  return \%gts;
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

sub _create_IndividualGenotypeFeatures {
  my $self = shift;
  my $individuals = shift;
  my $raw_gts = shift;
  my $vf = shift;
  
  # $vf could be a StructuralVariationFeature if it is generated by parsing the VCF
  return [] unless check_ref($vf, 'Bio::EnsEMBL::Variation::VariationFeature');
  
  my @genotypes;
  
  $self->{_gta} ||= $self->use_db ? $self->adaptor->db->get_IndividualGenotypeFeatureAdaptor : undef;
  
  foreach my $ind(@$individuals) {
    next unless defined($raw_gts->{$ind->{_raw_name}});
    
    my $raw_gt = $raw_gts->{$ind->{_raw_name}};
    my $phased = ($raw_gt =~ /\|/ ? 1 : 0);
    
    my @bits = split(/\||\/|\\/, $raw_gt);
    
    # reverse complement alleles if VF is on -ve strand
    if(($vf->{strand} || $vf->seq_region_strand || 1) < 0) {
      reverse_comp($$_) for @bits
    }
    
    # create genotype objects
    push @genotypes, Bio::EnsEMBL::Variation::IndividualGenotypeFeature->new_fast({
      _variation_id => $vf->{_variation_id},
      variation_feature => $vf,
      individual => $ind,
      genotype => \@bits,
      phased => $phased,
      adaptor => $self->{_gta},
      start => $vf->start,
      end => $vf->end,
      strand => $vf->seq_region_strand,
      slice => $vf->slice,
    });
  }
  
  return \@genotypes;
}

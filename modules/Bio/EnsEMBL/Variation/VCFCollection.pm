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

#
# Ensembl module for Bio::EnsEMBL::Variation::VCFCollection
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::VCFCollection

=head1 SYNOPSIS
  my $reg = 'Bio::EnsEMBL::Registry';

  # get VCF Collection Adaptor
  my $vca = $reg->get_adaptor('human', 'variation', 'vcfcollection');

  # iterate over collections
  foreach my $c(@{$vca->fetch_all})

    # get individuals
    my $individuals = $c->get_all_Individuals;

    # get genotypes for a VariationFeature
    my $gts = $c->get_all_IndividualGenotypeFeatures_by_VariationFeature($vf);
  }

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
use Bio::EnsEMBL::Variation::Individual;
use Bio::EnsEMBL::Variation::Population;

our %TYPES = (
  'remote' => 1,
  'local'  => 1,
);


=head2 new

  Arg [-ID]:                     string - identifier for this collection
  Arg [-TYPE]:                   string - "local" or "remote"
  Arg [-FILENAME_TEMPLATE]:      string
  Arg [-CHROMOSOMES]:            arrayref of strings
  Arg [-INDIVIDUAL_PREFIX]:      string
  Arg [-POPULATION_PREFIX]:      string
  Arg [-INDIVIDUAL_POPULATIONS]: hashref - { 'ind1': ['pop1','pop2'] }
  Arg [-ADAPTOR]:                Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor

  Example    : my $collection = Bio::EnsEMBL::Variation::VCFCollection->new(
                 -id                      => 'test',
                 -type                    => 'local',
                  -filename_template      => '/path/to/vcfs/test_###CHR###.vcf.gz',
                  -chromosomes            => [1, 2, 3],
                  -individual_populations => {
                    'ind1': ['pop1','pop2'],
                    'ind2': ['pop3']
                  }
               );

  Description: Constructor.  Instantiates a new VCFCollection object.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::VCFCollection
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my ($id, $type, $filename_template, $chromosomes, $ind_prefix, $pop_prefix, $ind_pops, $populations, $assembly, $source, $adaptor) = rearrange([qw(ID TYPE FILENAME_TEMPLATE CHROMOSOMES INDIVIDUAL_PREFIX POPULATION_PREFIX INDIVIDUAL_POPULATIONS POPULATIONS ASSEMBLY SOURCE ADAPTOR)], @_);
  
  throw("ERROR: No id defined for collection") unless $id;
  throw("ERROR: Collection type $type invalid") unless $type && defined($TYPES{$type});
 
  if( defined  $source && !$source->isa('Bio::EnsEMBL::Variation::Source')) {
    throw("Bio::EnsEMBL::Variation::Source argument expected");
  }

  my %collection = (
    adaptor => $adaptor,
    id => $id,
    type => $type,
    individual_prefix => $ind_prefix,
    population_prefix => $pop_prefix,
    populations => $populations,
    chromosomes => $chromosomes,
    filename_template => $filename_template,
    assembly  => $assembly,
    source => $source,
    _use_db => 1,
    _raw_populations => $ind_pops,
  );
  
  bless(\%collection, $class);
  
  return \%collection;
}


=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor $adaptor (optional)
               Set the adaptor for this VCFCollection
  Example    : my $adaptor = $collection->adaptor()
  Description: Getter/Setter for the adaptor.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub adaptor {
  my $self = shift;
  $self->{adaptor} = shift if @_;
  return $self->{adaptor};
}


=head2 id

  Arg [1]    : string $id (optional)
               The new value to set the ID attribute to
  Example    : my $id = $collection->id()
  Description: Getter/Setter for the observed count of this allele
               within its associated population.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub id {
  my $self = shift;
  $self->{id} = shift if @_;
  return $self->{id};
}


=head2 type

  Arg [1]    : string $type (optional)
               The new value to set the type attribute to
  Example    : my $type = $collection->type()
  Description: Getter/Setter for the type of this collection
               ('local' or 'remote').
  Returntype : string
  Exceptions : invalid type
  Caller     : general
  Status     : Stable

=cut

sub type {
  my $self = shift;
  
  if(@_) {
    my $type = shift;
    throw("ERROR: Collection type $type invalid") unless $type && defined($TYPES{$type});
    $self->{type} = shift if @_;
  }
  
  return $self->{type};
}


=head2 individual_prefix

  Arg [1]    : string $individual_prefix (optional)
               The new value to set the individual_prefix attribute to
  Example    : my $individual_prefix = $collection->individual_prefix()
  Description: Getter/Setter for the individual_prefix of this collection.
               This property can be useful to align individual names from
               the VCF header with entries in the individual database table.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub individual_prefix {
  my $self = shift;
  $self->{individual_prefix} = shift if @_;
  return $self->{individual_prefix} || '';
}


=head2 population_prefix

  Arg [1]    : string $population_prefix (optional)
               The new value to set the population_prefix attribute to
  Example    : my $population_prefix = $collection->population_prefix()
  Description: Getter/Setter for the population_prefix of this collection.
               This property can be useful to align population names from
               the config file with entries in the population database table.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub population_prefix {
  my $self = shift;
  $self->{population_prefix} = shift if @_;
  return $self->{population_prefix} || '';
}


=head2 filename_template

  Arg [1]    : string $filename_template (optional)
               The new value to set the filename_template attribute to
  Example    : my $filename_template = $collection->filename_template()
  Description: Getter/Setter for the filename template of this collection.
               The wildcard string '###CHR###' can be used in this template
               and will be replaced with the chromosome name when reading,
               allowing a collection to consist of e.g. one VCF file per
               chromosome.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub filename_template {
  my $self = shift;
  $self->{filename_template} = shift if @_;
  return $self->{filename_template};
}


=head2 list_chromosomes

  Example    : my $chrs = $collection->list_chromosomes()
  Description: Get list of chromosome names covered by this collection.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_chromosomes {
  return $_[0]->{chromosomes};
}


=head2 use_db

  Arg [1]    : int $use_db (optional)
               The new value to set the use_db attribute to.
  Example    : my $use_db = $collection->use_db()
  Description: Getter/Setter for the use_db attribute of this collection.
               If set to 1 (default), the API will attempt to retrieve
               Individual and Population objects from the database, using
               individual_prefix and population_prefix as appropriate.
               If set to 0, the API will create "fake" Individual and
               Population objects. Fake objects will also be created if
               the DB fetch fails for an individual.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub use_db {
  my $self = shift;
  $self->{_use_db} = shift if @_;
  return $self->{_use_db};
}


=head2 get_all_Individuals

  Example    : my $inds = $collection->get_all_Individuals()
  Description: Get all individual objects that will have genotypes in
               this collection.
  Returntype : arrayref of Bio::EnsEMBL::Variation::Individual
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Individuals {
  my $self = shift;
  
  if(!defined($self->{individuals})) {
    
    # we should only need to get individuals from one chromosome's VCF
    my $chr = $self->list_chromosomes->[0];
    my $vcf = $self->_get_vcf_by_chr($chr);
    
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
        Bio::EnsEMBL::Variation::Individual->new_fast({
          name            => $prefix.$ind_name,
          adaptor         => $ia,
          type_individual => 'outbred',
          display         => 'UNDISPLAYABLE',
          dbID            => --($self->{_individual_id}),
        });
      
      # store the raw name to easily match to data returned from other methods
      $ind->{_raw_name} = $ind_name;
      
      push @{$self->{individuals}}, $ind;
    }
  }
  
  return $self->{individuals};
}


=head2 get_all_Populations

  Example    : my $pops = $collection->get_all_Populations()
  Description: Get all population objects for the individuals in this
               collection.
  Returntype : arrayref of Bio::EnsEMBL::Variation::Population
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Populations {
  my $self = shift;

  if(!exists($self->{populations}) || !defined $self->{populations}->[0] ) {
    my $hash = $self->_get_Population_Individual_hash;

    if( $self->use_db) {
      my $pa = $self->adaptor->db->get_PopulationAdaptor;
      $self->{populations} = $pa->fetch_all_by_dbID_list([keys %$hash]);
    }
  }
  
  return $self->{populations};
}


=head2 has_Population
  Arg[1]     : string $pop OR Bio::EnsEMBL::Variation::Population $pop
  Example    : my $has_pop = $collection->has_Population($pop)
  Description: Returns true if this collection contains this population
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub has_Population {
  my $self = shift;
  my $pop = shift;
  
  my $name = ref($pop) eq '' ? $pop : $pop->name;
  
  return grep {$name eq $_} @{$self->_get_all_population_names};
}


=head2 get_all_IndividualGenotypeFeatures_by_VariationFeature

  Arg[1]     : Bio::EnsEMBL::Variation::VariationFeature $vf
  Arg[2]     : (optional) Bio::EnsEMBL::Variation::Population OR
               Bio::EnsEMBL::Variation::Individual
  Example    : my $gts = $collection->get_all_IndividualGenotypeFeatures_by_VariationFeature($vf)
  Description: Get all IndividualGenotypeFeatures for a given
               VariationFeature object
  Returntype : arrayref of Bio::EnsEMBL::Variation::IndividualGenotypeFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_IndividualGenotypeFeatures_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  my $sample = shift;
  
  # seek to record for VariationFeature
  return [] unless $self->_seek_by_VariationFeature($vf);
  
  my $vcf = $self->_current();
  
  my $individuals = $self->_limit_Individuals($self->get_all_Individuals, $sample);
  my @individual_names = map {$_->{_raw_name}} @$individuals;
  
  return $self->_create_IndividualGenotypeFeatures(
    $individuals,
    $vcf->get_individuals_genotypes(\@individual_names),
    $vf
  );
}


=head2 get_all_IndividualGenotypeFeatures_by_Slice

  Arg[1]     : Bio::EnsEMBL::Slice $slice
  Arg[2]     : (optional) Bio::EnsEMBL::Variation::Population OR
               Bio::EnsEMBL::Variation::Individual
  Example    : my $gts = $collection->get_all_IndividualGenotypeFeatures_by_Slice($slice)
  Description: Get all IndividualGenotypeFeatures for a given
               genomic region represented by a slice
  Returntype : arrayref of Bio::EnsEMBL::Variation::IndividualGenotypeFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_IndividualGenotypeFeatures_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $sample = shift;
  
  return [] unless $self->_seek_by_Slice($slice);
  
  my $vcf = $self->_current();
  
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
        $vf = $tmp_vf if(grep {$tmp_vf->variation_name eq $_ || $tmp_vf->variation_name eq 'ss'.$_} @{$vcf->get_IDs});
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

sub source {
  my $self = shift;
  
  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Source')) {
      throw("Bio::EnsEMBL::Variation::Source argument expected");
    }
    $self->{'source'} = shift;
  }
  
  return $self->{'source'};
}

sub source_name{
  my $self = shift;
  my $source = $self->{source};
  return unless defined $source;
  
  $source->name(@_) if(@_);
  return $source->name;
}

sub source_url{
  my $self = shift;
  my $source = $self->{source};
  return unless defined $source;

  $source->url(@_) if(@_);
  return $source->url;
}

sub assembly{
  my $self = shift;
  return $self->{assembly};
}

## INTERNAL METHODS
###################

# This method is called in LDFeatureContainerAdaptor::_fetch_by_Slice_VCF.
# It bypasses a lot of the object-creation overhead of
# get_all_IndividualGenotypeFeatures_by_Slice so should return quite a bit faster
sub _get_all_LD_genotypes_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $sample = shift;
  
  return {} unless $self->_seek_by_Slice($slice);
  
  my $vcf = $self->_current();
  
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
    
    # adjust alleles for non-SNVs
    if(defined($vf->{class_SO_term}) && $vf->{class_SO_term} ne 'SNV') {
      my $vcf = $self->_current();
      my @alleles = (($vcf->get_reference),@{$vcf->get_alternatives});
      my %first_char = map {substr($_, 0, 1)} @alleles;
      
      # only do this if the first base is the same in all alleles
      if(scalar keys %first_char == 1) {
        $_ = substr($_, 1) || '-' for @bits;
      }
    }
    
    # create genotype objects
    push @genotypes, Bio::EnsEMBL::Variation::IndividualGenotypeFeature->new_fast({
      _variation_id     => $vf->{_variation_id},
      variation_feature => $vf,
      individual        => $ind,
      genotype          => \@bits,
      phased            => $phased,
      adaptor           => $self->{_gta},
      start             => $vf->start,
      end               => $vf->end,
      strand            => $vf->seq_region_strand,
      slice             => $vf->slice,
    });
  }
  
  return \@genotypes;
}

sub _get_vcf_by_chr {
  my $self = shift;
  my $chr = shift;
  
  if(!exists($self->{files}) || !exists($self->{files}->{$chr})) {
    my $obj;
    
    # check we have this chromosome
    if(my $chrs = $self->list_chromosomes) {
      return unless grep {$chr eq $_} @$chrs;
    }
    
    my $file = $self->filename_template;
    
    $file =~ s/\#\#\#CHR\#\#\#/$chr/;
    $obj = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($file);
    
    $self->{files}->{$chr} = $obj;
  }
  
  return $self->{files}->{$chr};
}

sub _current {
  my $self = shift;
  $self->{current} = shift if @_;
  return $self->{current};
}

sub _seek {
  my $self = shift;
  my ($c, $s, $e) = @_;
  
  my $vcf = $self->_get_vcf_by_chr($c);
  return unless $vcf;
  
  # set current to the correct VCF
  $self->_current($vcf);
  
  # now seek
  $vcf->seek($c, $s, $e);
  
  return $self->_current;
}

sub _seek_by_Slice {
  my $self = shift;
  my $slice = shift;
  
  my $vcf = $self->_seek($slice->seq_region_name, $slice->start, $slice->end);
  return unless $vcf;
  
  $vcf->next();
  
  return defined($vcf->{record});
}

sub _seek_by_VariationFeature {
  my $self = shift;
  my $vf = shift;
  
  my $vcf = $self->_seek($vf->seq_region_name, $vf->seq_region_start - 2, $vf->seq_region_end + 2);
  return unless $vcf;
  
  # compare IDs
  my $count = 0;
  
  while($count++ < 10 && $vcf->next()) {
    last if !$vcf->{record};
    
    # if it has an ID, we can use that
    last if(grep {$vf->variation_name eq $_ || $vf->variation_name eq 'ss'.$_} @{$vcf->get_IDs});
    
    # otherwise compare coords
    last if $vcf->get_start == $vf->seq_region_start;
    
    # if we've gone too far, quit out
    return 0 if $vcf->get_start > $vf->seq_region_start + 1;
  }
  
  return defined($vcf->{record});
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

sub _get_all_population_names {
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

1;

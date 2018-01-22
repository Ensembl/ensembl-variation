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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $sa = $reg->get_adaptor('human', 'core', 'slice');
  $lda = $reg->get_adaptor('human', 'variation', 'ldfeaturecontainer');
  $vfa = $reg->get_adaptor('human', 'variation', 'variationfeature');

  # Get a LDFeatureContainer for a region
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

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Variation::LDFeatureContainer;
use vars qw(@ISA);
use Data::Dumper;

use POSIX;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use constant MAX_SNP_DISTANCE => 100_000;
use constant MIN_R2 => 0.0;
use constant MIN_D_PRIME => 0.0;

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

our $VCF_BINARY_FILE  = '';
our $BINARY_FILE      = '';
our $TMP_PATH         = '';

sub max_snp_distance {
  my $self = shift;
  return $self->{'max_snp_distance'} = shift if(@_);
  return $self->{'max_snp_distance'};
}

=head2 min_r2
  Arg [1]    : $min_r2 (optional)
               The new value to set minimum r2 to. Only return results whose r2 is greater than or equal to min r2.
  Example    : $min_r2 = $ld_feature_container_adaptor->min_r2()
  Description: Getter/Setter for min r2 value
  Returntype : Floating point
  Exceptions : None
  Caller     : General
  Status     : Stable
=cut

sub min_r2 {
  my ($self, $r2) = @_;
  if (defined $r2) {
    $self->{'r2'} = $r2;
  }
  return $self->{'r2'} || MIN_R2;
}

=head2 min_d_prime
  Arg [1]    : $min_d_prime (optional)
               The new value to set minimum d_prime to. Only return results whose d_prime is greater than or equal to min d_prime.
  Example    : $min_d_prime = $ld_feature_container_adaptor->min_d_prime()
  Description: Getter/Setter for min d_prime value
  Returntype : Floating point
  Exceptions : None
  Caller     : General
  Status     : Stable
=cut

sub min_d_prime {
  my ($self, $d_prime) = @_;
  if (defined $d_prime) {
    $self->{'d_prime'} = $d_prime;
  }
  return $self->{'d_prime'} || MIN_D_PRIME;
}

sub vcf_executable {
  my $self = shift;
  $VCF_BINARY_FILE = shift if @_;
  unless( $VCF_BINARY_FILE ) {
    my $binary_name = 'ld_vcf';
    ($VCF_BINARY_FILE) = grep {-e $_} map {"$_/$binary_name"} split /:/,$ENV{'PATH'};
  }
  return $VCF_BINARY_FILE; 
}

sub temp_path {
  my $self = shift;
  $TMP_PATH = shift if @_;
  $TMP_PATH ||= '/tmp'; 
  return $TMP_PATH;
}

=head2 fetch_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               
               The slice to fetch genes on. Assuming it is always correct (in the top level)
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Population $population. Population where 
                we want to select the LD information
  Example    : $ldFeatureContainer = $ldfeaturecontainer_adaptor->fetch_by_Slice($slice);
  Description: Overwrites superclass method to add the name of the slice to the LDFeatureContainer.
  Returntype : Bio::EnsEMBL::Variation::LDFeatureContainer
  Exceptions : thrown on bad argument
  Caller     : general
  Status     : Stable

=cut
sub fetch_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $population = shift;

  if(!ref($slice)) {
    throw('Bio::EnsEMBL::Slice arg or listref of Bio::EnsEMBL::Slice expected');
  }
  
  my $use_vcf = $self->db->use_vcf;
  throw("LD computation requires genotypes from VCF files. Set use_vcf to 1. See ensembl-variation/C_code/README.txt\n") unless $use_vcf;

  my @slice_objects = ();
  my $slice_name = "";

  if (ref $slice eq 'ARRAY') {
    foreach (@$slice) {
      if (!$_->isa('Bio::EnsEMBL::Slice')) {
        throw('Bio::EnsEMBL::Slice arg expected');
      }
      push @slice_objects, $_;
      $slice_name .= "_".$_->name;
    }
  } else {
    if (!$slice->isa('Bio::EnsEMBL::Slice')) {
      throw('Bio::EnsEMBL::Slice arg expected');
    }
    push @slice_objects, $slice;
    $slice_name = $slice->name;
  }

  # check cache
  my $key = join("_",
    $slice_name,
    ($population ? $population->dbID : ""),
    $use_vcf,
    $self->{_vf_pos} || 0,
    $self->min_r2,
    $self->min_d_prime,
    join("-", sort {$a <=> $b} keys %{$self->{_pairwise} || {}})
  );

  return $self->{_cached} if $self->{_cached} && $self->{_cached_key} eq $key;

  my $vcf_container = $self->_fetch_by_Slice_VCF($slice, $population);
  # cache before returning
  $self->{_cached} = $vcf_container;
  $self->{_cached_key} = $key;

  return $vcf_container;
}

=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL:Variation::Variation $v
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Population $pop
  Arg [3]    : (optional) int $max_distance
  Example    : my $ldFeatureContainers = $ldFetureContainerAdaptor->fetch_all_by_Variation($v);
  Description: Retrieves listref of LDFeatureContainers for a given variant.
               If optional population is supplied, values are only returned for that population.
               $max_distance between variant pairs defaults to 100kb
  Returntype : reference to Bio::EnsEMBL::Variation::LDFeatureContainer
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Variation {
  my $self = shift;
  my $v = shift;

  if (!ref($v) || !$v->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation arg expected');
  }

  my $vfs = $v->get_all_VariationFeatures();
  throw('Could not retrieve VariationFeatures (locations) for the given Variation. Include failed variants to return variants with multiple mappings.') if (scalar @$vfs == 0);

  my @containers = ();
  foreach my $vf (@$vfs) {
    my $ldfc = $self->fetch_by_VariationFeature($vf, @_);
    push @containers, $ldfc;
  }
  return \@containers;
}

=head2 fetch_by_VariationFeature

  Arg [1]    : Bio::EnsEMBL:Variation::VariationFeature $vf
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Population $pop
  Arg [3]    : (optional) int $max_distance
  Example    : my $ldFeatureContainer = $ldFetureContainerAdaptor->fetch_by_VariationFeature($vf);
  Description: Retrieves LDFeatureContainer for a given variation feature.
               If optional population is supplied, values are only returned for that population.
               $max_distance between variant pairs defaults to 100kb
  Returntype : reference to Bio::EnsEMBL::Variation::LDFeatureContainer
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_VariationFeature {
  my $self = shift;
  my $vf  = shift;
  my $pop = shift;
  my $distance = shift;

  if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Bio::EnsEMBL::Variation::VariationFeature arg expected');
  }

  if (!$vf->slice->is_reference) {
    warning('Variation feature is not located on the reference sequence but either on a patch or haplotype region.');
    return undef;
  }

  if(!defined($vf->dbID()) && !$vf->isa('Bio::EnsEMBL::Variation::VCFVariationFeature')) {
    throw("VariationFeature arg must have defined dbID");
  }
  
  # cache the position so objs_from_sth picks it up later to filter
  $self->{_vf_pos} = $vf->seq_region_start;
  $self->{_vf_name} = $vf->variation_name;
  
  # fetch by slice using expanded feature slice
  my $max_snp_distance = $distance || $self->{max_snp_distance} || MAX_SNP_DISTANCE;
  my $ldFeatureContainer = $self->fetch_by_Slice($vf->feature_Slice->expand($max_snp_distance, $max_snp_distance), $pop);
  
  # delete the cached pos
  delete $self->{_vf_pos};
  delete $self->{_vf_name};
  
  $ldFeatureContainer->name($vf->dbID);
  
  return $ldFeatureContainer;
}

=head2 fetch_by_VariationFeatures

  Arg [1]    : Listref of Bio::EnsEMBL:Variation::VariationFeature args
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Population $pop
  Example    : my $ldFeatureContainer = $ldFetureContainerAdaptor->fetch_by_VariationFeatures([$vf1, $vf2]);
  Description: Retrieves LDFeatureContainer for a given set of variation features. If optional population is supplied, values are only returned for that population.
  Returntype : reference to Bio::EnsEMBL::Variation::LDFeatureContainer
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_VariationFeatures {
  my $self = shift;
  my $vfs  = shift;
  my $population = shift;

  $DB::single = 1;

  my @slice_objects = ();
  if (!ref($vfs)) {
    throw('Listref of Bio::EnsEMBL::Variation::VariationFeature args expected');
  }
  foreach my $vf (@$vfs) {
    if (!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
      throw('Bio::EnsEMBL::Variation::VariationFeature arg expected');
    }
    if (!defined($vf->dbID())) {
      throw("VariationFeature arg must have defined dbID");
    }
    push @slice_objects, $vf->feature_Slice->expand(1, 1);
  }
  # cache positions
  foreach my $vf (@$vfs) {
    $self->{_pairwise}->{$vf->seq_region_start} = 1;
    $self->{_pairwise_vf_name}->{$vf->variation_name} = 1;
  }  
 
  # fetch by slice using expanded feature slice
  my $ldFeatureContainer = $self->fetch_by_Slice(\@slice_objects, $population);
  
  $ldFeatureContainer->name($vfs->[0]->dbID);
  
  return $ldFeatureContainer;
}

sub _fetch_by_Slice_VCF {
  my $self = shift;
  my $slice = shift;
  my $population = shift;
  my $vca = $self->db->get_VCFCollectionAdaptor();
  
  # fetch genotypes
  my $genotypes = {};

  # create hash mapping positions to variant names
  my %pos2name;

  my $bin = $self->vcf_executable;
  throw("Binary file not found. See ensembl-variation/C_code/README.txt\n") unless $bin;

  my $min_r2 = $self->min_r2;
  my $min_d_prime = $self->min_d_prime;
  my $container;
  my $collections = $vca->fetch_all;

  # get populations
  my @populations = $population ? ($population) : map {@{$_->get_all_Populations}} @$collections;

  foreach my $population (@populations) {
    foreach my $vc (@$collections) {
      my $sample_string = '';
      # skip this collection if it doesn't have the population we want
      if (defined($population)) {
        next unless $vc->has_Population($population);
        my $prefix = $vc->sample_prefix();
        $sample_string = join(",",
          map {$_ =~ s/^$prefix//; $_}
          map {$_->name}
          @{$population->get_all_Samples}
        );
      }

      my $cmd;
      my @files = ();
      my @regions = ();
      my @slices = ();
      if (ref($slice) eq 'ARRAY') { 
        push @slices, @$slice;
      } else {
        push @slices, $slice; 
      }
      foreach my $slice (@slices) {
        my $vcf_file = $vc->_get_vcf_filename_by_chr($slice->seq_region_name);
        throw("ERROR: Can't get VCF file\n") unless $vcf_file;
        push @files, $vcf_file;
        my $loc_string = sprintf("%s:%i-%i", $slice->seq_region_name, $slice->start, $slice->end);
        push @regions, $loc_string;
      }
      my $files_arg = join(',', @files); 
      my $regions_arg = join(',', @regions);
      my $number_of_files = scalar @files;
      $cmd = "$bin -f $files_arg -r $regions_arg -s $number_of_files -l $sample_string";

      if ($self->{_vf_name}) {
        $cmd .= " -v " . $self->{_vf_name};
      }
      # run LD binary and open as pipe
      open LD, "$cmd |"  or die "$!";

      # now create the container from the output of the LD binary
      my %feature_container = ();

      my $population_id = $population ? $population->dbID : 1;
   
      while(<LD>){
        my %ld_values = ();
        #get the ouput into the hashes
        chomp;
        my (
          $null,
          $ld_region_id,
          $ld_region_start,
          $id1,
          $ld_region_end,
          $id2,
          $r2,
          $d_prime,
          $sample_count
        ) = split /\s/;

        # filter by r2 and d_prime values
        next if ($r2 < $min_r2 || $d_prime < $min_d_prime);
        # skip entries unrelated to selected vf if doing fetch_all_by_VariationFeature
        if (defined($self->{_vf_pos})) {
          next unless $ld_region_start == $self->{_vf_pos} || $ld_region_end == $self->{_vf_pos};
        }
        # skip entries unrelated to selected vf if doing fetch_all_by_VariationFeature, exclude co-located variants with same location but different alleles: eg C/T and C/-
        if (defined($self->{_vf_name})) {
          next unless $id1 eq $self->{_vf_name} || $id2 eq $self->{_vf_name};
        }

        # skip entries for pairwise computation that don't match input variation feature loactions
        if (defined $self->{_pairwise}) {
          next unless ($self->{_pairwise}->{$ld_region_start} && $self->{_pairwise}->{$ld_region_end});
          next unless ($self->{_pairwise_vf_name}->{$id1} && $self->{_pairwise_vf_name}->{$id2});
        } 

        $ld_values{'d_prime'} = $d_prime;
        $ld_values{'r2'} = $r2;
        $ld_values{'sample_count'} = $sample_count;

        $id1 =~ s/\;.+//;
        $id2 =~ s/\;.+//;
        $pos2name{$ld_region_start} = $id1;
        $pos2name{$ld_region_end} = $id2;
        $feature_container{$ld_region_start . '-' . $ld_region_end}->{$population_id} = \%ld_values;
      }

      # Close the file handle per iteration, don't reuse the
      # glob without closing it.
      close LD;

      my $c = Bio::EnsEMBL::Variation::LDFeatureContainer->new(
        '-adaptor' => $self,
        '-ldContainer'=> \%feature_container,
        '-name' => '',
        '-slices' => [$slice],
      );
      $c->{'_vf_name'} = $self->{'_vf_name'};
      $c->{'_pop_ids'} = {$population_id => 1};

      if($container) {
        $self->_merge_containers($container, $c);
      }
      else {
        $container = $c;
      }
    }
  }

  $container->{pos2name} = \%pos2name if $container;
  delete $self->{_pairwise};
  delete $self->{_pairwise_vf_name};

  if (!$container) {
    warning('The population is not represented in the configured VCF file for fetching genotypes for LD computation.');
    return  Bio::EnsEMBL::Variation::LDFeatureContainer->new('-adaptor' => $self, '-ldContainer' => {}, 'name' => '', '-slices' => []);
  }
  return $container;
}

sub _merge_containers {
  my $self = shift;
  my $c1 = shift;
  my $c2 = shift;
  # merge VFs
  $c1->{variationFeatures}->{$_} ||= $c2->{variationFeatures}->{$_} for keys %{$c2->{variationFeatures} || {}};
  
  # merge pop IDs
  $c1->{_pop_ids}->{$_} ||= $c2->{_pop_ids}->{$_} for keys %{$c2->{_pop_ids} || {}};

  # merge pos2name
  $c1->{pos2name}->{$_} ||= $c2->{pos2name}->{$_} for keys %{$c2->{pos2name} || {}};

  # if both have pos2vf, merge
  if($c1->{pos2vf} && $c2->{pos2vf}) {
    $c1->{pos2vf}->{$_} ||= $c2->{pos2vf}->{$_} for keys %{$c2->{pos2vf} || {}};
  }
  # otherwise we have to delete and rely on the container to lazy load
  elsif($c1->{pos2vf}) {
    delete $c1->{pos2vf};
  }
  elsif($c2->{pos2vf}) {
    delete $c2->{pos2vf};
  }

  # otherwise we have to remove it
  
  # merge ldContainer
  foreach my $pair(keys %{$c2->{ldContainer} || {}}) {
    if($c1->{ldContainer}->{$pair}) {
      $c1->{ldContainer}->{$pair}->{$_} ||= $c2->{ldContainer}->{$pair}->{$_} for keys %{$c2->{ldContainer}->{$pair}};
    }
    else {
      $c1->{ldContainer}->{$pair} = $c2->{ldContainer}->{$pair};
    }
  }

  return $c1;
}

sub get_populations_hash_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $sth = $self->prepare(qq{SELECT population_id, name FROM population WHERE display = 'LD';});
  $sth->execute;
  my %results = map {$_->[0] => $_->[1]} @{$sth->fetchall_arrayref()};
  return \%results;
}

1;

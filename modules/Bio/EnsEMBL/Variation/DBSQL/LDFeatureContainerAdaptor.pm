=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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
use FileHandle;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use constant MAX_SNP_DISTANCE => 100_000;

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

#our $MAX_SNP_DISTANCE = 100000;
our $VCF_BINARY_FILE  = '';
our $BINARY_FILE      = '';
our $TMP_PATH         = '';

sub max_snp_distance {
  my $self = shift;
  return $self->{'max_snp_distance'} = shift if(@_);
  return $self->{'max_snp_distance'};
}

sub executable {
  my $self = shift;
  $BINARY_FILE = shift if @_;
  unless( $BINARY_FILE ) {
    my $binary_name = 'calc_genotypes';
    ($BINARY_FILE) = grep {-e $_} map {"$_/$binary_name"} split /:/,$ENV{'PATH'};
  }
  return $BINARY_FILE; 
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
  
  my @slice_objects = ();
  my $slice_name = "";
  my $use_vcf = $self->db->use_vcf;
  throw("LD computation requires genotypes from VCF files. Set use_vcf to 1 or 2. See ensembl-variation/C_code/README.txt\n") unless $use_vcf;

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
    join("-", sort {$a <=> $b} keys %{$self->{_pairwise} || {}})
  );
  return $self->{_cached} if $self->{_cached} && $self->{_cached_key} eq $key;
  my @genotypes = ();
  my $vcf_container;

  # use VCF?
  if ($use_vcf) {
    $vcf_container = $self->_fetch_by_Slice_VCF($slice, $population);

    if($use_vcf > 1 || ($population && $vcf_container)) {

      # cache before returning
      $self->{_cached} = $vcf_container;
      $self->{_cached_key} = $key;

      return $vcf_container;
    }
  } 
  
  my $siblings = {};
  #when there is no population selected, return LD in the HapMap and PerlEgen populations
  my $in_str;
  
  #if a population is passed as an argument, select the LD in the region with the population
  if ($population) {

    if (!ref($population) || !$population->isa('Bio::EnsEMBL::Variation::Population')) {
      throw('Bio::EnsEMBL::Variation::Population arg expected');
    }
    my $population_id = $population->dbID;

    $in_str = " = $population_id";
  }
  else {
    $in_str = $self->_get_LD_populations($siblings);
  }

  my $ldFeatureContainer = Bio::EnsEMBL::Variation::LDFeatureContainer->new(
      '-adaptor' => $self,
      '-ldContainer'=> {},
      '-name' => $slice_objects[0]->name,
      '-variationFeatures' => {}
    );
  
  $ldFeatureContainer = $self->_merge_containers($vcf_container, $ldFeatureContainer) if $vcf_container;

  # cache
  $self->{_cached} = $ldFeatureContainer;
  $self->{_cached_key} = $key;

  return $ldFeatureContainer;
}

=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL:Variation::Variation $v
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Population $pop
  Example    : my $ldFeatureContainers = $ldFetureContainerAdaptor->fetch_all_by_Variation($v);
  Description: Retrieves listref of LDFeatureContainers for a given variant. If optional population is supplied, values are only returned for that population.
  Returntype : reference to Bio::EnsEMBL::Variation::LDFeatureContainer
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Variation {
  my $self = shift;
  my $v = shift;
  my $pop = shift;

  if (!ref($v) || !$v->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation arg expected');
  }

  my $vfs = $v->get_all_VariationFeatures();
  throw('Could not retrieve VariationFeatures (locations) for the given Variation. Include failed variants to return variants with multiple mappings.') if (scalar @$vfs == 0);

  my @containers = ();
  foreach my $vf (@$vfs) {
    my $ldfc = $self->fetch_by_VariationFeature($vf, $pop);
    push @containers, $ldfc;
  }
  return \@containers;
}

=head2 fetch_by_VariationFeature

  Arg [1]    : Bio::EnsEMBL:Variation::VariationFeature $vf
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Population $pop
  Example    : my $ldFeatureContainer = $ldFetureContainerAdaptor->fetch_by_VariationFeature($vf);
  Description: Retrieves LDFeatureContainer for a given variation feature. If optional population is supplied, values are only returned for that population.
  Returntype : reference to Bio::EnsEMBL::Variation::LDFeatureContainer
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_VariationFeature {
  my $self = shift;
  my $vf  = shift;
  my $pop = shift;

  if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Bio::EnsEMBL::Variation::VariationFeature arg expected');
  }

  if (!$vf->slice->is_reference) {
    warning('Variation feature is not located on the reference sequence but either on a patch or haplotype region.');
    return undef;
  }

  if(!defined($vf->dbID())) {
    throw("VariationFeature arg must have defined dbID");
  }
  
  # cache the position so objs_from_sth picks it up later to filter
  $self->{_vf_pos} = $vf->seq_region_start;
  $self->{_vf_name} = $vf->variation_name;
  
  # fetch by slice using expanded feature slice
  my $max_snp_distance = $self->{max_snp_distance} || MAX_SNP_DISTANCE;
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
  my $pop = shift;

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
  my $ldFeatureContainer = $self->fetch_by_Slice(\@slice_objects, $pop);
  
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

  my $container;
  my $collections = $vca->fetch_all;

  # get populations
  my @populations = $population ? ($population) : map {@{$_->get_all_Populations}} @$collections;

  foreach my $population(@populations) {
    
    foreach my $vc(@$collections) {
      
      my $sample_string = '';
   
      # skip this collection if it doesn't have the population we want
      if(defined($population)) {
        next unless $vc->has_Population($population);

        my $prefix = $vc->sample_prefix();
       
        $sample_string = join(",",
          map {$_ =~ s/^$prefix//; $_}
          map {$_->name}
          @{$population->get_all_Samples}
        );
      }

      my $cmd;

      # two slices
      if(ref($slice) eq 'ARRAY') {
        my $vcf_file_1 = $vc->_get_vcf_filename_by_chr($slice->[0]->seq_region_name);
        my $vcf_file_2 = $vc->_get_vcf_filename_by_chr($slice->[1]->seq_region_name);

        throw("ERROR: Can't get VCF file\n") unless $vcf_file_1;
        throw("ERROR: Can't get VCF file\n") unless $vcf_file_2;

        my $loc_string_1 = sprintf("%s:%i-%i", $slice->[0]->seq_region_name, $slice->[0]->start, $slice->[0]->end);
        my $loc_string_2 = sprintf("%s:%i-%i", $slice->[1]->seq_region_name, $slice->[1]->start, $slice->[1]->end);

        $cmd = "$bin -f $vcf_file_1 -r $loc_string_1 -g $vcf_file_2 -s $loc_string_2 -l $sample_string";
      }
      # one slice
      else {
        my $vcf_file = $vc->_get_vcf_filename_by_chr($slice->seq_region_name);

        throw("ERROR: Can't get VCF file\n") unless $vcf_file;

        my $loc_string = sprintf("%s:%i-%i", $slice->seq_region_name, $slice->start, $slice->end);

        $cmd = "$bin -f $vcf_file -r $loc_string -l $sample_string";
      }
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


#for a given population, gets all samples that are children (have father or mother)
sub _get_siblings {
  my $self = shift;
  my $population_id = shift;
  my $siblings = shift;

  my $sth_sample = $self->db->dbc->prepare(qq{
    SELECT s.sample_id
    FROM sample s, individual i, sample_population sp
    WHERE sp.sample_id = s.sample_id
    AND s.individual_id = i.individual_id
    AND sp.population_id = ? 
    AND i.father_individual_id IS NOT NULL
    AND i.mother_individual_id IS NOT NULL
  });
  
  my ($sample_id);
  $sth_sample->execute($population_id);
  $sth_sample->bind_columns(\$sample_id);
  
  while ($sth_sample->fetch){
    # store population and sample since some samples are shared between populations
    $siblings->{$population_id.'-'.$sample_id}++;
  }
}

sub _get_LD_populations {
  my $self = shift;
  my $siblings = shift;
  my ($pop_id,$population_name);
  my $sth = $self->db->dbc->prepare(qq{SELECT population_id, name FROM population WHERE display = 'LD'});

  $sth->execute();
  $sth->bind_columns(\$pop_id,\$population_name);
  
  #get all the children that we do not want in the genotypes
  my @pops;
  while($sth->fetch){
    if($population_name =~ /CEU|YRI|MEX/){
      $self->_get_siblings($pop_id,$siblings);
    }
    push @pops, $pop_id;
  }
    
  my $in_str = " IN (" . join(',', @pops). ")";
	
  return $in_str if (defined $pops[0]);
  return '' if (!defined $pops[0]);
}

sub get_populations_hash_by_Slice {
  my $self = shift;
  my $slice = shift;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }

  my $pop_list = $self->_get_LD_populations();

  my ($sr, $slice_start, $slice_end) = ($slice->get_seq_region_id, $slice->start, $slice->end);

  my %results;

  my $sth = $self->prepare(qq{SELECT population_id, name FROM population WHERE population_id $pop_list;});
  $sth->execute;

  %results = map {$_->[0] => $_->[1]} @{$sth->fetchall_arrayref()};
  return \%results;

}

1;

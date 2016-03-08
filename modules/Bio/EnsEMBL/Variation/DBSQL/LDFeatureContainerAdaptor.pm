=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
    ($BINARY_FILE) = grep {-e $_} map {"$_/calc_genotypes"} split /:/,$ENV{'PATH'};
  }
  return $BINARY_FILE; 
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

  if (!ref($slice)) {
    throw('Bio::EnsEMBL::Slice arg or listref of Bio::EnsEMBL::Slice expected');
  }
  
  my @slice_objects = ();

  if (ref $slice eq 'ARRAY') {
    foreach (@$slice) {
      if (!$_->isa('Bio::EnsEMBL::Slice')) {
        throw('Bio::EnsEMBL::Slice arg expected');
      }
      push @slice_objects, $_;
    }
  } else {
    if (!$slice->isa('Bio::EnsEMBL::Slice')) {
      throw('Bio::EnsEMBL::Slice arg expected');
    }
    push @slice_objects, $slice;
  }

  my @genotypes = ();

  # use VCF?
  if ($self->db->use_vcf) {
    if (!defined($population)) {
      foreach my $slice (@slice_objects) {
        push @genotypes, $self->_fetch_by_Slice_VCF($slice);
      }
      if ($self->db->use_vcf > 1) {
        my $ldFeatureContainer = $self->_ld_calc(\@genotypes);
        $ldFeatureContainer->name($slice_objects[0]->name());
        return $ldFeatureContainer;
      }
    }
    elsif (grep {$_} map {$_->has_Population($population)} @{$self->db->get_VCFCollectionAdaptor->fetch_all}) {
      foreach my $slice (@slice_objects) {    
        push @genotypes, $self->_fetch_by_Slice_VCF($slice, $population);
      }
      my $ldFeatureContainer = $self->_ld_calc(\@genotypes);
      $ldFeatureContainer->name($slice_objects[0]->name());
      return $ldFeatureContainer;
    }
  } 
  
  my $siblings = {};
  #when there is no population selected, return LD in the HapMap and PerlEgen populations
  my $in_str = $self->_get_LD_populations($siblings);
  
  #if a population is passed as an argument, select the LD in the region with the population
  if ($population) {
    if (!ref($population) || !$population->isa('Bio::EnsEMBL::Variation::Population')) {
      throw('Bio::EnsEMBL::Variation::Population arg expected');
    }
    my $population_id = $population->dbID;
    $in_str = " = $population_id";
  }

  if ($in_str eq '') {
    #there is no population, not a human specie or not passed as an argument, return the empy container
    my $empty_container = Bio::EnsEMBL::Variation::LDFeatureContainer->new(
      '-ldContainer'=> {},
      '-name' => $slice_objects[0]->name,
      '-variationFeatures' => {}
    );
    return $empty_container;
  }

  foreach my $slice (@slice_objects) {
    push @genotypes, $self->_fetch_by_Slice_DB($slice, $in_str, $siblings);
  }

  my $ldFeatureContainer = $self->_ld_calc(\@genotypes);
  $ldFeatureContainer->name($slice_objects[0]->name());
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

  if(!defined($vf->dbID())) {
    throw("VariationFeature arg must have defined dbID");
  }
  
  # cache the position so objs_from_sth picks it up later to filter
  $self->{_vf_pos} = $vf->seq_region_start;
  
  # fetch by slice using expanded feature slice
  my $max_snp_distance = $self->{max_snp_distance} || MAX_SNP_DISTANCE;
  my $ldFeatureContainer = $self->fetch_by_Slice($vf->feature_Slice->expand($max_snp_distance, $max_snp_distance), $pop);
  
  # delete the cached pos
  delete $self->{_vf_pos};
  
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
  }  
 
  # fetch by slice using expanded feature slice
  my $ldFeatureContainer = $self->fetch_by_Slice(\@slice_objects, $pop);
  
  $ldFeatureContainer->name($vfs->[0]->dbID);
  
  return $ldFeatureContainer;
}

sub get_populations_by_Slice {
  my $self = shift;
  
  my $population_hash = $self->get_populations_hash_by_Slice(@_);
  return [values(%{$population_hash})];
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
  my $pop_threshold = 20;	# number of samples required
  my $gen_threshold = 3;	# number of genotypes per sample required

  # just get the population list if it's really too long
  if($slice->length > 10000000) {
	
    my $sth = $self->prepare(qq{SELECT population_id, name FROM population WHERE population_id $pop_list;});
    $sth->execute;
	
    %results = map {$_->[0] => $_->[1]} @{$sth->fetchall_arrayref()};
  }

  # do a guesstimate for long slices, otherwise it takes too long
  elsif($slice->length > 5000000) {
	
    my $sth = $self->prepare(qq{
      SELECT distinct(c.sample_id), p.name
      FROM compressed_genotype_region c, sample_population ip, population p, sample s, individual i
      WHERE c.sample_id = sp.sample_id
      AND sp.population_id = p.population_id
      AND c.sample_id = s.sample_id
      AND s.individual_id = i.individual_id
      AND c.seq_region_id = ?
      AND c.seq_region_start <= ?
      AND c.seq_region_end >= ?
      AND i.father_individual_id is NULL AND i.mother_individual_id is NULL
      AND (p.population_id $pop_list)
    });
	
    $sth->execute($sr, $slice_end, $slice_start);
	
    my %counts = ();
	
    while(my $row = $sth->fetchrow_arrayref()) {
      my ($sample_id, $pop_id, $pop_name) = @$row;
      $results{$pop_id} = $pop_name;
      $counts{$pop_id}++;
    }
	
    #Delete the populations that don't have enough samples
    delete @results{ grep {$counts{$_} <= $pop_threshold} keys(%counts)};
  }
  
  else {

    my $sth = $self->prepare(qq{
      SELECT p.population_id, p.name, c.sample_id, c.seq_region_start, c.seq_region_end, c.genotypes 
      FROM compressed_genotype_region c, sample_population sp, population p, sample s, individual i
      WHERE c.sample_id = sp.sample_id
      AND sp.population_id = p.population_id
      AND c.sample_id = s.sample_id
      AND s.individual_id = i.individual_id
      AND c.seq_region_id = ?
      AND c.seq_region_start <= ?
      AND c.seq_region_end >= ?
      AND i.father_individual_id is NULL AND i.mother_individual_id is NULL
      AND (p.population_id $pop_list)
    });
	
    $sth->execute($sr, $slice_end, $slice_start);
	
    my (%enough, %counts, %sample_pop, %counts_pop);
	
    my $row_count = 0;
	
    while(my $row = $sth->fetchrow_arrayref()) {
      my ($population_id, $population_name, $sample_id, $start, $end, $genotypes) = @$row;
	  
      $row_count++;
	  
      next if $enough{$sample_id};
	  
      $results{$population_id} = $population_name;
      $sample_pop{$population_id} = $population_name;
	  
      # if the row is only partially within the slice
      if($start < $slice_start || $end > $slice_end) {
		
        my @genotypes = unpack("(www)*", $genotypes);
        my $snp_start = $start;
		
        while( my( $variation_id, $gt_code, $gap ) = splice @genotypes, 0, 3 ) {
          if(($snp_start >= $slice_start) && ($snp_start <= $slice_end)) {		
            $counts{$sample_id}++;
          }
		  
          $snp_start += $gap + 1 if defined $gap;
          last if $snp_start > $slice_end;
        }
      }
	  
      # if the row is fully within the slice
      else {
        $counts{$sample_id} += (((length($genotypes) - 2) / 4) + 1);
      }
	  
      $enough{$sample_id} = 1 if $counts{$sample_id} >= $gen_threshold;
      $counts_pop{$population_id}++ if $counts{$sample_id} >= $gen_threshold;
    }
	
    delete @results{grep {$counts_pop{$_} <= $pop_threshold} keys %counts_pop};
  }
  
  return \%results;
}

sub _fetch_by_Slice_DB {
  my $self = shift;
  my $slice = shift;
  my $in_str = shift;
  my $siblings = shift;
  my $sth = $self->prepare(qq{
    SELECT c.sample_id, c.seq_region_id, c.seq_region_start, c.seq_region_end, c.genotypes, sp.population_id
    FROM compressed_genotype_region c, sample_population sp
    WHERE  sp.sample_id = c.sample_id
    AND   sp.population_id $in_str
    AND   c.seq_region_id = ?
    AND   c.seq_region_start >= ? and c.seq_region_start <= ?
    AND   c.seq_region_end >= ?
    ORDER BY c.seq_region_id, c.seq_region_start
  }, {mysql_use_result => 1});

  $sth->bind_param(1,$slice->get_seq_region_id,SQL_INTEGER);
  $sth->bind_param(2,$slice->start - MAX_SNP_DISTANCE,SQL_INTEGER) if ($slice->start - MAX_SNP_DISTANCE >= 1);
  $sth->bind_param(2,1,SQL_INTEGER) if ($slice->start - MAX_SNP_DISTANCE < 1);
  $sth->bind_param(3,$slice->end,SQL_INTEGER);
  $sth->bind_param(4,($slice->start < 1 ? 1 : $slice->start),SQL_INTEGER);
  $sth->execute();
  $self->_objs_from_sth($sth, $slice, $siblings);
}

sub _fetch_by_Slice_VCF {
  my $self = shift;
  my $slice = shift;
  my $population = shift;
  my $vca = $self->db->get_VCFCollectionAdaptor();
  
  # fetch genotypes
  my $genotypes = {};
    
  # create hash giving populations for each sample
  my %pops;

  # create hash mapping positions to variant names
  my %pos2name;
  
  foreach my $vc(@{$vca->fetch_all}) {
    
    # skip this collection if it doesn't have the population we want
    if(defined($population)) {
      next unless $vc->has_Population($population);
    }
    
    # get "raw" genotypes; comes back as a hash like $hash->{$pos}->{$ind_name} = $gt
    # doing this saves constructing objects we don't need e.g. Genotypes, VariationFeatures
    my ($vc_genotypes, $vc_pos2name) = @{$vc->_get_all_LD_genotypes_by_Slice($slice, $population)};
 
   
    # copy them to main $genotypes hash
    foreach my $p(keys %$vc_genotypes) {
      $genotypes->{$p}->{$_} = $vc_genotypes->{$p}->{$_} for keys %{$vc_genotypes->{$p}};
    }

    # and copy position map
    $pos2name{$_} = $vc_pos2name->{$_} for keys %$vc_pos2name;
    
    # get Population->Sample hash; we need to trim and transpose this
    my $hash = $vc->_get_Population_Sample_hash();
    
    # copy hash before deleting from it
    $hash = { %$hash };
    
    if(defined($population)) {
      delete $hash->{$_} for grep {$_ != $population->dbID} keys %$hash;
    }
    
    # get all samples
    my %sample_dbID_name = map {$_->dbID => ($_->{_raw_name} || $_->name)} @{$vc->get_all_Samples()};
    
    # populate transposed hash
    foreach my $pop_id(keys %$hash) {
      $pops{$sample_dbID_name{$_}}{$pop_id} = 1 for keys %{$hash->{$pop_id}};
    }
  }
  return $self->_objs_from_sth_vcf($genotypes, $slice, \%pops, \%pos2name);
}

sub _merge_containers {
  my $self = shift;
  my $c1 = shift;
  my $c2 = shift;
  
  # merge VFs
  $c1->{variationFeatures}->{$_} ||= $c2->{variationFeatures}->{$_} for keys %{$c2->{variationFeatures}};
  
  # merge pop IDs
  $c1->{_pop_ids}->{$_} ||= $c2->{_pop_ids}->{$_} for keys %{$c2->{_pop_ids}};
  
  # merge ldContainer
  foreach my $pair(keys %{$c2->{ldContainer}}) {
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

#reads one line from the compress_genotypes table, uncompress the data, and writes it to the different hashes: one containing the number of bases for the variation and the other with the actual genotype information we need to print in the file
sub _store_genotype{
  my $self = shift;
  my $sample_information = shift;
  my $alleles_variation = shift;
  my $sample_id = shift;
  my $seq_region_start = shift;
  my $genotype = shift;
  my $population_id = shift;
  my $slice = shift;
	
  my $snp_start = $seq_region_start;
  my ($slice_start, $slice_end) = ($slice->start, $slice->end);
	
  my @genotypes = unpack("(www)*", $genotype);
  while( my( $variation_id, $gt_code, $gap ) = splice @genotypes, 0, 3 ) {
    my $gt = $self->_genotype_from_code($gt_code);
		
    if(
    defined $gt &&
    ($snp_start >= $slice_start) &&
    ($snp_start <= $slice_end) &&
    $gt->[0] =~ /^[ACGT]$/ &&
    $gt->[1] =~ /[ACGT]$/
    ) {		
      $alleles_variation->{$snp_start}->{$population_id}->{$gt->[0]}++;
      $alleles_variation->{$snp_start}->{$population_id}->{$gt->[1]}++;
			
      $sample_information->{$population_id}->{$snp_start}->{$sample_id}->{allele_1} = $gt->[0];
      $sample_information->{$population_id}->{$snp_start}->{$sample_id}->{allele_2} = $gt->[1];
    }
		
    $snp_start += $gap + 1  if defined $gap;
  }
}

sub _genotype_from_code {
  my $self = shift;
  my $code = shift;
	
  return $self->{'_genotype_codes'}->{$code} if defined($self->{'_genotype_codes'}) && defined($self->{'_genotype_codes'}->{$code});
	
  if(!defined($self->{'_gtca'})) {
    $self->{'_gtca'} = $self->db->get_GenotypeCodeAdaptor;
  }
  my $gtca = $self->{'_gtca'};
	
  my $gtc = $gtca->fetch_by_dbID($code);
  $self->{'_genotype_codes'}->{$code} = $gtc ? $gtc->genotype : undef;
	
  return $self->{'_genotype_codes'}->{$code};
}

#
# Converts the genotype into the required format for the calculation of the pairwise_ld value: AA, Aa or aa
# From the Allele table, will select the alleles and compare to the alleles in the genotype
#

sub _convert_genotype{
  my $self = shift;
  my $alleles_variation = shift; #reference to the hash containing the alleles for the variation present in the genotypes
  my $sample_information = shift; #reference to a hash containing the values to be written to the file
  my @alleles_ordered; #the array will contain the alleles ordered by apparitions in the genotypes (only 2 values possible)
    
  @alleles_ordered = sort({$alleles_variation->{$b} <=> $alleles_variation->{$a}} keys %{$alleles_variation});
    
  #let's convert the allele_1 allele_2 to a genotype in the AA, Aa or aa format, where A corresponds to the major allele and a to the minor
  foreach my $sample_id (keys %{$sample_information}){
    #if both alleles are different, this is the Aa genotype
    if ($sample_information->{$sample_id}{allele_1} ne $sample_information->{$sample_id}{allele_2}){
      $sample_information->{$sample_id}{genotype} = 'Aa';
    }
    #when they are the same, must find out which is the major
    else{	    
      if ($alleles_ordered[0] eq $sample_information->{$sample_id}{allele_1}){
        #it is the major allele
        $sample_information->{$sample_id}{genotype} = 'AA';
      }
      else{
        $sample_information->{$sample_id}{genotype} = 'aa';
      }
	    
    }
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


sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;
  my $slice = shift;
  my $siblings = shift;
  
  my %alleles_variation = (); #will contain a record of the alleles in the variation. A will be the major, and a the minor. When more than 2 alleles
  #, the genotypes for that variation will be discarded
  my %sample_information = (); #to store the relation of snps->samples

  my ($sample_id, $seq_region_id, $seq_region_start, $seq_region_end, $genotypes, $population_id);
  my $previous_seq_region_id = 0;
  
  $sth->bind_columns(\$sample_id, \$seq_region_id, \$seq_region_start, \$seq_region_end, \$genotypes, \$population_id);
  while($sth->fetch()) {
    #only print genotypes without parents genotyped
    if (!exists $siblings->{$population_id . '-' . $sample_id}){ #necessary to use the population_id
      $self->_store_genotype(\%sample_information, \%alleles_variation, $sample_id, $seq_region_start, $genotypes, $population_id, $slice);
      $previous_seq_region_id = $seq_region_id;
    }
  }
  
  $sth->finish();
  return {
    'alleles_variation' => \%alleles_variation,
    'sample_information' => \%sample_information,
    'slice' => $slice,
  }; 
}

sub _objs_from_sth_vcf {
  my $self = shift;
  my $gts = shift;
  my $slice = shift;
  my $pops = shift;
  my $pos2name = shift;
  
  my (%alleles_variation, %sample_information);
  
  # create an artificial numerical ID for samples
  my $current_sample_id = 1;
  my %sample_ids;
  
  foreach my $snp_start(keys %$gts) {
    foreach my $sample(sort keys %{$gts->{$snp_start}}) {
      next unless $gts->{$snp_start}->{$sample};
      
      my @gt = split(/\||\//, $gts->{$snp_start}->{$sample});

      next unless scalar @gt == 2;
      next unless grep {defined($_)} @gt;
      
      # get sample ID
      if(!exists($sample_ids{$sample})) {
        $sample_ids{$sample} = $current_sample_id++;
      }
      my $sample_id = $sample_ids{$sample};
      
      foreach my $pop_id(keys %{$pops->{$sample}}) {
        $alleles_variation{$snp_start}->{$pop_id}->{$gt[0]}++;
        $alleles_variation{$snp_start}->{$pop_id}->{$gt[1]}++;
    
        $sample_information{$pop_id}->{$snp_start}->{$sample_id}->{allele_1} = $gt[0];
        $sample_information{$pop_id}->{$snp_start}->{$sample_id}->{allele_2} = $gt[1];
      }
    }
  }
  return {
    'alleles_variation' => \%alleles_variation,
    'sample_information' => \%sample_information,
    'slice' => $slice,
    'pos2name' => $pos2name,
  }; 
}

sub _ld_calc {
  my $self = shift;
  my $genotype_hashes = shift; 
  my $use_vcf = $self->db->use_vcf;
  my $alleles_variation = {};
  my $sample_information = {};
  my $pos2name = {};
  my $pos2vf = {};
  my @slices;
  my $vfa;

  # what we're doing here is copying data from component hashes to a merged hash
  # this is because we can receive more than one genotype hash
  foreach my $genotype_hash (@$genotype_hashes) { 

    # this contains allele frequency data
    my $next_alleles_variation = $genotype_hash->{'alleles_variation'};
    foreach my $snp_start (keys %$next_alleles_variation) {
      foreach my $pop_id (keys %{$next_alleles_variation->{$snp_start}}) {
        my $hash = $next_alleles_variation->{$snp_start}->{$pop_id};
        $alleles_variation->{$snp_start}->{$pop_id} = $hash; 
      }
    }

    # this contains the actual genotypes
    my $next_sample_information = $genotype_hash->{'sample_information'};
    foreach my $pop_id (keys %$next_sample_information) {
      foreach my $snp_start (keys %{$next_sample_information->{$pop_id}}) {
        my $hash = $next_sample_information->{$pop_id}->{$snp_start};
        $sample_information->{$pop_id}->{$snp_start} = $hash;
      }
    }

    # we'll store the slices and send them to the container
    # this way the container can lazy-load the VariationFeature objects
    push @slices, $genotype_hash->{'slice'};
    
    # the pos2name hash maps positions to variant names
    # the VCF API sends these through as part of its genotype fetch
    if(my $next_pos2name = $genotype_hash->{'pos2name'}) {
      $pos2name->{$_} = $next_pos2name->{$_} for keys %$next_pos2name;
    }
    # But the database API has no such feature
    # this means we need to look up the VariationFeatures here.
    # But since we can do that here, we can then also send pos2vf through
    # to the container so no lazy-loading needs to be done
    else {
      my $slice = $genotype_hash->{'slice'};
      $vfa ||= $self->db->get_VariationFeatureAdaptor();
      my $vfs = $vfa->fetch_all_by_Slice($slice);
      my $region_Slice = $slice->seq_region_Slice();

      # the sort here "favours" variants from dbSNP since it is assumed dbSNP will have source_id = 1
      # otherwise we'd get co-located COSMIC IDs overwriting the dbSNP ones
      foreach my $vf(map {$_->transfer($region_Slice)} sort {$b->{_source_id} <=> $a->{_source_id}} @{$vfs}) {
        $pos2name->{$vf->start} = $vf->variation_name;
        $pos2vf->{$vf->start} = $vf;
      }
    }
  }
  
  my %_pop_ids;

  my $previous_seq_region_id = 0;
  
  #we have to print the variations
  my (%in_files, %in_file_names);

  foreach my $snp_start (sort{$a<=>$b} keys %$alleles_variation){
    foreach my $population (keys %{$alleles_variation->{$snp_start}}){
      my $fh;
	  
      # create file handles in hash
      if(!defined($in_files{$population})) {
        $fh = FileHandle->new;
        my $f_name = $self->temp_path."/".sprintf( "ld%08x%08x%08x%08x", $population, $$, time, rand( 0x7fffffff));		
        $in_file_names{$population} = $f_name;
        $fh->open(">".$f_name.".in");
        $in_files{$population} = $fh;
      }
	  
      else {
        $fh = $in_files{$population};
      }
	  
      #if the variation has 2 alleles, print all the genotypes to the file
      if (keys %{$alleles_variation->{$snp_start}{$population}} == 2){
        $self->_convert_genotype($alleles_variation->{$snp_start}{$population},$sample_information->{$population}{$snp_start});
        foreach my $sample_id (keys %{$sample_information->{$population}{$snp_start}}){
          print $fh join("\t",$previous_seq_region_id, $snp_start, $snp_start, $population, $sample_id, $sample_information->{$population}{$snp_start}{$sample_id}{genotype})."\n" || warn $!;

        }
      }
    }
  }
  
  # close file handles and check file sizes
  foreach my $key(keys %in_files) {
    my $f = $in_files{$key};
    $f->close;
    my $file = $in_file_names{$key} . '.in';
    if (-z $file) { # file is empty
      unlink($in_file_names{$key}.'.in');
      delete $in_file_names{$key};
      delete $in_files{$key};
    }
  }

  my @cmd = qw(calc_genotypes);

  #open the pipe between processes if the binary file exists in the PATH
  my $bin = $self->executable;
  if( ! $bin ) {
    warning("Binary file calc_genotypes not found. Please, read the ensembl-variation/C_code/README.txt file if you want to use LD calculation\n");
    goto OUT;
  }
 
  # run LD binary
  `$bin <$_\.in >$_\.out` for values %in_file_names;

  # now create the container from the output of the LD binary
  my %feature_container = ();
 
  foreach my $file (values %in_file_names) {
    open OUT, "$file.out";
    while(<OUT>){
      my %ld_values = ();
	  
      #     936	965891	164284	166818	0.628094	0.999996	120 
      #get the ouput into the hashes
      chomp;
      my ($population_id,$ld_region_id,$ld_region_start,$ld_region_end,$r2,$d_prime,$sample_count) = split /\s/;
      # skip entries unrelated to selected vf if doing fetch_all_by_VariationFeature
      if (defined($self->{_vf_pos})) {
        next unless $ld_region_start == $self->{_vf_pos} || $ld_region_end == $self->{_vf_pos};
      }
      # skip entries for pairwise computation that don't match input variation feature loactions
      if (defined $self->{_pairwise}) {
        next unless ($self->{_pairwise}->{$ld_region_start} && $self->{_pairwise}->{$ld_region_end});
      } 

      $ld_values{'d_prime'} = $d_prime;
      $ld_values{'r2'} = $r2;
      $ld_values{'sample_count'} = $sample_count;
	  
      if (!defined $pos2name->{$ld_region_start} || !defined $pos2name->{$ld_region_end}){
        next; #problem to fix in the compressed genotype table: some of the positions seem to be wrong
      }
    
      $feature_container{$ld_region_start . '-' . $ld_region_end}->{$population_id} = \%ld_values;
	  
      $_pop_ids{$population_id} = 1;	  
    }
    close OUT || die "Could not close filehandle: $!\n";
	
    unlink( "$file.in" );
    unlink( "$file.out" );
  }
  OUT:
  my $t = Bio::EnsEMBL::Variation::LDFeatureContainer->new(
    '-ldContainer'=> \%feature_container,
    '-name' => '',
    '-pos2name' => $pos2name,
    '-pos2vf' => $pos2vf,
    '-slices' => \@slices,
  );

  $t->{'_pop_ids'} =\%_pop_ids;
  delete $self->{_pairwise};
  
  return $t;      
}

1;

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

=head1 NAME

Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer

=head1 SYNOPSIS

=head1 DESCRIPTION

A container for TranscriptHaplotype objects

=cut

package Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use Bio::EnsEMBL::Variation::CDSHaplotype;
use Bio::EnsEMBL::Variation::ProteinHaplotype;
use Bio::EnsEMBL::Variation::TranscriptDiplotype;

use Digest::MD5 qw(md5_hex);
use Scalar::Util qw(weaken);



=head2 new

  Arg [-TRANSCRIPT]: Bio::EnsEMBL::Transcript
  Arg [-GENOTYPES]:  arrayref of Bio::EnsEMBL::Variation::SampleGenotypeFeature (non-reference genotypes)
  Arg [-SAMPLES]:    arrayref of Bio::EnsEMBL::Variation::Sample (all samples)
  Arg [-DB]:         Bio::EnsEMBL::Variation::DBSQL::DBAdaptor

  Example    : my $thc = Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer->new(
                  -TRANSCRIPT => $transcript,
                  -GENOTYPES  => $genotypes,
                  -SAMPLES    => $samples,
                  -DB         => $db
               );

  Description: Constructor.  Instantiates a new TranscriptHaplotypeContainer object.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::TranscriptHaplotypeContainer
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($transcript, $gts, $samples, $db) = rearrange([qw(TRANSCRIPT GENOTYPES SAMPLES DB)], @_);
  
  # check what we've been given looks sensible
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'Transcript');

  # only works for coding transcripts ATM
  return undef unless $transcript->translation;

  assert_ref($gts, 'ARRAY', 'Genotypes listref');
  assert_ref($gts->[0], 'Bio::EnsEMBL::Variation::SampleGenotypeFeature', 'First member of genotypes listref') if scalar @$gts;

  if($samples) {
    assert_ref($samples, 'ARRAY', 'Samples listref');
    assert_ref($samples->[0], 'Bio::EnsEMBL::Variation::Sample', 'First member of samples listref');
  }
  
  my $self = {
    _transcript => $transcript,
    _sample_genotype_features => $gts || [],
    _db => $db,
    _samples => $samples,
    transcript_id => $transcript->stable_id,
  };
  bless($self, $class);
  
  $self->_init();
  
  return $self;
}


=head2 transcript

  Example    : my $tr = $thc->transcript()
  Description: Getter/setter for the Transcript object used as the reference
               for this container
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub transcript {
  my $self = shift;
  $self->{_transcript} = shift if @_;
  return $self->{_transcript};
}


=head2 get_all_Samples

  Example    : my $samples = $thc->get_all_Samples()
  Description: Get list of samples in this container
  Returntype : listref of Bio::EnsEMBL::Variation::Sample
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Samples {
  my $self = shift;

  if(!$self->{_samples}) {
    $self->{_samples} = [values %{{map {$_->sample->name => $_->sample()} @{$self->get_all_SampleGenotypeFeatures}}}];
  }

  return $self->{_samples};
}


=head2 get_all_Populations

  Example    : my $pops = $thc->get_all_Populations()
  Description: Get list of populations in this container
  Returntype : listref of Bio::EnsEMBL::Variation::Population
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Populations {
  my $self = shift;

  if(!exists($self->{_population_objects})) {
    $self->{_population_objects} = [values %{{map {$_->dbID => $_} map {@{$_->get_all_Populations}} @{$self->get_all_Samples}}}];
  }

  return $self->{_population_objects};
}


=head2 get_TranscriptHaplotype_by_name

  Example    : my $th = $thc->get_TranscriptHaplotype_by_name()
  Description: Fetch a TranscriptHaplotype by name
  Returntype : Bio::EnsEMBL::Variation::TranscriptHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_TranscriptHaplotype_by_name {
  my $self = shift;
  my $name = shift;
  my @matched = grep {$_->name eq $name} @{$self->get_all_TranscriptHaplotypes};
  return scalar @matched ? $matched[0] : undef;
}


=head2 get_all_TranscriptHaplotypes

  Example    : my @ths = @{$thc->get_all_CDSHaplotypes()}
  Description: Get all CDS and protein haplotypes for this container
  Returntype : arrayref of Bio::EnsEMBL::Variation::TranscriptHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_TranscriptHaplotypes {
  return [@{$_[0]->get_all_CDSHaplotypes}, @{$_[0]->get_all_ProteinHaplotypes}];
}


=head2 get_all_CDSHaplotypes

  Example    : my @chs = @{$thc->get_all_CDSHaplotypes()}
  Description: Get all CDS haplotypes for this container
  Returntype : arrayref of Bio::EnsEMBL::Variation::CDSHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_CDSHaplotypes {
  return [values %{$_[0]->{cds_haplotypes}}];
}


=head2 get_all_ProteinHaplotypes

  Example    : my @phs = @{$thc->get_all_ProteinHaplotypes()}
  Description: Get all protein haplotypes for this container
  Returntype : arrayref of Bio::EnsEMBL::Variation::ProteinHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_ProteinHaplotypes {
  return [values %{$_[0]->{protein_haplotypes}}];
}


=head2 get_all_TranscriptHaplotypes_by_Sample

  Arg[1]     : Bio::EnsEMBL::Variation::Sample $sample
  Example    : my @ths = @{$thc->get_all_TranscriptHaplotypes_by_Sample()}
  Description: Get all CDS and protein haplotypes for a specific sample
  Returntype : arrayref of Bio::EnsEMBL::Variation::TranscriptHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_TranscriptHaplotypes_by_Sample {
  my $self = shift;
  my $sample = shift;
  return [grep {$_->{samples}->{$sample->name}} @{$self->get_all_TranscriptHaplotypes}];
}


=head2 get_all_CDSHaplotypes_by_Sample

  Arg[1]     : Bio::EnsEMBL::Variation::Sample $sample
  Example    : my @chs = @{$thc->get_all_CDSHaplotypes_by_Sample()}
  Description: Get all CDS haplotypes for a specified sample
  Returntype : arrayref of Bio::EnsEMBL::Variation::CDSHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_CDSHaplotypes_by_Sample {
  my $self = shift;
  my $sample = shift;
  return [grep {$_->{samples}->{$sample->name}} @{$self->get_all_CDSHaplotypes}];
}


=head2 get_all_ProteinHaplotypes_by_Sample

  Arg[1]     : Bio::EnsEMBL::Variation::Sample $sample
  Example    : my @phs = @{$thc->get_all_ProteinHaplotypes_by_Sample()}
  Description: Get all protein haplotypes for a specified sample
  Returntype : arrayref of Bio::EnsEMBL::Variation::ProteinHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_ProteinHaplotypes_by_Sample {
  my $self = shift;
  my $sample = shift;
  return [grep {$_->{samples}->{$sample->name}} @{$self->get_all_ProteinHaplotypes}];
}


=head2 get_all_most_frequent_CDSHaplotypes

  Example    : my @chs = @{$thc->get_all_most_frequent_CDSHaplotypes()}
  Description: Get all CDS haplotypes that have an observed count equal
               to the maximum observed count
  Returntype : arrayref of Bio::EnsEMBL::Variation::CDSHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_most_frequent_CDSHaplotypes {
  my $self = shift;
  
  if(!exists($self->{_most_frequent_cds})) {
    $self->{_most_frequent_cds} = $self->_get_most_frequent($self->get_all_CDSHaplotypes);
  }
  return $self->{_most_frequent_cds};
}


=head2 get_all_most_frequent_ProteinHaplotypes

  Example    : my @phs = @{$thc->get_all_most_frequent_ProteinHaplotypes()}
  Description: Get all protein haplotypes that have an observed count equal
               to the maximum observed count
  Returntype : arrayref of Bio::EnsEMBL::ProteinHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_most_frequent_ProteinHaplotypes {
  my $self = shift;
  
  if(!exists($self->{_most_frequent_protein})) {
    $self->{_most_frequent_protein} = $self->_get_most_frequent($self->get_all_ProteinHaplotypes);
  }
  return $self->{_most_frequent_protein};
}


=head2 get_all_CDSDiplotypes

  Example    : my @cds = @{$thc->get_all_CDSDiplotypes()}
  Description: Get all uniquely observed CDS diplotypes - a diplotype is the
               pair of haplotypes observed in a given sample/individual
  Returntype : arrayref of Bio::EnsEMBL::Variation::TranscriptDiplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_CDSDiplotypes {
  my $self = shift;
  return $self->_generic_get_all_diplotypes('CDS', @_);
}


=head2 get_all_ProteinDiplotypes

  Example    : my @pds = @{$thc->get_all_ProteinDiplotypes()}
  Description: Get all uniquely observed protein diplotypes - a diplotype is the
               pair of haplotypes observed in a given sample/individual
  Returntype : arrayref of Bio::EnsEMBL::Variation::TranscriptDiplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_ProteinDiplotypes {
  my $self = shift;
  return $self->_generic_get_all_diplotypes('Protein', @_);
}

=head2 total_haplotype_count

  Example    : my $count = $thc->total_haplotype_count
  Description: Get total observed haplotype count
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub total_haplotype_count {
  my $self = shift;
  
  if(!exists($self->{total_haplotype_count})) {
    my $count = 0;
    $count += $_ for values %{$self->{_counts}};
    
    # divide by 2 only if we have protein haplotypes
    $count /= 2 if scalar keys %{$self->{protein_haplotypes}};
    
    $self->{total_haplotype_count} = $count;
  }
  
  return $self->{total_haplotype_count};
}


=head2 total_population_counts

  Example    : my %counts = %{$thc->total_population_counts()}
  Description: Get observed haplotype counts by population name
  Returntype : hashref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub total_population_counts {
  my $self = shift;
  
  if(!exists($self->{total_population_counts})) {
    my $counts = {};
    
    my $hash = $self->_get_sample_population_hash();
    
    foreach my $sample(keys %$hash) {
      $counts->{$_} += (scalar keys %{$self->{protein_haplotypes}} ? 2 : 1) for keys %{$hash->{$sample}};
    }
    
    $self->{total_population_counts} = $counts;
  }
  
  return $self->{total_population_counts};
}


=head2 get_all_SampleGenotypeFeatures

  Example    : my @gts = @{$thc->get_all_SampleGenotypeFeatures()}
  Description: Get all SampleGenotypeFeature objects associated with
               this container
  Returntype : arrayref of Bio::EnsEMBL::Variation::SampleGenotypeFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_SampleGenotypeFeatures {
  my $self = shift;
  return $self->{_sample_genotype_features};
}


=head2 db

  Example    : my $db = $thc->db()
  Description: Getter/setter for the DB object for this container
  Returntype : Bio::EnsEMBL::DBAdaptor
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub db {
  my $self = shift;
  $self->{_db} = shift if @_;
  return $self->{_db};
}



###################
# PRIVATE METHODS #
###################

## Get/set all the unique VariationFeatures from the SampleGenotypeFeature objects
sub _variation_features {
  my $self = shift;
  my $vfs = shift;
  
  if($vfs) {
    $self->{_variation_features} = $vfs;
  }
  elsif(!exists($self->{_variation_features})) {
    $self->{_variation_features} = [map {$_->variation_feature} @{$self->get_all_SampleGenotypeFeatures}];
  }
  
  return $self->{_variation_features};
}

## Fetch a haplotype by its hex, dont care about type
sub _get_TranscriptHaplotype_by_hex {
  return $_[0]->{cds_haplotypes}->{$_[1]} || $_[0]->{protein_haplotypes}->{$_[1]} || undef;
}

## find the most frequent haplotype amongst a given arrayref
sub _get_most_frequent {
  my $self = shift;
  my $haplos = shift;
  
  # get hexes of all the haplotypes
  my @hexes = map {$_->_hex} @$haplos;
  
  # find the maximum count from this list of hexes
  my $max_count = (sort {$a <=> $b} map {$self->{_counts}->{$_} || 0} @hexes)[-1];
  
  return [grep {$self->{_counts}->{$_->_hex} == $max_count} @$haplos];
}

## Generates a two-level hash giving the populations for each sample
sub _get_sample_population_hash {
  my $self = shift;
  
  if(!exists($self->{_sample_population_hash})) {
    my $hash = {};

    foreach my $sample(@{$self->get_all_Samples}) {
      $hash->{$sample->name}->{$_->name} = 1 for @{$sample->get_all_Populations};
    }
    
    $self->{_sample_population_hash} = $hash;
  }
  
  return $self->{_sample_population_hash};
}

## generic sub used by get_all_ProteinDiplotypes and get_all_CDSDiplotypes
## there's nothing particular to CDS or Protein diplotypes
## so all logic is shared
sub _generic_get_all_diplotypes {
  my ($self, $type) = @_;

  my $obj_key = '_'.$type.'_diplotypes';

  if(!exists($self->{$obj_key})) {

    # what we want is a hash of haplotype pairs
    # keyed on individual
    my %by_sample = ();

    my $method = 'get_all_'.$type.'Haplotypes';

    foreach my $ht(@{$self->$method}) {
      foreach my $sample(keys %{$ht->{samples}}) {

        # add for each time it is observed, i.e. add twice if it is homozygous
        push @{$by_sample{$sample}}, $ht for 1..$ht->{samples}->{$sample};
      }
    }

    # now do some aggregation
    my %diplotypes;

    foreach my $sample(keys %by_sample) {

      # sort them by hex so we're always consistent in the ordering
      # enables uniquification
      my @haplos = sort {$a->_hex cmp $b->_hex} @{$by_sample{$sample}};

      # create a unique key by joining the hexes
      my $hex_pair = join("_", map {$_->_hex} @haplos);

      # create a hashref that will become an object
      $diplotypes{$hex_pair} ||= {
        -haplotypes => \@haplos
      };

      # add samples as we come across them
      push @{$diplotypes{$hex_pair}->{-samples}}, $sample;
    }

    # now convert to objects
    $self->{$obj_key} = [map {Bio::EnsEMBL::Variation::TranscriptDiplotype->new(%$_)} values %diplotypes];
  }

  return $self->{$obj_key};
}

## we need this to get diplotype frequencies
sub _total_sample_count {
  my $self = shift;

  if(!exists($self->{_total_sample_count})) {
    $self->{_total_sample_count} = scalar @{$self->get_all_Samples};
  }

  return $self->{_total_sample_count};
}

## gets total population counts of diplotypes
## since there will be one diplotype per sample, this basically
## just returns a hashref of sample counts per population
sub _total_diplotype_population_counts {
  my $self = shift;
  
  if(!exists($self->{_total_diplotype_population_counts})) {
    my $counts = {};
    
    my $hash = $self->_get_sample_population_hash();
    
    foreach my $sample(keys %$hash) {
      $counts->{$_}++ for keys %{$hash->{$sample}};
    }
    
    $self->{_total_diplotype_population_counts} = $counts;
  }
  
  return $self->{_total_diplotype_population_counts};
}

## Creates all the TranscriptHaplotype objects for this container
sub _init {
  my $self = shift;
  
  my $tr = $self->transcript;
  my $vfs = $self->_variation_features;

  my $is_protein_coding = $tr->translation ? 1 : 0;

  # cache reference sequences on transcript object
  $tr->{cds} = $tr->{_variation_effect_feature_cache}->{translateable_seq} ||= $tr->translateable_seq;

  # find out if the transcript has a conventional stop codon
  my $last_codon = substr($tr->{cds}, -3, 3);
  my $has_stop_codon = ($last_codon eq 'TAA' || $last_codon eq 'TAG' || $last_codon eq 'TGA') ? 1 : 0;

  # only append '*' to the sequence if the CDS has a stop codon
  # otherwise it looks like * has disappeared every time from transcripts with no stop codon
  $tr->{protein} = $is_protein_coding ? ($tr->{_variation_effect_feature_cache}->{peptide} ||= $tr->translation->seq).($has_stop_codon ? '*' : '') : undef;
  
  # get mappings of all variation feature coordinates
  my $mappings = $self->_get_mappings;
  
  # remove any variation features that didn't get a mapping
  my @new_vars = grep {$_->{_cds_mapping}} @$vfs;
  $self->_variation_features(\@new_vars);
  
  # group genotypes by sample and sort
  my %by_sample;
  push @{$by_sample{$_->sample->name}}, $_ for @{$self->get_all_SampleGenotypeFeatures};
  
  # do the transcript magic
  foreach my $sample(keys %by_sample) {
    my $obj = $self->_mutate_sequences($self->_filter_and_sort_genotypes($by_sample{$sample}), $sample);
    
    # store unique alt seqs so we only align each once
    foreach my $type(qw(cds protein)) {

      next if $type eq 'protein' && !$is_protein_coding;
      
      # iterate over haplotypes
      # usually there will be 2, though for e.g. male samples on an chrX transcript there will be 1
      foreach my $i(0..$#{$obj->{$type}}) {
        
        my $seq = $obj->{$type}->[$i];
        my $hex = md5_hex($seq);
        
        # try to fetch 
        my $haplotype = $self->_get_TranscriptHaplotype_by_hex($hex);
        
        # if it doesn't exist yet, create new object
        if(!$haplotype) {          
          my %hash = (
            -container => $self,
            -type      => $type,
            -seq       => $seq,
            -hex       => $hex,
            -indel     => $obj->{indel}->[$i]
          );
          
          if($type eq 'cds') {
            $haplotype = Bio::EnsEMBL::Variation::CDSHaplotype->new(%hash);
          }
          else {
            $haplotype = Bio::EnsEMBL::Variation::ProteinHaplotype->new(%hash);
          }
          
          $self->{$type.'_haplotypes'}->{$hex} = $haplotype;
          
          weaken($haplotype->{_container});
        }
        
        # increment counts
        $haplotype->{samples}->{$sample}++;
        $self->{_counts}->{$hex}++;
        
        # store mapping between hashes of cds and protein
        if($is_protein_coding) {
          my $other_type = $type eq 'cds' ? 'protein' : 'cds';
          my $other_hex  = md5_hex($obj->{$other_type}->[$i]);
          
          $haplotype->{other_hexes}->{$other_hex} = 1;
        }
      }
    }
  }

  $self->_add_reference_haplotypes();
  
  # clear cached translated sequences
  delete($self->{_translations});
}

sub _add_reference_haplotypes {
  my $self = shift;

  # get list of samples that have alt GTs
  my %have_alt_gts = map {$_->sample->name => 1} @{$self->get_all_SampleGenotypeFeatures};

  # now get the remaining samples and assign them reference status
  my @ref_samples = grep {!$have_alt_gts{$_->name}} @{$self->get_all_Samples};

  if(scalar @ref_samples) {
    my $tr = $self->transcript;

    my $is_protein_coding = $tr->translation ? 1 : 0;

    # get/create reference CDS haplotype
    my $ref_cds_seq = $tr->{cds};
    my $ref_cds_hex = md5_hex($ref_cds_seq);

    my $ref_cds_haplotype = $self->_get_TranscriptHaplotype_by_hex($ref_cds_hex);

    if(!$ref_cds_haplotype) {
      $ref_cds_haplotype = Bio::EnsEMBL::Variation::CDSHaplotype->new(
        -container => $self,
        -type      => 'cds',
        -seq       => $ref_cds_seq,
        -hex       => $ref_cds_hex,
        -indel     => 0
      );

      $self->{'cds_haplotypes'}->{$ref_cds_hex} = $ref_cds_haplotype;
      
      weaken($ref_cds_haplotype->{_container});
    }

    # get/create reference protein haplotype
    my $ref_protein_haplotype;
    my $ref_protein_hex;

    if($is_protein_coding) {
      my $ref_protein_seq = $tr->{protein};
      $ref_protein_hex = md5_hex($ref_protein_seq);

      $ref_protein_haplotype = $self->_get_TranscriptHaplotype_by_hex($ref_protein_hex);

      if(!$ref_protein_haplotype) {
        $ref_protein_haplotype = Bio::EnsEMBL::Variation::ProteinHaplotype->new(
          -container => $self,
          -type      => 'protein',
          -seq       => $ref_protein_seq,
          -hex       => $ref_protein_hex,
          -indel     => 0,
        );

        $self->{'protein_haplotypes'}->{$ref_protein_haplotype} = $ref_protein_haplotype;
        
        weaken($ref_protein_haplotype->{_container});
      }

      # update other haplotype hex
      $ref_cds_haplotype->{other_hexes}->{$ref_protein_hex} = 1;
      $ref_protein_haplotype->{other_hexes}->{$ref_cds_hex} = 1;
    }

    # increment counts
    my $sample_ploidy  = $self->_sample_ploidy();
    my $default_ploidy = $self->_default_ploidy();

    for(@ref_samples) {
      my $ploidy = $sample_ploidy->{$_->name} || $default_ploidy;

      $ref_cds_haplotype->{samples}->{$_->name}     += $ploidy;
      $self->{_counts}->{$ref_cds_hex}              += $ploidy;

      if($is_protein_coding) {
        $ref_protein_haplotype->{samples}->{$_->name} += $ploidy;
        $self->{_counts}->{$ref_protein_hex}          += $ploidy;
      }
    }
  }
}

## get the ploidy of each sample
## this will vary e.g. for different genders in an X chromosome region
## we do this by fetching ALL genotypes for one variant
sub _sample_ploidy {
  my $self = shift;

  if(!exists($self->{_sample_ploidy})) {
    my $gts = $self->get_all_SampleGenotypeFeatures;

    if($gts && scalar @$gts) {
      if(my $v = $gts->[0]->variation) {        
        %{$self->{_sample_ploidy}} =
          map {$_->sample->name => scalar @{$_->genotype}}
          @{$v->get_all_SampleGenotypes};
      }
    }
    
    if(!$self->{_sample_ploidy}) {
      my $default_ploidy = $self->_default_ploidy();
      %{$self->{_sample_ploidy}} = map {$_->name => $default_ploidy} @{$self->get_all_Samples};
    }
  }
  
  return $self->{_sample_ploidy};
}

## get default ploidy
## fallback for if we can't get per-sample ploidy
sub _default_ploidy {
  my $self = shift;

  if(!exists($self->{_ploidy})) {
    if(my $db = $self->db) {
      if(my $va = $db->get_VariationAdaptor) {
        $self->{_ploidy} = $va->ploidy;
      }
    }

    $self->{_ploidy} = 2 unless defined($self->{_ploidy});
  }

  return $self->{_ploidy};
}

## Uses a transcript's mapper to get genomic->CDS coord mappings for each
## VariationFeature, with the results stored on $vf
sub _get_mappings {
  my $self = shift;
  
  my $tr = $self->transcript();
  my $tr_strand = $tr->strand;
  my $mapper = $tr->{_variation_effect_feature_cache}->{mapper} ||= $tr->get_TranscriptMapper;
  
  foreach my $vf(@{$self->_variation_features}) {
    next if $vf->{_cds_mapping} || $vf->{_cds_mapping_failed};

    my ($vf_start, $vf_end) = $vf->{slice} ? ($vf->seq_region_start, $vf->seq_region_end) : ($vf->start, $vf->end);
    
    my @mapped = $mapper->genomic2cds($vf_start, $vf_end, $tr_strand);

    # clean mapping
    if(scalar @mapped == 1) {
      if($mapped[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
        $vf->{_cds_mapping_failed} = 'gap';
      }
      else {
        $vf->{_cds_mapping} = $mapped[0];
      }
    }

    else {

      # complex mapping, we can't deal with this ATM
      if(grep {$_->isa('Bio::EnsEMBL::Mapper::Coordinate')} @mapped) {
        warn "WARNING: genomic coord ".$vf_start."-".$vf_end." possibly maps across coding/non-coding boundary in ".$tr->stable_id."\n";

        $vf->{_cds_mapping_failed} = 'complex';
      }

      else {
        $vf->{_cds_mapping_failed} ||= 'multiple_gaps';
      }
    }
  }
}

## this subroutine filters out genotypes with no mappings
## and sorts them into reverse order ready for applying to the sequence
sub _filter_and_sort_genotypes {
  my $self = shift;
  my $gts = shift;

  return [
    map  { $_->[0] }
    sort { $b->[1] <=> $a->[1] }
    map  { [$_, $_->variation_feature->{_cds_mapping}->start] }
    grep { $_->variation_feature->{_cds_mapping} }
    @$gts
  ];
}

## Applies a set of SampleGenotypeFeatures to the transcript sequence
## Returns a hashref that contains the info necessary to construct a
## TranscriptHaplotype object
sub _mutate_sequences {
  my $self = shift;
  my $gts = shift;
  my $sample_name = shift;
  
  my $fingerprint = $self->_fingerprint_gts($gts);
  
  if(!exists($self->{_mutated_by_fingerprint}) || !exists($self->{_mutated_by_fingerprint}->{$fingerprint})) {
  
    my $tr = $self->transcript;
    my $tr_strand = $tr->strand;
  
    my $codon_table;
  
    if($tr->{_variation_effect_feature_cache} && $tr->{_variation_effect_feature_cache}->{codon_table}) {
      $codon_table = $tr->{_variation_effect_feature_cache}->{codon_table};
    }
    else {
      my $attrib = $tr->slice->get_all_Attributes('codon_table')->[0];
      $codon_table = $attrib ? $attrib->value : 1;
    }
  
    my $mutated = {};

    # get sample ploidy
    my $ploidy = $self->_sample_ploidy()->{$sample_name} || $self->_default_ploidy();
  
    for my $hap(0..($ploidy - 1)) {
      my $seq = $tr->{_variation_effect_feature_cache}->{translateable_seq} ||= $tr->translateable_seq;
      my $ref_seq = $seq;
    
      # flag if indel observed, we need to align seqs later
      my $indel = 0;
    
      # iterate through, they are already sorted into reverse order by sequence mapping
      foreach my $gt(@$gts) {
        my $vf = $gt->variation_feature;
        my ($s, $e) = ($vf->{start}, $vf->{end});
      
        my $mapping = $vf->{_cds_mapping};
      
        my $genotype = $gt->genotype->[$hap];
      
        # remove del characters
        $genotype =~ s/\-//g;
      
        # reverse complement sequence?
        reverse_comp(\$genotype) if $tr_strand ne $vf->strand;
      
        my $replace_len = ($mapping->end - $mapping->start) + 1;
        $indel = 1 if length($genotype) != $replace_len;
      
        substr($seq, $mapping->start - 1, $replace_len, $genotype);
      }
    
      # now translate
      my $protein = $self->_get_translation($seq, $codon_table) if $tr->translation;
    
      # remove anything beyond a stop
      $protein =~ s/\*.+/\*/;

      # remove stop
      # $protein =~ s/\*$//;
    
      push @{$mutated->{cds}}, $seq;
      push @{$mutated->{protein}}, $protein if $protein;
      push @{$mutated->{indel}}, $indel;
    }
  
    $self->{_mutated_by_fingerprint}->{$fingerprint} = $mutated;
  }
  
  return $self->{_mutated_by_fingerprint}->{$fingerprint};
}

## gets an md5 hex for an samples VFs/genotypes combination
## this means we don't have to run _mutate_sequences more than
## once for essentially the same sequence
sub _fingerprint_gts {
  my $self = shift;
  my $gts = shift;
  
  my @raw;
  
  foreach my $gt(@$gts) {
    my $s = $gt->genotype_string;
    my $vf = $gt->variation_feature;
    
    my $v = $vf->dbID || $vf->{start}.'-'.$vf->{end};
    
    push @raw, sprintf('%s:%s', $v, $s);
  }
  
  return md5_hex(join('_', @raw));
}

## Gets a translation - these are cached on $self by a hex of the seq to avoid
## doing the translation more than once for each distinct CDS sequence
sub _get_translation {
  my $self = shift;
  my $seq = shift;
  my $codon_table = shift;
  
  my $hex = md5_hex($seq);
  
  if(!exists($self->{_translations}) || !exists($self->{_translations}->{$hex})) {
    $self->{_translations}->{$hex} = Bio::Seq->new(-seq => $seq)->translate(undef, undef, undef, $codon_table)->seq;
  }
  
  return $self->{_translations}->{$hex};
}

## cache ProteinFunctionPredictionMatrices
## used by ProteinHaplotypes to get SIFT/PolyPhen scores for diffs
sub _get_prediction_matrices {
  my $self = shift;
  
  if(!exists($self->{_prediction_matrices})) {
    my $matrices = {};
    
    my $tr = $self->transcript;
    
    # retrive from DB
    if(my $db = $self->db) {
      my $pfpma = $db->get_ProteinFunctionPredictionMatrixAdaptor();
      my $peptide = $tr->{_variation_effect_feature_cache}->{peptide} ||= $tr->translation->seq;
      my $tr_md5 = md5_hex($peptide);
      
      foreach my $tool(qw(sift polyphen)) {
        my $method = sprintf('fetch_%s_predictions_by_translation_md5', $tool);
        
        $matrices->{$tool} = $pfpma->$method($tr_md5);
      }
    }
    
    # or can retrieve from cache if transcript has come from VEP cache
    elsif($matrices = $tr->{_variation_effect_feature_cache}->{protein_function_predictions}) {
      $matrices->{polyphen} = $matrices->{polyphen_humvar} if $matrices->{polyphen_humvar};
      delete $matrices->{polyphen_humdiv} if $matrices->{polyphen_humdiv};
    }
    
    $self->{_prediction_matrices} = $matrices;
  }
  
  return $self->{_prediction_matrices};
}

## Generates TranscriptVariation objects representing consequences.
## Only used when our VariationFeatures have been created from e.g. a VCF file
sub _get_unique_VariationFeatures_with_consequences {
  my $self = shift;
  
  my $tr = $self->transcript;
  my $vfs = $self->_variation_features;
  
  # uniquify the VFs - this is because there will be multiple VFs with the same
  # coords and alleles due to the way the VEP's VCF parser works
  my %unique_vfs;
  
  foreach my $vf(@$vfs) {
    my @alleles = split /\//, $vf->{allele_string};
    my $ref_allele = shift @alleles;
    
    my ($s, $e) = ($vf->{start}, $vf->{end});
    
    my $key = join('-', $s, $e, $ref_allele);
    
    if(!defined($unique_vfs{$key})) {
      my $vf_copy = { %$vf };
      bless $vf_copy, ref($vf);
      delete $vf_copy->{$_} for qw(_line phased sample genotype non_variant hom_ref);
      $vf_copy->{ref_allele} = $ref_allele;
      $unique_vfs{$key} = $vf_copy;
    }
    
    $unique_vfs{$key}->{alleles}->{$_} = 1 for @alleles;
  }
  
  # we only need a fake TVA
  my $tva = Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor->new_fake;
  
  # merge VF alleles and get consequences
  foreach my $vf(values %unique_vfs) {
    $vf->{allele_string} = join("/", $vf->{ref_allele}, keys %{$vf->{alleles}});
    $vf->{non_variant} = 1 unless scalar keys %{$vf->{alleles}};
    delete $vf->{$_} for qw(alleles ref_allele);
    
    my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
      -transcript        => $tr,
      -variation_feature => $vf,
      -adaptor           => $tva,
      -no_ref_check      => 1,
      -no_transfer       => 1
    );
    
    # prefetching stuff here prevents doing loads at the
    # end and makes progress reporting more useful
    $tv->_prefetch_for_vep;
    
    $vf->add_TranscriptVariation($tv);
  }
  
  return [sort {$a->{start} <=> $b->{start}} values %unique_vfs];
}

# map raw diff strings to variation features
# NB not all diffs will have a single VF if there are compound effects
# resolved to a single raw diff by the sequence alignment
sub _get_raw_diff_to_VariationFeature_hash {
  my $self = shift;
  
  if(!exists($self->{_raw_to_vf})) {
    my $vfs = $self->_variation_features;
  
    my $tr_strand = $self->transcript->strand;
  
    my %raw_to_vf = ();
  
    foreach my $vf(@$vfs) {
      if(my $mapping = $vf->{_cds_mapping}) {
        my $vf_strand = $vf->strand;
        
        my @alleles = split('/', $vf->allele_string);
        my $ref = shift @alleles;
      
        reverse_comp(\$ref) if $tr_strand ne $vf_strand;
      
        foreach my $alt(@alleles) {
        
          # insertion
          if($ref eq '-') {
            reverse_comp(\$alt) if $tr_strand ne $vf_strand;
          
            $raw_to_vf{$mapping->start.'ins'.$alt} = $vf;
          }
        
          # deletion
          elsif($alt eq '-') {
            $raw_to_vf{$mapping->start.'del'.$ref} = $vf;
          }
        
          # sub
          else {
            reverse_comp(\$alt) if $tr_strand ne $vf_strand;            
          
            $raw_to_vf{$mapping->start.$ref.'>'.$alt} = $vf;
          }
        }
      }
    }
  
    $self->{_raw_to_vf} = \%raw_to_vf;
  }
  
  return $self->{_raw_to_vf};
}

## Prefetches everything that would otherwise be lazy-loaded
## This is mainly used to populate hash keys that need to be present when
## creating JSON output
sub _prefetch_everything {
  my $self = shift;
  
  foreach my $haplo(@{$self->get_all_TranscriptHaplotypes}) {
    $haplo->$_ for qw(name count frequency get_all_diffs get_all_population_frequencies);
  }
}

## use this to tell the JSON serialiser to delete certain keys
sub _dont_export {
  my $self = shift;
  my $key = shift;

  $self->{_dont_export} ||= {};
  $self->{_dont_export}->{$key} = 1 if $key;

  return $self->{_dont_export};
}

## Convert this object to a hash that can be written as JSON.
## Basically just prefetches everything that would otherwise be lazy loaded,
## deletes "private" keys starting with "_", and presents the haplotypes in
## order of frequency
sub TO_JSON {
  my $self = shift;

  # prefetch counts, predictions etc
  $self->_prefetch_everything();
  
  # make a hash copy of self
  my %copy = %{$self};

  delete $copy{$_} for keys %{$self->_dont_export};
  
  # delete keys starting with "_"
  delete $copy{$_} for grep {$_ =~ /^\_/} keys %copy;
  
  # convert haplotype hashrefs to listrefs
  $copy{'cds_haplotypes'} = [map {$_->TO_JSON} sort {$b->count <=> $a->count} @{$self->get_all_CDSHaplotypes}];
  $copy{'protein_haplotypes'} = [map {$_->TO_JSON} sort {$b->count <=> $a->count} @{$self->get_all_ProteinHaplotypes}];
  
  return \%copy;
}

1;

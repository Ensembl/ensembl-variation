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

=head1 NAME

Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer

=head1 SYNOPSIS

=head1 DESCRIPTION

A container for TranscriptHaplotype objects

=cut

package Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::Variation::CDSHaplotype;
use Bio::EnsEMBL::Variation::ProteinHaplotype;

use Digest::MD5 qw(md5_hex);
use Scalar::Util qw(weaken);

=head2 new
=cut 

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $transcript = shift;
  my $gts = shift;
  my $db = shift;
  
  my $self = {
    _transcript => $transcript,
    _individual_genotype_features => $gts,
    _db => $db,
    transcript_id => $transcript->stable_id,
  };
  bless($self, $class);
  
  $self->_init();
  
  return $self;
}

sub transcript {
  my $self = shift;
  $self->{_transcript} = shift if @_;
  return $self->{_transcript};
}

sub get_all_IndividualGenotypeFeatures {
  my $self = shift;
  return $self->{_individual_genotype_features};
}

sub _variation_features {
  my $self = shift;
  my $vfs = shift;
  
  if($vfs) {
    $self->{_variation_features} = $vfs;
  }
  else {
    $self->{_variation_features} = [map {$_->variation_feature} @{$self->get_all_IndividualGenotypeFeatures}];
  }
  return $self->{_variation_features};
}

sub db {
  my $self = shift;
  $self->{_db} = shift if @_;
  return $self->{_db};
}

sub get_TranscriptHaplotype_by_hex {
  return $_[0]->{cds_haplotypes}->{$_[1]} || $_[0]->{protein_haplotypes}->{$_[1]} || undef;
}

sub get_TranscriptHaplotype_by_name {
  my $self = shift;
  my $name = shift;
  my @matched = grep {$_->name eq $name} @{$self->get_all_TranscriptHaplotypes};
  return scalar @matched ? $matched[0] : undef;
}

sub get_all_CDSHaplotypes {
  return [values %{$_[0]->{cds_haplotypes}}];
}

sub get_all_ProteinHaplotypes {
  return [values %{$_[0]->{protein_haplotypes}}];
}

sub get_all_TranscriptHaplotypes {
  return [@{$_[0]->get_all_CDSHaplotypes}, @{$_[0]->get_all_ProteinHaplotypes}];
}

sub get_all_CDSHaplotypes_by_Individual {
  my $self = shift;
  my $ind = shift;
  return [grep {$_->{individuals}->{$ind->name}} @{$self->get_all_CDSHaplotypes}];
}

sub get_all_ProteinHaplotypes_by_Individual {
  my $self = shift;
  my $ind = shift;
  return [grep {$_->{individuals}->{$ind->name}} @{$self->get_all_ProteinHaplotypes}];
}

sub get_all_TranscriptHaplotypes_by_Individual {
  my $self = shift;
  my $ind = shift;
  return [grep {$_->{individuals}->{$ind->name}} @{$self->get_all_TranscriptHaplotypes}];
}

sub get_all_most_frequent_CDSHaplotypes {
  my $self = shift;
  
  if(!exists($self->{_most_frequent_cds})) {
    $self->{_most_frequent_cds} = $self->_get_most_frequent($self->get_all_CDSHaplotypes);
  }
  return $self->{_most_frequent_cds};
}

sub get_all_most_frequent_ProteinHaplotypes {
  my $self = shift;
  
  if(!exists($self->{_most_frequent_protein})) {
    $self->{_most_frequent_protein} = $self->_get_most_frequent($self->get_all_ProteinHaplotypes);
  }
  return $self->{_most_frequent_protein};
}

sub _get_most_frequent {
  my $self = shift;
  my $haplos = shift;
  
  # get hexes of all the haplotypes
  my @hexes = map {$_->hex} @$haplos;
  
  # find the maximum count from this list of hexes
  my $max_count = (sort {$a <=> $b} map {$self->{_counts}->{$_} || 0} @hexes)[-1];
  
  return [grep {$self->{_counts}->{$_->hex} == $max_count} @$haplos];
}

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

sub total_population_counts {
  my $self = shift;
  
  if(!exists($self->{total_population_counts})) {
    my $counts = {};
    
    my $hash = $self->_get_individual_population_hash();
    
    foreach my $ind(keys %$hash) {
      $counts->{$_} += (scalar keys %{$self->{protein_haplotypes}} ? 2 : 1) for keys %{$hash->{$ind}};
    }
    
    $self->{total_population_counts} = $counts;
  }
  
  return $self->{total_population_counts};
}

## Generates a two-level hash giving the populations for each individual
sub _get_individual_population_hash {
  my $self = shift;
  
  if(!exists($self->{_individual_population_hash})) {
    my $hash = {};
    
    foreach my $ind(values %{{map {$_->individual->name => $_->individual()} @{$self->get_all_IndividualGenotypeFeatures}}}) {
      $hash->{$ind->name}->{$_->name} = 1 for @{$ind->get_all_Populations};
    }
    
    $self->{_individual_population_hash} = $hash;
  }
  
  return $self->{_individual_population_hash};
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

## Creates all the TranscriptHaplotype objects for this container
sub _init {
  my $self = shift;
  
  my $tr = $self->transcript;
  my $vfs = $self->_variation_features;
  
  # cache reference sequences on transcript object
  $tr->{cds} = $tr->{_variation_effect_feature_cache}->{translateable_seq} || $tr->translateable_seq;
  $tr->{protein} = ($tr->{_variation_effect_feature_cache}->{peptide} || $tr->translation->seq).'*';
  
  # get mappings of all variation feature coordinates
  my $mappings = $self->_get_mappings;
  
  # remove any variation features that didn't get a mapping
  my @new_vars = grep {defined($mappings->{$_->{start}.'-'.$_->{end}})} @$vfs;
  $self->_variation_features(\@new_vars);
  
  # group vfs by individual
  my %by_ind;
  push @{$by_ind{$_->individual->name}}, $_ for @{$self->get_all_IndividualGenotypeFeatures};
  
  # do the transcript magic
  foreach my $ind(keys %by_ind) {
    my $obj = $self->_mutate_sequences($by_ind{$ind});
    
    # store unique alt seqs so we only align each once
    foreach my $type(qw(cds protein)) {
      
      # iterate over haplotypes (should be 2)
      foreach my $i(0..$#{$obj->{$type}}) {
        
        my $seq = $obj->{$type}->[$i];
        my $hex = md5_hex($seq);
        
        # try to fetch 
        my $haplotype = $self->get_TranscriptHaplotype_by_hex($hex);
        
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
        $haplotype->{individuals}->{$ind}++;
        $self->{_counts}->{$hex}++;
        
        # store mapping between hashes of cds and protein
        my $other_type = $type eq 'cds' ? 'protein' : 'cds';
        my $other_hex  = md5_hex($obj->{$other_type}->[$i]);
        
        $haplotype->{other_hexes}->{$other_hex} = 1;
      }
    }
  }
  
  # clear cached translated sequences
  delete($self->{_translations});
}

## Uses a transcript's mapper to get genomic->CDS coord mappings for each
## VariationFeaturem, with the results keyed on $vf->{start} rather than
## $vf->seq_region_start as this saves some time otherwise spent diving into
## core code
sub _get_mappings {
  my $self = shift;
  
  if(!exists($self->{_mappings})) {
    my $tr = $self->transcript();
    
    # get unique coords
    my %coords = map {$_->{start}.'-'.$_->{end} => $_} @{$self->_variation_features};
    
    # get transcript mapper
    my $mapper = $tr->{_variation_effect_feature_cache}->{mapper} || $tr->get_TranscriptMapper;
    
    my %mappings;
    
    foreach my $coord(keys %coords) {
      my $vf = $coords{$coord};
      my ($s, $e) = $vf->{slice} ? ($vf->seq_region_start, $vf->seq_region_end) : split '-', $coord;
      
      my @mapped = $mapper->genomic2cds($s, $e, $tr->strand);
      
      # clean mapping
      if(scalar @mapped == 1 && !$mapped[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
        $mappings{$coord} = $mapped[0];
      }
      
      # complex mapping, we can't deal with this ATM
      elsif(scalar @mapped > 1) {
        warn "WARNING: genomic coord $coord possibly maps across coding/non-coding boundary in ".$tr->stable_id."\n";
      }
    }
    
    $self->{_mappings} = \%mappings;
  }
  
  return $self->{_mappings};
}

## Applies a set of IndividualGenotypeFeatures to the transcript sequence
## Returns a hashref that contains the info necessary to construct a
## TranscriptHaplotype object
sub _mutate_sequences {
  my $self = shift;
  my $gts = shift;
  
  my $tr = $self->transcript;
  my $mappings = $self->_get_mappings;
  
  my $codon_table = $tr->{_variation_effect_feature_cache}->{codon_table} || $tr->slice->get_all_Attributes('codon_table')->[0]->value;
  
  my $return = {};
  
  for my $hap(0,1) {
    my $seq = $tr->{_variation_effect_feature_cache}->{translateable_seq};
    my $ref_seq = $seq;
    
    # flag if indel observed, we need to align seqs later
    my $indel = 0;
    
    # iterate through in reverse order
    foreach my $gt(reverse @$gts) {
      my $vf = $gt->variation_feature;
      my ($s, $e) = ($vf->{start}, $vf->{end});
      
      my $mapping = $mappings->{$s.'-'.$e};
      next unless $mapping;
      
      my $genotype = $gt->genotype->[$hap];
      
      # remove del characters
      $genotype =~ s/\-//g;
      
      # reverse complement sequence?
      reverse_comp(\$genotype) if $tr->strand ne $vf->strand;
      
      my $replace_len = ($mapping->end - $mapping->start) + 1;
      $indel = 1 if length($genotype) != $replace_len;
      
      substr($seq, $mapping->start - 1, $replace_len, $genotype);
    }
    
    # now translate
    my $protein = $self->_get_translation($seq, $codon_table);
    
    # remove anything beyond a stop
    #$protein =~ s/\*.+/\*/;
    
    push @{$return->{cds}}, $seq;
    push @{$return->{protein}}, $protein;
    push @{$return->{indel}}, $indel;
  }
  
  return $return;
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

## Get e.g. SIFT and PolyPhen predictions for all mutations observed
## They are then cached and can be used to annotate the diffs retrieved from
## each individual haplotype
sub _get_missense_predictions {
  my $self = shift;
  
  if(!exists($self->{_missense_predictions})) {
    
    # get SIFT and PolyPhen scores by mutation
    my $missense_preds = {};
    
    # declare $vfs here otherwise it goes out of scope and the VariationFeature
    # objects can disappear, weird side effect of weaken()?
    my ($vfs, @tvs);
    my $tr = $self->transcript();
    
    # fetch from database?
    if($vfs->[0]->{dbID}) {
      
      $vfs = $self->_variation_features;
      
      # get a TranscriptVariationAdaptor via the VariationFeature's adaptor
      my $tva = $vfs->[0]->adaptor->db->get_TranscriptVariationAdaptor();
      
      # its faster to fetch them all in one go as this only hits the DB once
      @tvs = @{$tva->fetch_all_by_VariationFeatures($vfs, [$tr])};
    }
    
    # otherwise fetch from cache on object
    else {
      
      # get VF consequences
      $vfs = $self->_get_unique_VariationFeatures_with_consequences;
      
      @tvs = map {@{$_->get_all_TranscriptVariations([$tr])}} @$vfs;
    }
    
    foreach my $tva(map {@{$_->get_all_alternate_TranscriptVariationAlleles}} @tvs) {
      next unless grep {$_->SO_term eq 'missense_variant'} @{$tva->get_all_OverlapConsequences};
      #next unless $tva->pep_allele_string;
      my $mut = $tva->transcript_variation->translation_start.$tva->pep_allele_string;
      $mut =~ s/\//\>/;
      
      $missense_preds->{$mut} ||= {
        sift_prediction => $tva->sift_prediction,
        sift_score => $tva->sift_score,
        polyphen_prediction => $tva->polyphen_prediction,
        polyphen_score => $tva->polyphen_score,
        #consequences => [map {$_->SO_term} @{$tva->get_all_OverlapConsequences}],
      };
      
      foreach my $key(keys %{$missense_preds->{$mut}}) {
        delete $missense_preds->{$mut}->{$key} unless defined $missense_preds->{$mut}->{$key};
        $missense_preds->{$mut}->{$key} =~ s/ /\_/ if ref($missense_preds->{$mut}->{$key}) eq 'SCALAR';
      }
    }
    
    $self->{_missense_predictions} = $missense_preds;
  }
  
  return $self->{_missense_predictions};
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
      delete $vf_copy->{$_} for qw(_line phased individual genotype non_variant hom_ref);
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
  
  # delete keys starting with "_"
  delete $copy{$_} for grep {$_ =~ /^\_/} keys %copy;
  
  # convert haplotype hashrefs to listrefs
  $copy{'cds_haplotypes'} = [sort {$b->count <=> $a->count} @{$self->get_all_CDSHaplotypes}];
  $copy{'protein_haplotypes'} = [sort {$b->count <=> $a->count} @{$self->get_all_ProteinHaplotypes}];
  
  return \%copy;
}

1;

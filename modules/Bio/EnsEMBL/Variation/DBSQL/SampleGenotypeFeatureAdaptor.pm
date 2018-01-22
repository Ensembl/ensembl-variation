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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::SampleGenotypeFeatureAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::SampleGenotypeFeatureAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $sga = $reg->get_adaptor('human', 'variation', 'samplegenotype');

  #returns all genotypes in a certain Slice

  $genotypes = $sga->fetch_by_Slice($slice);



=head1 DESCRIPTION

This adaptor provides database connectivity for SampleGenotypeFeature objects.
SampleGenotypeFeatures may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

package Bio::EnsEMBL::Variation::DBSQL::SampleGenotypeFeatureAdaptor;

use strict;
use warnings;

use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::SampleGenotypeFeature;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor);


=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation $variation
  Example    : my $var = $variation_adaptor->fetch_by_name( "rs1121" )
               $sgtypes = $sgtype_adaptor->fetch_all_by_Variation( $var )
  Description: Retrieves a list of sample genotypes for the given Variation.
               If none are available an empty listref is returned.
  Returntype : listref Bio::EnsEMBL::Variation::SampleGenotype 
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub fetch_all_by_Variation {
  my $self = shift;
  my $variation = shift;
  my $sample = shift;

  my $res;
  if (!ref($variation) || !$variation->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation argument expected');
  }

  if (!defined($variation->dbID())) {
    warning("Cannot retrieve genotypes for variation without set dbID");
    return [];
  }

  $self->{_variation_id} = $variation->dbID;

  # foreach of the hitting variation Features, get the Genotype information
  foreach my $vf (@{$variation->get_all_VariationFeatures}){

    # get the feature slice for this VF
    my $fs = $vf->feature_Slice();

    # if the feature slice is start > end
    if($fs->start > $fs->end) {

      # get a new slice with the start and end the right way round
      # otherwise the call won't pick any variations up
      my $new_fs = $fs->{'adaptor'}->fetch_by_region($fs->coord_system->name,$fs->seq_region_name,$fs->end,$fs->start);
      $fs = $new_fs;
    }

    map {$_->variation($variation); push @{$res}, $_} @{$self->fetch_all_by_Slice($fs, $sample)};
  }

  delete $self->{_variation_id};

  return $res;
}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL:Slice $slice
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Sample $sample
  Example    : my @SampleGenotypesFeatures = @{$ca->fetch_all_by_Slice($slice)};
  Description: Retrieves all SampleGenotypeFeature features for a given slice for
               a certain sample (if provided). 
  Returntype : reference to list Bio::EnsEMBL::Variation::SampleGenotypeFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $sample = shift;

  my @results;
  my $features;
  my $constraint;

  #if passed inividual or population, add constraint
  if (defined $sample && defined $sample->dbID){
    my $sample_str;

    if ($sample->isa("Bio::EnsEMBL::Variation::Population")) {
      my $samples = $sample->get_all_Samples;
      my @list;
      push @list, $_->dbID foreach @$samples;
      $sample_str = (@list > 1)  ? " IN (".join(',',@list).")"   :   ' = \''.$list[0].'\'';
      $constraint = " c.sample_id $sample_str";
    }
    else {
      $constraint = ' c.sample_id = ' . $sample->dbID;
    }
  }

  my $use_vcf = $self->db->use_vcf();
  if ($use_vcf) { 
    @$features = map {@{$_->get_all_SampleGenotypeFeatures_by_Slice($slice, $sample)}} @{$self->db->get_VCFCollectionAdaptor->fetch_all() || []};
  }
  if ($use_vcf <= 1) {
    push @$features, @{$self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint)};
  }

  my $seq_region_slice = $slice->seq_region_Slice;
  foreach my $sampleFeature (@{$features}){
    if ($sampleFeature->start > 0 && ($slice->end - $slice->start + 1) >= $sampleFeature->end){	
      if ($sampleFeature->slice->strand == -1){ #ignore the different strand transformation
        # Position will change if the strand is negative so change the strand to 1 temporarily
        $sampleFeature->slice->{'strand'} = 1;
        my $newFeature = $sampleFeature->transfer($seq_region_slice); 
        $sampleFeature->slice->{'strand'} = -1;
        $newFeature->slice->{'strand'} = -1;
        $newFeature->variation($sampleFeature->variation);
        push @results, $newFeature;
      }
      else {
        push @results, $sampleFeature->transfer($seq_region_slice);
      }
    } 
  } 
  return \@results;
}

=head2 fetch_all_unique_by_Slice
  Arg [1]    : Bio::EnsEMBL:Slice $slice
  Arg [2]    : Bio::EnsEMBL::Variation::Sample $sample
  Arg [3]    : Bio::EnsEMBL::Variation::Population $population
  Example    : my @sample_genotype_features = @{$sgf_adaptor->fetch_all_unique_by_Slice($slice, $sample, $population)};
  Description: Retrieves all SampleGenotypeFeatures for a given slice that are unique
               for the given sample in the given population.
               The allowed samples and populations are restricted to data that comes from VCF files. 
               Use the Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor::fetch_all_vcf_Populations to retrieve all populations
               that are represented in VCF files.
  Returntype : reference to list Bio::EnsEMBL::Variation::SampleGenotypeFeature
  Exceptions : throw on bad and missing argument
  Caller     : general
  Status     : Stable
=cut

sub fetch_all_unique_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $sample = shift;
  my $population = shift;

  if (!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw("Bio::EnsEMBL::Slice arg expected");
  }
  if (!ref($sample) || !$sample->isa('Bio::EnsEMBL::Variation::Sample')) {
    throw("Bio::EnsEMBL::Variation::Sample arg expected");
  }
  if (!ref($population) || !$population->isa('Bio::EnsEMBL::Variation::Population')) {
    throw("Bio::EnsEMBL::Variation::Population arg expected");
  }

  my $use_vcf = $self->db->use_vcf();
  if (!$use_vcf) {
    warning('You need to set use_vcf: $sample_genotype_feature_adaptor->db->use_vcf(1)');
    return [];
  }
  my $non_ref_only = 0;
  my @features = map {@{$_->get_all_SampleGenotypeFeatures_by_Slice($slice, $population, $non_ref_only, $sample, 'unique')}} @{$self->db->get_VCFCollectionAdaptor->fetch_all() || []};
  return \@features;
}

=head2 fetch_all_differences_by_Slice
  Arg [1]    : Bio::EnsEMBL:Slice $slice
  Arg [2]    : Bio::EnsEMBL::Variation::Sample $sample
  Arg [3]    : Bio::EnsEMBL::Variation::Population $population
  Example    : my @sample_genotype_features = @{$sgf_adaptor->fetch_all_differences_by_Slice($slice, $sample, $population)};
  Description: Retrieves all SampleGenotypeFeatures for a given slice and sample. It also stores genotypes for samples that differ
               from the given sample's genotype and are also different from the reference.
               The allowed samples and populations are restricted to data that comes from VCF files. 
               Use the Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor::fetch_all_vcf_Populations to retrieve all populations
               that are represented in VCF files.
  Returntype : reference to list Bio::EnsEMBL::Variation::SampleGenotypeFeature
  Exceptions : throw on bad and missing argument
  Caller     : general
  Status     : Stable
=cut

sub fetch_all_differences_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $sample = shift;
  my $population = shift;

  if (!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw("Bio::EnsEMBL::Slice arg expected");
  }
  if (!ref($sample) || !$sample->isa('Bio::EnsEMBL::Variation::Sample')) {
    throw("Bio::EnsEMBL::Variation::Sample arg expected");
  }
  if (!ref($population) || !$population->isa('Bio::EnsEMBL::Variation::Population')) {
    throw("Bio::EnsEMBL::Variation::Population arg expected");
  }

  my $use_vcf = $self->db->use_vcf();
  if (!$use_vcf) {
    warning('You need to set use_vcf: $sample_genotype_feature_adaptor->db->use_vcf(1)');
    return [];
  }
  my $non_ref_only = 0;
  my @features = map {@{$_->get_all_SampleGenotypeFeatures_by_Slice($slice, $population, $non_ref_only, $sample, 'differences')}} @{$self->db->get_VCFCollectionAdaptor->fetch_all() || []};
  return \@features;
}

sub _tables {
    return (['compressed_genotype_region','c']);
}

sub _columns {
    return qw(c.sample_id c.seq_region_id c.seq_region_start c.seq_region_end c.seq_region_strand c.genotypes);
}

sub _write_columns {
    return $_[0]->_columns;
}

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  my $sa = $self->db()->dnadb()->get_SliceAdaptor();

  my (@results, %slice_hash, %sr_name_hash, %sr_cs_hash, %sample_hash, %gt_code_hash);

  my ($sample_id, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $genotypes);

  $sth->bind_columns(\$sample_id, \$seq_region_id, \$seq_region_start, \$seq_region_end, \$seq_region_strand, \$genotypes);

	my $asm_cs;
	my $cmp_cs;
	my $asm_cs_vers;
	my $asm_cs_name;
	my $cmp_cs_vers;
	my $cmp_cs_name;
	if($mapper) {
		$asm_cs = $mapper->assembled_CoordSystem();
		$cmp_cs = $mapper->component_CoordSystem();
		$asm_cs_name = $asm_cs->name();
		$asm_cs_vers = $asm_cs->version();
		$cmp_cs_name = $cmp_cs->name();
		$cmp_cs_vers = $cmp_cs->version();
	}

	my $dest_slice_start;
	my $dest_slice_end;
	my $dest_slice_strand;
	my $dest_slice_length;
	if($dest_slice) {
		$dest_slice_start  = $dest_slice->start();
		$dest_slice_end    = $dest_slice->end();
		$dest_slice_strand = $dest_slice->strand();
		$dest_slice_length = $dest_slice->length();
	}

	FEATURE: while($sth->fetch()) {
		#get the slice object
		my $slice = $slice_hash{"ID:".$seq_region_id};
		
		if(!$slice) {
			$slice = $sa->fetch_by_seq_region_id($seq_region_id);
			$slice_hash{"ID:".$seq_region_id} = $slice;
			$sr_name_hash{$seq_region_id} = $slice->seq_region_name();
			$sr_cs_hash{$seq_region_id} = $slice->coord_system();
		}
		
		my $orig_start = $seq_region_start;
		
		#
		# remap the feature coordinates to another coord system
		# if a mapper was provided
		#
		if($mapper) {
			my $sr_name = $sr_name_hash{$seq_region_id};
			my $sr_cs   = $sr_cs_hash{$seq_region_id};
			
			($sr_name,$seq_region_start,$seq_region_end,$seq_region_strand) =
			$mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
			$seq_region_strand, $sr_cs);
			
			#skip features that map to gaps or coord system boundaries
			next FEATURE if(!defined($sr_name));
			
			#get a slice in the coord system we just mapped to
			if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
				$slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"} ||=
				$sa->fetch_by_region($cmp_cs_name, $sr_name,undef, undef, undef,
				$cmp_cs_vers);
			} else {
				$slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"} ||=
				$sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef,
				$asm_cs_vers);
			}
		}
		
		#
		# If a destination slice was provided convert the coords
		# If the dest_slice starts at 1 and is foward strand, nothing needs doing
		#
		if($dest_slice) {
			if($dest_slice_start != 1 || $dest_slice_strand != 1) {
				if($dest_slice_strand == 1) {
					$seq_region_start = $seq_region_start - $dest_slice_start + 1;
					$seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
				} else {
					my $tmp_seq_region_start = $seq_region_start;
					#		$seq_region_start = $dest_slice_end - $seq_region_end + 1;
					$seq_region_start = $seq_region_start - $dest_slice_start + 1;
					$seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
					#		$seq_region_strand *= -1;
				}
				
				#throw away features off the end of the requested slice
				if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
					next FEATURE;
				}
			}
			$slice = $dest_slice;
		}	
		
		#my @genotypes = unpack '(ww)*', $genotypes;
		my @genotypes = unpack '(www)*', $genotypes;
		my $snp_start = $seq_region_start;
		
		#while( my( $gt_code, $gap ) = splice @genotypes, 0, 2 ) {
		while( my( $variation_id, $gt_code, $gap ) = splice @genotypes, 0, 3 ) {
			
      # if fetching by variation ID, skip those without
			if(defined($self->{_variation_id})) {
				if($variation_id != $self->{_variation_id}) {
					$snp_start += $gap + 1 if defined $gap;
					next;
				}
			}
      
      # if fetching by slice, skip those outside slice range
      if($dest_slice && ($snp_start < 1 || $snp_start > $dest_slice_length)) {
        $snp_start += $gap + 1 if defined $gap;
        next;
      }
			
			my $sgtype  = Bio::EnsEMBL::Variation::SampleGenotypeFeature->new_fast({
				'start'   => $snp_start,
				'end'     => $snp_start,
				'strand'  => $seq_region_strand,
				'slice'   => $slice,
				'gt_code' => $gt_code,
				'adaptor' => $self,
				'_variation_id' => $variation_id
			});
			
			$sample_hash{$sample_id} ||= [];
			push @{$sample_hash{$sample_id}}, $sgtype;
			
			$gt_code_hash{$gt_code} ||= [];
			push @{$gt_code_hash{$gt_code}}, $sgtype;
			
			push @results, $sgtype;
			$snp_start += $gap + 1 if defined $gap;
		}
	}
	
	# get all sample in one query (faster)
	# and add to already created genotypes	
	my $sample_adpt = $self->db()->get_SampleAdaptor();
	my $samples = $sample_adpt->fetch_all_by_dbID_list([keys %sample_hash]);
	
	foreach my $s (@$samples) {
		foreach my $sgty (@{$sample_hash{$s->dbID()}}) {
			$sgty->{sample} = $s;
		}
	}
	
	# get all genotypes from codes
	my $gtca = $self->db->get_GenotypeCodeAdaptor();
	my $gtcs = $gtca->fetch_all_by_dbID_list([keys %gt_code_hash]);
	
	foreach my $gtc(@$gtcs) {
    foreach my $sgty(@{$gt_code_hash{$gtc->dbID}}) {
      $sgty->{genotype} = $gtc->genotype;
      $sgty->{phased}   = $gtc->phased;
    }
	}
	
	# unique sort the results on sample and position (we don't care if GTs disagree)
	my %tmp_hash = map {$_->{sample}."_".$_->{start}."_".$_->{end} => $_} @results;
	
	return [values %tmp_hash];
}

1;

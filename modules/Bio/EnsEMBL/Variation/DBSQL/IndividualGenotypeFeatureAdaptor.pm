=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeFeatureAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeFeatureAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $iga = $reg->get_adaptor("human","variation","individualgenotype");

  #returns all genotypes in a certain Slice

  $genotypes = $iga->fetch_by_Slice($slice);



=head1 DESCRIPTION

This adaptor provides database connectivity for IndividualGenotypeFeature objects.
IndividualGenotypeFeatures may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

package Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeFeatureAdaptor;

use strict;
use warnings;

use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::IndividualGenotypeFeature;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor);


=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation $variation
  Example    : my $var = $variation_adaptor->fetch_by_name( "rs1121" )
               $igtypes = $igtype_adaptor->fetch_all_by_Variation( $var )
  Description: Retrieves a list of individual genotypes for the given Variation.
               If none are available an empty listref is returned.
  Returntype : listref Bio::EnsEMBL::Variation::IndividualGenotype 
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub fetch_all_by_Variation {
    my $self = shift;
    my $variation = shift;
	my $individual = shift;

    my $res;
    if(!ref($variation) || !$variation->isa('Bio::EnsEMBL::Variation::Variation')) {
		throw('Bio::EnsEMBL::Variation::Variation argument expected');
    }

    if(!defined($variation->dbID())) {
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
		
		map {$_->variation($variation); push @{$res}, $_} @{$self->fetch_all_by_Slice($fs, $individual)};
	}
	
	delete $self->{_variation_id};
	
    return $res;
}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL:Slice $slice
  Arg [2]    : (optional) Bio::EnsEMBL::Variation::Individual $individual
  Example    : my @IndividualGenotypesFeatures = @{$ca->fetch_all_by_Slice($slice)};
  Description: Retrieves all IndividualGenotypeFeature features for a given slice for
               a certain individual (if provided). 
  Returntype : reference to list Bio::EnsEMBL::Variation::IndividualGenotype
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice{
    my $self = shift;
    my $slice = shift;
    my $individual = shift;
	
    my @results;
    my $features;
    my $constraint;
	
	#if passed inividual or population, add constraint
	if (defined $individual && defined $individual->dbID){
		my $instr;
		
		if($individual->isa("Bio::EnsEMBL::Variation::Population")) {
		  my $inds = $individual->get_all_Individuals;
		  my @list;
		  push @list, $_->dbID foreach @$inds;
		  $instr = (@list > 1)  ? " IN (".join(',',@list).")"   :   ' = \''.$list[0].'\'';
		  $constraint = " c.individual_id $instr";
		}
		else {
		  $constraint = ' c.individual_id = ' . $individual->dbID;
		}
	}
	
	$features = $self->SUPER::fetch_all_by_Slice_constraint($slice,$constraint);
	
	
	my $seq_region_slice = $slice->seq_region_Slice;

	foreach my $indFeature (@{$features}){
		if ($indFeature->start > 0 && ($slice->end-$slice->start +1) >= $indFeature->end){	
			if ($indFeature->slice->strand == -1){ #ignore the different strand transformation
			
				# Position will change if the strand is negative so change the strand to 1 temporarily
				$indFeature->slice->{'strand'} = 1;
				my $newFeature = $indFeature->transfer($seq_region_slice); 
				$indFeature->slice->{'strand'} = -1;
				$newFeature->slice->{'strand'} = -1;
				$newFeature->variation($indFeature->variation);
				push @results, $newFeature;
			}
			else{
				push @results,$indFeature->transfer($seq_region_slice);
			}
		}
	}
	
    return \@results;
    
}

sub _tables{
    return (['compressed_genotype_region','c']);
}

sub _columns{
    return qw(individual_id seq_region_id seq_region_start seq_region_end seq_region_strand genotypes);
}

sub _write_columns{
    return $_[0]->_columns;
}

sub _objs_from_sth{
    my ($self, $sth, $mapper, $dest_slice) = @_;
	
    #
    # This code is ugly because an attempt has been made to remove as many
    # function calls as possible for speed purposes.  Thus many caches and
    # a fair bit of gymnastics is used.
    #
    
    my $sa = $self->db()->dnadb()->get_SliceAdaptor();
    
    my (@results, %slice_hash, %sr_name_hash, %sr_cs_hash, %individual_hash, %gt_code_hash);

	my ($individual_id, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $genotypes);
	
	$sth->bind_columns(
		\$individual_id, \$seq_region_id, \$seq_region_start,
		\$seq_region_end, \$seq_region_strand, \$genotypes
	);

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
			
			if(defined($self->{_variation_id})) {
				if($variation_id != $self->{_variation_id}) {
					$snp_start += $gap + 1 if defined $gap;
					next;
				}
			}
			
			my $igtype  = Bio::EnsEMBL::Variation::IndividualGenotypeFeature->new_fast({
				'start'   => $snp_start,
				'end'     => $snp_start,
				'strand'  => $seq_region_strand,
				'slice'   => $slice,
				'gt_code' => $gt_code,
				'adaptor' => $self,
				'_variation_id' => $variation_id
			});
			
			$individual_hash{$individual_id} ||= [];
			push @{$individual_hash{$individual_id}}, $igtype;
			
			$gt_code_hash{$gt_code} ||= [];
			push @{$gt_code_hash{$gt_code}}, $igtype;
			
			push @results, $igtype;
			$snp_start += $gap + 1 if defined $gap;
		}
	}
	
	# get all individual in one query (faster)
	# and add to already created genotypes	
	my $ia = $self->db()->get_IndividualAdaptor();
	my $inds = $ia->fetch_all_by_dbID_list([keys %individual_hash]);
	
	foreach my $i (@$inds) {
		foreach my $igty (@{$individual_hash{$i->dbID()}}) {
			$igty->{individual} = $i;
		}
	}
	
	# get all genotypes from codes
	my $gtca = $self->db->get_GenotypeCodeAdaptor();
	my $gtcs = $gtca->fetch_all_by_dbID_list([keys %gt_code_hash]);
	
	foreach my $gtc(@$gtcs) {
		foreach my $igty(@{$gt_code_hash{$gtc->dbID}}) {
			$igty->{genotype} = $gtc->genotype;
		}
	}
	
	# unique sort the results on individual and position (we don't care if GTs disagree)
	my %tmp_hash = map {$_->{individual}."_".$_->{start}."_".$_->{end} => $_} @results;
	
	return [values %tmp_hash];
}

1;

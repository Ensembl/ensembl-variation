#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor
#
# Copyright (c) 2005 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor

=head1 SYNOPSIS

  $db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  $iga = $db->get_IndividualGenotypeAdaptor();

  #returns all genotypes in a certain Slice

  $genotypes = $iga->fetch_by_Slice($slice);



=head1 DESCRIPTION

This adaptor provides database connectivity for IndividualGenotype objects.
IndividualGenotypes may be retrieved from the Ensembl variation database by
several means using this module.

=head1 AUTHOR - Daniel Rios

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor;

use strict;
use warnings;

use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::IndividualGenotype;
use Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor;

@ISA = qw(Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor);

use Data::Dumper;

=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation $variation
  Example    : my $var = $variation_adaptor->fetch_by_name( "rs1121" )
               $igtypes = $igtype_adaptor->fetch_all_by_Variation( $var )
  Description: Retrieves a list of individual genotypes for the given Variation.
               If none are available an empty listref is returned.
  Returntype : listref Bio::EnsEMBL::Variation::IndividualGenotype 
  Exceptions : none
  Caller     : general

=cut


sub fetch_all_by_Variation {
    my $self = shift;
    my $variation = shift;

    my $res;
    if(!ref($variation) || !$variation->isa('Bio::EnsEMBL::Variation::Variation')) {
	throw('Bio::EnsEMBL::Variation::Variation argument expected');
    }

    if(!defined($variation->dbID())) {
	warning("Cannot retrieve genotypes for variation without set dbID");
	return [];
    }	
    my $vfa = $self->db->get_VariationFeatureAdaptor();
    if (!$vfa){
	throw("Cannot retrieve genotypes for variation without adaptor set");
	return [];
    }
    my $variation_features = $vfa->fetch_all_by_Variation($variation);
    #foreach of the hitting variation Features, get the Genotype information
    foreach my $vf (@{$variation_features}){
	map {$_->variation($variation); push @{$res}, $_} @{$self->fetch_all_by_Slice($vf->feature_Slice)};
    }
    #and include the genotypes from the multiple genotype table
    $self->_multiple(1);
    push @{$res}, @{$self->SUPER::fetch_all_by_Variation($variation)};
    return $res;
}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL:Slice $slice
  Example    : my @IndividualGenotypesFeatures = @{$ca->fetch_all_by_Slice($slice)};
  Description: Retrieves all IndividualGenotypeFeature features for a given slice. 
  Returntype : reference to list Bio::EnsEMBL::Variation::IndividualGenotype
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub fetch_all_by_Slice{
    my $self = shift;
    my $slice = shift;
    my @results;
    if (!$self->_multiple){
	my $features = $self->SUPER::fetch_all_by_Slice($slice);
	#need to check the feature is within the Slice
	foreach my $indFeature (@{$features}){
	    if ($indFeature->start > 0 && ($slice->end-$slice->start +1) >= $indFeature->end){
		push @results,$indFeature->transfer($slice->seq_region_Slice);
	    }
	}
    }
    else{
	#and include the genotypes from the multiple genotype table
	push @results, @{$self->SUPER::fetch_all_by_Slice($slice)};
    }
    return \@results;
    
}


sub _tables{
    my $self = shift;
    return (['compressed_genotype_single_bp','c']) if (!$self->_multiple);
    return $self->SUPER::_tables if ($self->_multiple);
}

sub _columns{
    my $self = shift;

    return qw(sample_id seq_region_id seq_region_start seq_region_end seq_region_strand genotypes) if (!$self->_multiple);
    return $self->SUPER::_columns if ($self->_multiple);
}

sub _objs_from_sth{
    my ($self, $sth, $mapper, $dest_slice) = @_;
    
    return $self->SUPER::_objs_from_sth($sth,$mapper,$dest_slice) if ($self->_multiple);
    #
    # This code is ugly because an attempt has been made to remove as many
    # function calls as possible for speed purposes.  Thus many caches and
    # a fair bit of gymnastics is used.
    #
    
    my $sa = $self->db()->dnadb()->get_SliceAdaptor();
    
    my @results;
    my %slice_hash;
    my %sr_name_hash;
    my %sr_cs_hash;
    my %individual_hash;

  my ($sample_id, $seq_region_id, $seq_region_start,
      $seq_region_end, $seq_region_strand, $genotypes);

  $sth->bind_columns(\$sample_id, \$seq_region_id,
                     \$seq_region_start, \$seq_region_end, \$seq_region_strand,
                     \$genotypes);

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
		$seq_region_start = $dest_slice_end - $seq_region_end + 1;
		$seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
		$seq_region_strand *= -1;
	    }
	    
	    #throw away features off the end of the requested slice
	    if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
		next FEATURE;
	    }
	}
	$slice = $dest_slice;
    }
    #create the different Features for all the Genotypes compressed in the genotype field
    my $blob = substr($genotypes,2);
    my @genotypes = unpack("naa" x (length($blob)/4),$blob);

    unshift @genotypes, substr($genotypes,1,1); #add the second allele of the first genotype
    unshift @genotypes, substr($genotypes,0,1); #add the first allele of the first genotype
    unshift @genotypes, 0; #the first SNP is in the position indicated by the seq_region1
    my ($snp_start, $allele_1, $allele_2);

    for (my $i=0; $i < @genotypes -1;$i+=3){
	#number of gaps
	if ($i == 0){
	    $snp_start = $seq_region_start; #first SNP is in the beginning of the region
	}
	else{
	    $snp_start += $genotypes[$i] +1;
	}
	#genotype
	$allele_1 = $genotypes[$i+1];
	$allele_2 = $genotypes[$i+2];

	my $igtype = Bio::EnsEMBL::Variation::IndividualGenotype->new_fast({
	    'start'    => $snp_start,
	    'end'      => $snp_start,
	    'strand'   => $seq_region_strand,
	    'slice'    => $slice,	    
	    'allele1'  => $allele_1,
	    'allele2' => $allele_2
	});
	$individual_hash{$sample_id} ||= [];
	push @{$individual_hash{$sample_id}}, $igtype;
	push @results, $igtype;
    }

}
     # get all individual in one query (faster)
     # and add to already created genotypes
     my @ind_ids = keys %individual_hash;
     
     my $ia = $self->db()->get_IndividualAdaptor();
     my $inds = $ia->fetch_all_by_dbID_list(\@ind_ids);
     
     foreach my $i (@$inds) {
	 foreach my $igty (@{$individual_hash{$i->dbID()}}) {
	     $igty->individual($i);
	 }
     }

     return \@results;
}

sub _default_where_clause  {
    my $self = shift;
    return '' if (!$self->_multiple);
    return $self->SUPER::_default_where_clause if ($self->_multiple);
}

1;

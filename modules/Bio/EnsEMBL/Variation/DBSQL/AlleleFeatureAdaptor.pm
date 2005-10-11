#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::AlleleFeatureAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::AlleleFeatureAdaptor

=head1 SYNOPSIS

  $vdb = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);
  $db  = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  # tell the variation database where core database information can be
  # be found
  $vdb->dnadb($db);

  $afa = $vdb->get_AlleleFeatureAdaptor();
  $sa  = $db->get_SliceAdaptor();

  # Get a VariationFeature by its internal identifier
  $af = $afa->fetch_by_dbID(145);

  # get all AlleleFeatures in a region
  $slice = $sa->fetch_by_region('chromosome', 'X', 1e6, 2e6);
  foreach $af (@{$afa->fetch_all_by_Slice($slice)}) {
    print $af->start(), '-', $af->end(), ' ', $af->allele(), "\n";
  }


=head1 DESCRIPTION

This adaptor provides database connectivity for AlleleFeature objects.
Genomic locations of alleles in samples can be obtained from the 
database using this adaptor.  See the base class BaseFeatureAdaptor for more information.

=head1 AUTHOR - Daniel Rios

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::AlleleFeatureAdaptor;

use Bio::EnsEMBL::Variation::AlleleFeature;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Sequence qw(expand);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(ambiguity_code);

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');


=head2 from_IndividualSlice
    
    Arg[0]      : (optional) int $number
    Example     : $afa->from_IndividualSlice(1);
    Description : Getter/Setter to know wether the Adaptor has been called
                  from and IndividualSlice (and must create the objects from the
		  genotype table) or from the StrainSlice (and must create the objects
                  from the allele table)
    ReturnType  : int
    Exceptions  : none
    Caller      : general    

=cut

sub from_IndividualSlice{
    my $self = shift;

    return $self->{'from_IndividualSlice'} = shift if(@_);
    return $self->{'from_IndividualSlice'};    
}

=head2 fetch_all_by_Slice_Population

   Arg[0]      : Bio::EnsEMBL::Slice $slice
   Arg[1]      : Bio::EnsEMBL::Variation::Population $population
   Example     : my $vf = $vfa->fetch_all_by_Slice_Individual($slice,$population);   
   Description : Gets all the VariationFeatures in a certain Slice for a given
                 Population
   ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
   Exceptions  : thrown on bad arguments
   Caller      : general
   
=cut

sub fetch_all_by_Slice_Population{
    my $self = shift;
    my $slice = shift;
    my $population = shift;

    if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
	throw('Bio::EnsEMBL::Slice arg expected');
    }

    if(!ref($population) || !$population->isa('Bio::EnsEMBL::Variation::Population')) {
	throw('Bio::EnsEMBL::Variation::Population arg expected');
    }
    if(!defined($population->dbID())) {
	throw("Population arg must have defined dbID");
    }
    
    my $constraint = "a.sample_id = " . $population->dbID;
    #call the method fetch_all_by_Slice_constraint with the population constraint
    return $self->fetch_all_by_Slice_constraint($slice,$constraint);    
}


=head2 fetch_all_IndividualSlice

   Arg[0]      : Bio::EnsEMBL::Slice $slice
   Arg[1]      : Bio::EnsEMBL::Variation::Population $population
   Example     : my $vf = $vfa->fetch_all_by_Slice($slice,$population);   
   Description : Gets all the VariationFeatures in a certain Slice for a given
                 Population
   ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
   Exceptions  : thrown on bad arguments
   Caller      : general
   
=cut

sub fetch_all_IndividualSlice{
    my $self = shift;
    my $slice = shift;

    if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
	throw('Bio::EnsEMBL::Slice arg expected');
    }


    #indicate to select the data from the single_bp table
    $self->_multiple_bp(0);
    #call the method fetch_all_by_Slice_constraint with the population constraint for the single_bp
    my $allele_features = $self->fetch_all_by_Slice($slice);
    
    #remove from the cache the result of the query (would not run the second one)
#    my $key = uc(join(':', $slice->name, $constraint));

#    delete $self->{'_slice_feature_cache'}->{$key};
    #indicate to select the data from the multiple_bp table
    $self->_multiple_bp(1);
    #and add the same information for the multiple_bp table

#    my $allele_features_multiple = $self->fetch_all_by_Slice_constraint($slice,$constraint);
#    push @{$allele_features},$allele_features_multiple if (@{$allele_features_multiple} > 0);

    #and include then the result of both queries
#    $self->{'_slice_feature_cache'}->{$key} = $allele_features;
    return $allele_features;
}

sub _tables{    
    my $self = shift;

    if ($self->from_IndividualSlice){	
	return (['variation_feature','vf'], ['individual_genotype_single_bp','ig']) if (!$self->_multiple_bp());
	return (['variation_feature','vf'], ['individual_genotype_multiple_bp','ig']) if ($self->_multiple_bp());
    }
    else{
	return (['variation_feature','vf'],   ['allele','a']);
    }

}

sub _columns{
    my $self = shift;

    return ('vf.variation_id' ,'ig.sample_id', 'CONCAT(ig.allele_1,"|",ig.allele_2) as alleles',
	      'vf.seq_region_id', 'vf.seq_region_start', 'vf.seq_region_end', 
	      'vf.seq_region_strand', 'vf.variation_name') if ($self->from_IndividualSlice());

    return qw(a.variation_id a.sample_id a.allele 
	      vf.seq_region_id vf.seq_region_start vf.seq_region_end 
	      vf.seq_region_strand vf.variation_name);
}

sub _default_where_clause{
    my $self = shift;
    return "ig.variation_id = vf.variation_id" if ($self->from_IndividualSlice());
    return "a.variation_id = vf.variation_id";
}

sub _objs_from_sth{
 my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  my $sa = $self->db()->dnadb()->get_SliceAdaptor();

  my @features;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my ($variation_id, $sample_id, $allele,$seq_region_id,
      $seq_region_start,$seq_region_end, $seq_region_strand, $variation_name );

  $sth->bind_columns(\$variation_id,\$sample_id,\$allele,
		     \$seq_region_id,\$seq_region_start,\$seq_region_end,\$seq_region_strand,
		     \$variation_name);

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
    $allele = ambiguity_code($allele) if ($self->from_IndividualSlice);
    push @features, Bio::EnsEMBL::Variation::AlleleFeature->new_fast(
								     {'start'    => $seq_region_start,
								      'end'      => $seq_region_end,
								      'strand'   => $seq_region_strand,
								      'slice'    => $slice,
								      'allele_string' => $allele,
								      'variation_name' => $variation_name,
								      'adaptor'  => $self,
								      '_variation_id' => $variation_id,
								      '_sample_id' => $sample_id});      
}
 return\@features;
}

#internal function used to determine wether selecting data from the individual_genotype_single_bp table or the multiple_bp
sub _multiple_bp{
    my $self = shift;
    $self->{'multiple_bp'} = shift if (@_);
    return $self->{'multiple_bp'};
}

1;

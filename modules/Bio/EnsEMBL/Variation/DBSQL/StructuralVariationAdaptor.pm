
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $sa = $reg->get_adaptor("human","core","slice");
  $sva = $reg->get_adaptor("human","variation","structuralvariation");
  
  # Get a StructuralVariation by its internal identifier
  $sv = $sva->fetch_by_dbID(145);

  # get all StructuralVariations in a region
  $slice = $sa->fetch_by_region('chromosome', 'X', 1e6, 2e6);
  foreach $sv (@{$sva->fetch_all_by_Slice($slice)}) {
    print $sv->start(), '-', $sv->end(), ' ', $svf->class(), "\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity for StructuralVariation objects.
Genomic locations of structural variations can be obtained from the database using this
adaptor.  See the base class BaseFeatureAdaptor for more information.

=head1 AUTHOR - William McLaren

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor;

use Bio::EnsEMBL::Variation::StructuralVariation;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');


# method used by superclass to construct SQL
sub _tables { return (['structural_variation', 'sv'],
		      [ 'source', 's']); }


sub _default_where_clause {
  my $self = shift;

  return 'sv.source_id = s.source_id';
}

sub _columns {
  return qw( sv.structural_variation_id sv.seq_region_id sv.seq_region_start
             sv.seq_region_end sv.seq_region_strand
             sv.variation_name s.name s.description sv.class
			 sv.bound_start sv.bound_end);
}



sub _objs_from_sth {
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

  my ($struct_variation_id, $seq_region_id, $seq_region_start,
      $seq_region_end, $seq_region_strand,
      $variation_name, $source_name, $source_description, $sv_class,
	  $bound_start, $bound_end);

  $sth->bind_columns(\$struct_variation_id, \$seq_region_id,
                     \$seq_region_start, \$seq_region_end, \$seq_region_strand,
                     \$variation_name, \$source_name, \$source_description,
					 \$sv_class, \$bound_start, \$bound_end);

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
	
    push @features, $self->_create_feature_fast('Bio::EnsEMBL::Variation::StructuralVariation',
    #push @features, Bio::EnsEMBL::Variation::StructuralVariation->new_fast(
    #if use new_fast, then do not need "-" infront of key, i.e 'start' => $seq_region_start,

      {'start'    => $seq_region_start,
       'end'      => $seq_region_end,
       'strand'   => $seq_region_strand,
       'slice'    => $slice,
       'variation_name' => $variation_name,
       'adaptor'  => $self,
       'dbID'     => $struct_variation_id,
       'source'   => $source_name,
	   'source_description' => $source_description,
	   'class'     => $sv_class,
	   'bound_start' => $bound_start,
	   'bound_end'   => $bound_end,});
  }

  return \@features;


}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$simple_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all simple features in 
               the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub list_dbIDs {
  my $self = shift;
  return $self->_list_dbIDs('structural_variation_feature');
}


=head2 get_all_synonym_sources

    Args[1]     : Bio::EnsEMBL::Variation::StructuralVariation vf
    Example     : my @sources = @{$vf_adaptor->get_all_synonym_sources($vf)};
    Description : returns a list of all the sources for synonyms of this
                  StructuralVariation
    ReturnType  : reference to list of strings
    Exceptions  : none
    Caller      : general
    Status      : At Risk
                : Variation database is under development.

=cut

sub get_all_synonym_sources{
    my $self = shift;
    my $vf = shift;
    my %sources;
    my @sources;

    if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::StructuralVariation')) {
	 throw("Bio::EnsEMBL::Variation::StructuralVariation argument expected");
    }
    
    if (!defined($vf->{'_variation_id'}) && !defined($vf->{'variation'})){
	warning("Not possible to get synonym sources for the StructuralVariation: you need to attach a Variation first");
	return \@sources;
    }
    #get the variation_id
    my $variation_id;
    if (defined ($vf->{'_variation_id'})){
	$variation_id = $vf->{'_variation_id'};
    }
    else{
	$variation_id = $vf->variation->dbID();
    }
    #and go to the varyation_synonym table to get the extra sources
    my $source_name;
    my $sth = $self->prepare(qq{SELECT s.name 
				FROM variation_synonym vs, source s 
				WHERE s.source_id = vs.source_id
			        AND   vs.variation_id = ?
			    });
    $sth->bind_param(1,$variation_id,SQL_INTEGER);
    $sth->execute();
    $sth->bind_columns(\$source_name);
    while ($sth->fetch){
	$sources{$source_name}++;
    }
    @sources = keys(%sources); 
    return \@sources;
}

=head2 fetch_by_name

    Args[1]     : string $name
    Example     : my $structural_variation = $sv_adaptor->fetch_by_name('esv263');
    Description : returns the structural variation with the given variation name (or undef if one isn't found).
                  If the name argument is undef this will be converted to NULL in the SQL statement generated.
    ReturnType  : Bio::EnsEMBL::Variation::StructuralVariation
    Exceptions  : thrown if there are multiple objects found with the same variation name
    Caller      : general
    Status      : At Risk
                : Variation database is under development.

=cut

sub fetch_by_name {
    my ($self, $name) = @_;
    
    my $constraint = sprintf('variation_name = %s', $self->dbc->db_handle->quote( $name, SQL_VARCHAR ) );
    my $objs = $self->generic_fetch($constraint);
    throw("Multiple structural variations found with the same name: '$name'") if @$objs > 1;
    return $objs->[0] if @$objs == 1;
}

1;

=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor
#
# Copyright (c) 2011 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
	$sa   = $reg->get_adaptor("human","core","slice");
  $svfa = $reg->get_adaptor("human","variation","structuralvariationfeature");
  $sva  = $reg->get_adaptor("human","variation","structuralvariation");

  # Get a StructuralVariationFeature by its internal identifier
  $svf = $svfa->fetch_by_dbID(145);

  # get all StructuralVariationFeatures in a region
  $slice = $sa->fetch_by_region('chromosome', 'X', 1e6, 2e6);
  foreach $svf (@{$svfa->fetch_all_by_Slice($slice)}) {
    print $svf->start(), '-', $svf->end(), ' ', $svf->allele_string(), "\n";
  }

  # fetch all genome hits for a particular structural variation
  $sv = $sva->fetch_by_name('esv1285');

  foreach $svf (@{$svfa->fetch_all_by_StructuralVariation($sv)}) {
    print $svf->seq_region_name(), $svf->seq_region_start(), '-',
          $svf->seq_region_end(),"\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity for StructuralVariationFeature objects.
Genomic locations of structural variations can be obtained from the database using this
adaptor. See the base class BaseFeatureAdaptor for more information.
By default, the 'fetch_all_by_...'-methods will not return variations
that have been flagged as failed in the Ensembl QC. This behaviour can be modified
by setting the include_failed_variations flag in Bio::EnsEMBL::Variation::DBSQL::DBAdaptor.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor;

use Bio::EnsEMBL::Variation::StructuralVariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Iterator;

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor', 'Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Example    : my $svfs = $svfa->fetch_all_by_Slice($slice);
  Description: Retrieves all strucutral variation features on the given Slice.
  Returntype : listref of Bio::EnsEMBL::StrucutralVariationFeatures
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my ($self, $slice) = @_;
  return $self->fetch_all_by_Slice_constraint($slice, '');
}


=head2 fetch_all_by_StructuralVariation

  Arg [1]    : Bio::EnsEMBL:Variation::StructuralVariation $var
  Example    : my @svfs = @{$svfa->fetch_all_by_StructuralVariation($var)};
  Description: Retrieves all structural variation features for a given structural variation. Most
               structural variations should only hit the genome once and only a return
               a single structural variation feature.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut


sub fetch_all_by_StructuralVariation {
  my $self = shift;
  my $var  = shift;

  if(!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::StructuralVariation')) {
    throw('Bio::EnsEMBL::Variation::StructuralVariation arg expected');
  }

  if(!defined($var->dbID())) {
    throw("StructuralVariation arg must have defined dbID");
  }

  return $self->generic_fetch("svf.structural_variation_id = ".$var->dbID());
}


=head2 fetch_Iterator_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Description: Returns a listref of strucutral variation features created 
               from the database which are on the Slice defined by $slice 
               and fulfill the SQL constraint defined by $constraint, using the iterator method.
  Returntype : listref of StructuralVariationFeatures
  Exceptions : thrown if $slice is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_Iterator_by_Slice_constraint {
    my ($self, $slice, $constraint) = @_;
    
    $self->{_iterator} = 1;
    
    my $iterator = $self->fetch_all_by_Slice_constraint($slice, $constraint);

    $self->{_iterator} = 0;
    
    return $iterator;
}


=head2 fetch_all_by_Slice_SO_term

  Arg [1]	   : Bio::EnsEMBL::Slice
  Arg [2]    : SO term (string)
  Example    : $slice = $slice_adaptor->fetch_by_region("chromosome", 1, 100000, 200000);
               $SO_term = 'copy_number_variation';
               @svfs = @{$svf_adaptor->fetch_all_by_Slice_SO_term($slice,$SO_term)};
  Description: Retrieves all structural variation features in a slice with a variant type 
	             (structural variation class) or an allele type (supporting structural variation class) 
							 corresponding to the SO term.
  Returntype : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Slice_SO_term {
	my $self    = shift;
	my $slice   = shift;
	my $SO_term = shift;
	
	if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
	
	my $aa = $self->db->get_AttributeAdaptor;
	my $sv_class_id = $aa->attrib_id_for_type_value('SO_term',$SO_term);
	
	if (!defined($sv_class_id)) {
		warn "The SO term '$SO_term' has not been found";
		return [];
	}
	
	my $cols = join ",", $self->_columns();

	my $sth = $self->prepare(qq{
  SELECT DISTINCT $cols
	FROM structural_variation_feature svf, source s, supporting_structural_variation ssv
	WHERE svf.source_id = s.source_id
		AND svf.seq_region_id = ?
		AND svf.seq_region_end > ?
		AND svf.seq_region_start < ?
  	AND svf.structural_variation_id = ssv.structural_variation_id
		AND (svf.class_attrib_id = ? OR ssv.class_attrib_id = ?)
  });
  $sth->execute($slice->get_seq_region_id, $slice->start, $slice->end, $sv_class_id, $sv_class_id);
	
	my $result = $self->_objs_from_sth($sth);
	$sth->finish;
	
	return $result;
}


# method used by superclass to construct SQL
sub _tables { 
	my $self = shift;
    
	my @tables = ( ['structural_variation_feature', 'svf'], [ 'source', 's'] );
		
	return @tables;
}

sub _default_where_clause {
  my $self = shift;

  return 'svf.source_id = s.source_id';
}

sub _columns {
  return qw( svf.structural_variation_feature_id svf.seq_region_id svf.outer_start svf.seq_region_start 
						 svf.inner_start svf.inner_end svf.seq_region_end svf.outer_end svf.seq_region_strand 
						 svf.structural_variation_id svf.variation_name s.name s.version svf.class_attrib_id 
						 svf.allele_string);
}

sub _objs_from_sth {
    my ($self, $sth, $mapper, $dest_slice) = @_;

    #warn $sth->sql;

    # 
    # This code is ugly because an attempt has been made to remove as many
    # function calls as possible for speed purposes.  Thus many caches and
    # a fair bit of gymnastics is used.
    #

    my $sa = $self->db()->dnadb()->get_SliceAdaptor();

    my $aa = $self->db->get_AttributeAdaptor;

    my %slice_hash;
    my %sr_name_hash;
    my %sr_cs_hash;

    my ($structural_variation_feature_id, $seq_region_id, $outer_start, $seq_region_start, $inner_start, 
		    $inner_end, $seq_region_end, $outer_end, $seq_region_strand, $structural_variation_id,
        $variation_name, $source_name, $source_version, $class_attrib_id,$allele_string, $last_svf_id);

    $sth->bind_columns(\$structural_variation_feature_id, \$seq_region_id, \$outer_start, \$seq_region_start, 
		                   \$inner_start, \$inner_end, \$seq_region_end, \$outer_end, \$seq_region_strand, 
											 \$structural_variation_id, \$variation_name, \$source_name, \$source_version, 
											 \$class_attrib_id, \$allele_string);

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

    my $finished = 0;
    
    my $iterator = Bio::EnsEMBL::Utils::Iterator->new(sub{    
        
        return undef if $finished;

        FEATURE: while( $sth->fetch ) {
        
            # Skip if we are getting multiple rows because of the left join to failed variation
            next if (defined($last_svf_id) && $last_svf_id == $structural_variation_feature_id);
            $last_svf_id = $structural_variation_feature_id;
    
            #get the slice object
            my $slice = $slice_hash{"ID:".$seq_region_id};
            if(!$slice) {
                $slice = $sa->fetch_by_seq_region_id($seq_region_id);
                $slice_hash{"ID:".$seq_region_id} = $slice;
                $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
                $sr_cs_hash{$seq_region_id} = $slice->coord_system();
            }
            
            # remap the feature coordinates to another coord system
            # if a mapper was provided
            
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
            return $self->_create_feature_fast('Bio::EnsEMBL::Variation::StructuralVariationFeature',
            #push @features, Bio::EnsEMBL::Variation::StructuralVariationFeature->new_fast(
            #if use new_fast, then do not need "-" infront of key, i.e 'start' => $seq_region_start,
        
               {'outer_start'    => $outer_start,
							  'start'          => $seq_region_start,
								'inner_start'    => $inner_start,
								'inner_end'      => $inner_end,
                'end'            => $seq_region_end,
								'outer_end'      => $outer_end,
                'strand'         => $seq_region_strand,
                'slice'          => $slice,
                'variation_name' => $variation_name,
                'adaptor'        => $self,
                'dbID'           => $structural_variation_feature_id,
                'source'         => $source_name,
                'source_version' => $source_version,
                'structural_variation_id' => $structural_variation_id,
                'class_SO_term'  => $aa->attrib_value_for_id($class_attrib_id),
								'allele_string'  => $allele_string,
               }
            );
        }
        
        unless ($finished) {
            $sth->finish;
            $finished = 1;
        }
        
        return undef;
    });
    
    if ($self->{_iterator}) {
        return $iterator;
    }
    else {
			my $svfs = $iterator->to_arrayref;
      #warn "Got ".scalar(@$vfs). "VFs";
      return $svfs;
    }   
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @svfs = @{$svfa->list_dbIDs()};
  Description: Gets an array of internal ids for all simple features in 
               the current db
  Returntype : list of integers
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub list_dbIDs {
  my $self = shift;
  return $self->_list_dbIDs('structural_variation_feature');
}


1;

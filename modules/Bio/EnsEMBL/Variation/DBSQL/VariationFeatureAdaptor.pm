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
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $vfa = $reg->get_adaptor("human","variation","variationfeature");
  $sa = $reg->get_adaptor("human","core","slice");
  $va = $reg->get_adaptor("human","variation","variation");

  # Get a VariationFeature by its internal identifier
  $vf = $va->fetch_by_dbID(145);

  # Include the variations that have been flagged as failed in the fetch
  $vfa->db->include_failed_variations(1);
  
  # get all VariationFeatures in a region
  $slice = $sa->fetch_by_region('chromosome', 'X', 1e6, 2e6);
  foreach $vf (@{$vfa->fetch_all_by_Slice($slice)}) {
    print $vf->start(), '-', $vf->end(), ' ', $vf->allele_string(), "\n";
  }

  # fetch all genome hits for a particular variation
  $v = $va->fetch_by_name('rs56');

  foreach $vf (@{$vfa->fetch_all_by_Variation($v)}) {
    print $vf->seq_region_name(), $vf->seq_region_start(), '-',
          $vf->seq_region_end(),"\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity for VariationFeature objects.
Genomic locations of variations can be obtained from the database using this
adaptor. See the base class BaseFeatureAdaptor for more information.
By default, the 'fetch_all_by_...'-methods will not return variations
that have been flagged as failed in the Ensembl QC. This behaviour can be modified
by setting the include_failed_variations flag in Bio::EnsEMBL::Variation::DBSQL::DBAdaptor.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;

use Bio::EnsEMBL::Variation::Allele;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Iterator;
use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_validation_code);

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor', 'Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');
our $MAX_VARIATION_SET_ID = 64;

=head2 fetch_all

  Description: Returns a listref of all germline variation features
  Returntype : listref of VariationFeatures
  Status     : At risk

=cut

sub fetch_all {
    my $self = shift;
    my $constraint = 'vf.somatic = 0';
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_somatic

  Description: Returns a listref of all somatic variation features
  Returntype : listref of VariationFeatures
  Status     : At risk

=cut

sub fetch_all_somatic {
    my $self = shift;
    my $constraint = 'vf.somatic = 1';
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Description: Returns a listref of germline variation features created 
               from the database which are on the Slice defined by $slice 
               and fulfill the SQL constraint defined by $constraint.
  Returntype : listref of VariationFeatures
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_by_Slice_constraint {
    my ($self, $slice, $constraint) = @_;
    
    # by default, filter outsomatic mutations
    my $somatic_constraint = 'vf.somatic = 0';
    
    if ($constraint) {
        $constraint .= " AND $somatic_constraint";
    }
    else {
        $constraint = $somatic_constraint;
    }
    
    # Add the constraint for failed variations
    $constraint .= " AND " . $self->db->_exclude_failed_variations_constraint();
    
    return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}

sub fetch_all_by_Slice_constraint_with_Variations {
    my $self = shift;
    $self->{_get_variations} = 1;
    my $vfs = $self->fetch_all_by_Slice_constraint(@_);
    $self->{_get_variations} = 0;
    return $vfs;
}

sub fetch_all_by_Slice_constraint_with_TranscriptVariations {
    my $self = shift;
    $self->{_get_transcript_variations} = 1;
    my $vfs = $self->fetch_all_by_Slice_constraint(@_);
    $self->{_get_transcript_variations} = 0;
    return $vfs;
}

sub fetch_all_somatic_by_Slice_constraint_with_TranscriptVariations {
    my $self = shift;
    $self->{_get_transcript_variations} = 1;
    my $vfs = $self->fetch_all_somatic_by_Slice_constraint(@_);
    $self->{_get_transcript_variations} = 0;
    return $vfs;
}

=head2 fetch_all_somatic_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Description: Returns a listref of somatic variation features created 
               from the database which are on the Slice defined by $slice 
               and fulfill the SQL constraint defined by $constraint.
  Returntype : listref of VariationFeatures
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_somatic_by_Slice_constraint {
    my ($self, $slice, $constraint) = @_;
    
    my $somatic_constraint = 'vf.somatic = 1';
    
    if ($constraint) {
        $constraint .= " AND $somatic_constraint";
    }
    else {
        $constraint = $somatic_constraint;
    }
    
    # Add the constraint for failed variations
    $constraint .= " AND " . $self->db->_exclude_failed_variations_constraint();
    
    return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Example    : my $vfs = $vfa->fetch_all_by_Slice($slice);
  Description: Retrieves all variation features on the given Slice.
               NOTE: only germline variations will be returned, if you want 
               somatic mutations use the fetch_all_somatic_by_Slice method.
  Returntype : listref of Bio::EnsEMBL::VariationFeatures
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my ($self, $slice) = @_;
  return $self->fetch_all_by_Slice_constraint($slice, '');
}

=head2 fetch_all_somatic_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Example    : my $svfs = $vfa->fetch_all_somatic_by_Slice($slice);
  Description: Retrieves a list of variation features representing somatic mutations on the given Slice.
  Returntype : listref of Bio::EnsEMBL::VariationFeatures
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_somatic_by_Slice {
  my ($self, $slice) = @_;
  return $self->fetch_all_somatic_by_Slice_constraint($slice, '');
}

=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL:Variation::Variation $var
  Example    : my @vfs = @{$vfa->fetch_all_by_Variation($var)};
  Description: Retrieves all variation features for a given variation.  Most
               variations should only hit the genome once and only a return
               a single variation feature.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut


sub fetch_all_by_Variation {
  my $self = shift;
  my $var  = shift;

  if(!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation arg expected');
  }

  if(!defined($var->dbID())) {
    throw("Variation arg must have defined dbID");
  }

  return $self->generic_fetch("vf.variation_id = ".$var->dbID());
}

=head2 fetch_all_genotyped_by_Slice

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Example    : my @vfs = @{$vfa->fetch_all_genotyped_by_Slice($slice)};
  Description: Retrieves all variation features that have been gentoyped for a given slice.
               Most variations should only hit the genome once and only a return
               a single variation feature.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_genotyped_by_Slice{
    my $self = shift;
    my $slice = shift;

    my $constraint = "vf.flags & 1 AND vf.somatic = 0";
    #call the method fetch_all_by_Slice_constraint with the genotyped constraint
    return $self->fetch_all_by_Slice_constraint($slice,$constraint);
}

sub _internal_fetch_all_with_annotation_by_Slice{

	my $self = shift;
	my $slice = shift;
	my $v_source = shift;
	my $p_source = shift;
	my $annotation = shift;
	my $constraint = shift;
	
	if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
		throw('Bio::EnsEMBL::Slice arg expected');
	}
	
	my $extra_sql = '';
    
    if(defined $v_source) {
		$extra_sql .= qq{ AND s.name = '$v_source' };
    }
    
    if(defined $p_source) {
		$extra_sql .= qq{ AND ps.name = '$p_source' };
    }
    
    if(defined $annotation) {
		if($annotation =~ /^[0-9]+$/) {
		  $extra_sql .= qq{ AND p.phenotype_id = $annotation };
		}
		else {
		  $extra_sql .= qq{ AND (p.name = '$annotation' OR p.description LIKE '%$annotation%') };
		}
    }
    
    if ($constraint) {
        $extra_sql .= qq{ AND $constraint };
    }
    
    # Add the constraint for failed variations
    $extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();
    
    my $cols = join ",", $self->_columns();
    
    my $sth = $self->prepare(qq{
		SELECT $cols
		FROM (variation_feature vf, variation_annotation va,
		phenotype p, source s, source ps, study st) # need to link twice to source
		LEFT JOIN failed_variation fv ON (fv.variation_id = vf.variation_id)
		WHERE va.study_id = st.study_id
		AND st.source_id = ps.source_id
		AND vf.source_id = s.source_id
		AND vf.variation_id = va.variation_id
		AND va.phenotype_id = p.phenotype_id
		$extra_sql
		AND vf.seq_region_id = ?
		AND vf.seq_region_end > ?
		AND vf.seq_region_start < ?
		GROUP BY vf.variation_feature_id
    });
    
    $sth->execute($slice->get_seq_region_id, $slice->start, $slice->end);
    
    return $self->_objs_from_sth($sth, undef, $slice);
}

=head2 fetch_all_with_annotation_by_Slice

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Arg [2]    : $variation_feature_source [optional]
  Arg [3]    : $annotation_source [optional]
  Arg [4]    : $annotation_name [optional]
  Example    : my @vfs = @{$vfa->fetch_all_with_annotation_by_Slice($slice)};
  Description: Retrieves all germline variation features associated with annotations for
               a given slice.
               The optional $variation_feature_source argument can be used to
               retrieve only variation features from a paricular source.
               The optional $annotation source argument can be used to
               retrieve only variation features with annotations provided by
               a particular source.
               The optional $annotation_name argument can
               be used to retrieve only variation features associated with
               that annotation - this can also be a phenotype's dbID.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_with_annotation_by_Slice {
    my $self = shift;
    my ($slice, $v_source, $p_source, $annotation) = @_;
    my $constraint = 'vf.somatic = 0';
    return $self->_internal_fetch_all_with_annotation_by_Slice($slice, $v_source, $p_source, $annotation, $constraint);
}

=head2 fetch_all_somatic_with_annotation_by_Slice

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Arg [2]    : $variation_feature_source [optional]
  Arg [3]    : $annotation_source [optional]
  Arg [4]    : $annotation_name [optional]
  Example    : my @vfs = @{$vfa->fetch_all_somatic_with_annotation_by_Slice($slice)};
  Description: Retrieves all somatic variation features associated with annotations for
               a given slice.
               (see fetch_all_with_annotation_by_Slice documentation for description of
               the other parameters)
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_somatic_with_annotation_by_Slice {
    my $self = shift;
    my ($slice, $v_source, $p_source, $annotation) = @_;
    my $constraint = 'vf.somatic = 1';
    return $self->_internal_fetch_all_with_annotation_by_Slice($slice, $v_source, $p_source, $annotation, $constraint);
}

=head2 fetch_all_by_Slice_VariationSet

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Arg [2]    : Bio::EnsEMBL:Variation::VariationSet $set
  Example    : my @vfs =
@{$vfa->fetch_all_by_Slice_VariationSet($slice, $set)};
  Description: Retrieves all variation features in a slice that belong to a
			   given variation set.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_VariationSet{

  my $self = shift;
  my $slice = shift;
  my $set = shift;
	
  #$self->{_get_variations} = 1;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  if(!ref($set) || !$set->isa('Bio::EnsEMBL::Variation::VariationSet')) {
    throw('Bio::EnsEMBL::Variation::VariationSet arg expected');
  }
  
  #�Get the bitvalue for this set and its subsets
  my $bitvalue = $set->_get_bitvalue();
  
  # Add a constraint to only return VariationFeatures having the primary keys of the supplied VariationSet or its subsets in the variation_set_id column
  my $constraint = " vf.variation_set_id & $bitvalue ";
  
  #�Get the VariationFeatures by calling fetch_all_by_Slice_constraint
  my $vfs = $self->fetch_all_by_Slice_constraint($slice,$constraint);
  

  $self->{_get_variations} = 0;

  return $vfs;
}


=head2 fetch_all_by_Slice_Population

  Arg [1]	 : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::Variation::Population
  Arg [3]	 : $minimum_frequency (optional)
  Example    : $pop = $pop_adaptor->fetch_by_dbID(659);
			  $slice = $slice_adaptor->fetch_by_region("chromosome", 1, 1, 1e6);
              @vfs = @{$vf_adaptor->fetch_all_by_Slice_Population($pop,$slice)};
  Description: Retrieves all variation features in a slice which are stored for
			   a specified population. If $minimum_frequency is supplied, only
			   variations with a minor allele frequency (MAF) greater than
			   $minimum_frequency will be returned.
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Slice_Population {
  my $self = shift;
  
  my $slice = shift;
  
  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  
  my $pop = shift;
  
  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
	throw('Bio::EnsEMBL::Variation::Population arg expected');
  }
  
  # default to 5% frequency
  my $freq = shift;
  my $extra_sql = '';
  
  if(defined $freq) {
	
	# adjust frequency if given a percentage
	$freq /= 100 if $freq > 1;
	$extra_sql = qq{ AND (IF(a.frequency > 0.5, 1-a.frequency, a.frequency) > $freq) }
  }
  
  # Add the constraint for failed variations
  $extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();

  my $cols = join ",", $self->_columns();
  
  my $sth = $self->prepare(qq{
	SELECT $cols
	FROM (variation_feature vf, source s, allele a)
	LEFT JOIN failed_variation fv ON (fv.variation_id = vf.variation_id)
	WHERE vf.source_id = s.source_id
	AND vf.variation_id = a.variation_id
	AND a.sample_id = ?
	$extra_sql
	AND vf.seq_region_id = ?
	AND vf.seq_region_end >= ?
	AND vf.seq_region_start <= ?
    GROUP BY a.variation_id
  });
  
  $sth->execute($pop->dbID, $slice->get_seq_region_id, $slice->start, $slice->end);

##### ALTERNATIVE USING SUB-QUERY - WORKS FOR MULTIALLELIC BUT TAKES FOREVER!!!
#  my $cols = join ",", $self->_columns();
#  
#  my $sth = $self->prepare(qq{
#	SELECT $cols
#	FROM variation_feature vf, source s, allele a_outer
#	WHERE
#	  vf.source_id = s.source_id AND
#	  vf.variation_id = a_outer.variation_id AND
#	  a_outer.sample_id = ? AND
#	  a_outer.frequency > ? AND
#	  a_outer.frequency < 0.5 AND
#	  vf.seq_region_id = ? AND
#	  vf.seq_region_end >= ? AND
#	  vf.seq_region_start <= ? AND
#	  NOT EXISTS (
#		SELECT
#		  *
#		FROM
#		  allele a_inner
#		WHERE
#		  a_inner.variation_id = a_outer.variation_id AND
#		  a_inner.sample_id = a_outer.sample_id AND
#		  a_inner.frequency < a_outer.frequency
#	  );
#  });
  
  #$sth->execute($pop->dbID, $freq, $slice->get_seq_region_id, $slice->start, $slice->end);
  
  return $self->_objs_from_sth($sth, undef, $slice);
}

sub _internal_fetch_all_with_annotation {
    
    my ($self, $v_source, $p_source, $annotation, $constraint) = @_;
    
    my $extra_sql = '';
    
    if(defined $v_source) {
        $extra_sql .= qq{ AND s.name = '$v_source' };
    }
    
    if(defined $p_source) {
        $extra_sql .= qq{ AND ps.name = '$p_source' };
    }
    
    if(defined $annotation) {
        if($annotation =~ /^[0-9]+$/) {
          $extra_sql .= qq{ AND p.phenotype_id = $annotation };
        }
        else {
          $extra_sql .= qq{ AND (p.name = '$annotation' OR p.description LIKE '%$annotation%') };
        }
    }
    
    if ($constraint) {
        $extra_sql .= qq{ AND $constraint };
    }
    
    # Add the constraint for failed variations
    $extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();
    
    my $cols = join ",", $self->_columns();
    
    my $sth = $self->prepare(qq{
        SELECT $cols
        FROM (variation_feature vf, variation_annotation va,
        phenotype p, source s, source ps, study st) # need to link twice to source
	LEFT JOIN failed_variation fv ON (fv.variation_id = vf.variation_id)
        WHERE va.study_id = st.study_id
				AND st.source_id = ps.source_id 
        AND vf.source_id = s.source_id
        AND vf.variation_id = va.variation_id
        AND va.phenotype_id = p.phenotype_id
        $extra_sql
        GROUP BY vf.variation_feature_id
    });
    
    $sth->execute;
    
    return $self->_objs_from_sth($sth);
}

=head2 fetch_all_with_annotation

  Arg [1]    : $variation_feature_source [optional]
  Arg [2]    : $annotation_source [optional]
  Arg [3]    : $annotation_name [optional]
  Example    : my @vfs = @{$vfa->fetch_all_with_annotation('EGA', undef, 123)};
  Description: Retrieves all germline variation features associated with the given annotation
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Caller     : webcode
  Status     : Experimental

=cut

sub fetch_all_with_annotation {
    
    my ($self, $v_source, $p_source, $annotation, $constraint) = @_;
    
    my $somatic_constraint = 'vf.somatic = 0';
    
    if ($constraint) {
        $constraint .= " AND $somatic_constraint";
    }
    else {
        $constraint = $somatic_constraint;
    }
    
    return $self->_internal_fetch_all_with_annotation($v_source, $p_source, $annotation, $constraint);
}

=head2 fetch_all_somatic_with_annotation

  Arg [1]    : $variation_feature_source [optional]
  Arg [2]    : $annotation_source [optional]
  Arg [3]    : $annotation_name [optional]
  Example    : my @vfs = @{$vfa->fetch_all_somatic_with_annotation('COSMIC', undef, 807)};
  Description: Retrieves all somatic variation features associated with the given annotation
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Caller     : webcode
  Status     : Experimental

=cut

sub fetch_all_somatic_with_annotation {
    
    my ($self, $v_source, $p_source, $annotation, $constraint) = @_;
    
    my $somatic_constraint = 'vf.somatic = 1';
    
    if ($constraint) {
        $constraint .= " AND $somatic_constraint";
    }
    else {
        $constraint = $somatic_constraint;
    }
    
    return $self->_internal_fetch_all_with_annotation($v_source, $p_source, $annotation, $constraint);
}

sub fetch_Iterator_by_Slice_constraint {
    my ($self, $slice, $constraint) = @_;
    
    $self->{_iterator} = 1;
    
    my $iterator = $self->fetch_all_by_Slice_constraint($slice, $constraint);

    $self->{_iterator} = 0;
    
    return $iterator;
}

# method used by superclass to construct SQL
sub _tables { 
    my $self = shift;
    
    my @tables = (
        ['variation_feature', 'vf'],
		[ 'source', 's']
	);
	
	#�If we are including failed_variations, add that table
	push(@tables,['failed_variation', 'fv']) unless ($self->db->include_failed_variations());
	
	return @tables;
}

#�Add a left join to the failed_variation table
sub _left_join { 
    my $self = shift;
    
    # If we are including failed variations, skip the left join
    return () if ($self->db->include_failed_variations());
    return ([ 'failed_variation', 'fv.variation_id = vf.variation_id']); 
}

sub _default_where_clause {
  my $self = shift;

  return 'vf.source_id = s.source_id';
}

sub _columns {
  return qw( vf.variation_feature_id vf.seq_region_id vf.seq_region_start
             vf.seq_region_end vf.seq_region_strand vf.variation_id
             vf.allele_string vf.variation_name vf.map_weight s.name s.version vf.somatic 
             vf.validation_status vf.consequence_type vf.class_attrib_id);
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

    my @features;
    my %slice_hash;
    my %sr_name_hash;
    my %sr_cs_hash;

    my ($variation_feature_id, $seq_region_id, $seq_region_start,
      $seq_region_end, $seq_region_strand, $variation_id,
      $allele_string, $variation_name, $map_weight, $source_name, $source_version,
      $is_somatic, $validation_status, $consequence_type, $class_attrib_id, $last_vf_id);

    $sth->bind_columns(\$variation_feature_id, \$seq_region_id,
                     \$seq_region_start, \$seq_region_end, \$seq_region_strand,
                     \$variation_id, \$allele_string, \$variation_name,
                     \$map_weight, \$source_name, \$source_version, \$is_somatic, \$validation_status, 
                     \$consequence_type, \$class_attrib_id);

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
            next if (defined($last_vf_id) && $last_vf_id == $variation_feature_id);
            $last_vf_id = $variation_feature_id;
    
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
            my $validation_code;
            if (defined($validation_status)) {
                $validation_code = get_validation_code([split(',',$validation_status)]);
            }
            
            #my $overlap_consequences = $self->_variation_feature_consequences_for_set_number($consequence_type);
            
            my $overlap_consequences = [ map { $OVERLAP_CONSEQUENCES{$_} } split /,/, $consequence_type ];

            # consequence_type
            return $self->_create_feature_fast('Bio::EnsEMBL::Variation::VariationFeature',
            #push @features, Bio::EnsEMBL::Variation::VariationFeature->new_fast(
            #if use new_fast, then do not need "-" infront of key, i.e 'start' => $seq_region_start,
        
                {'start'    => $seq_region_start,
                'end'      => $seq_region_end,
                'strand'   => $seq_region_strand,
                'slice'    => $slice,
                'allele_string' => $allele_string,
                'variation_name' => $variation_name,
                'adaptor'  => $self,
                'dbID'     => $variation_feature_id,
                'map_weight' => $map_weight,
                'source'   => $source_name,
                'source_version' => $source_version,
                'is_somatic' => $is_somatic,
                'validation_code' => $validation_code,
                'overlap_consequences' => $overlap_consequences,
                '_variation_id' => $variation_id,
                'class_SO_term' => $aa->attrib_value_for_id($class_attrib_id),
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
        if ($self->{_get_variations}) {
            my $vfs = $iterator->to_arrayref;
            my @v_ids = map { $_->{_variation_id} } @$vfs;
            my $vs = $self->db->get_VariationAdaptor->fetch_all_by_dbID_list(\@v_ids);
            my %vs_by_id = map { $_->dbID => $_ } @$vs;
            #warn "Got variations";
            map { $_->variation( $vs_by_id{ $_->{_variation_id} }) } @$vfs;
            return $vfs;
        }
        if ($self->{_get_transcript_variations}) {
            my $vfs = $iterator->to_arrayref;
            return $vfs unless @$vfs;
            #warn "getting transcript variations";
            my $tvs = $self->db->get_TranscriptVariationAdaptor->fetch_all_by_VariationFeatures($vfs);
            for my $tv (@$tvs) {
                $tv->variation_feature->add_TranscriptVariation($tv);
                #$tv->variation_feature->{transcript_variations}->{$tv->transcript_stable_id} = $tv;
            }
            return $vfs;
        }
        else {
            my $vfs = $iterator->to_arrayref;
            #warn "Got ".scalar(@$vfs). "VFs";
            return $vfs;
        }
    }   
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
  return $self->_list_dbIDs('variation_feature');
}


=head2 get_all_synonym_sources

    Args[1]     : Bio::EnsEMBL::Variation::VariationFeature vf
    Example     : my @sources = @{$vf_adaptor->get_all_synonym_sources($vf)};
    Description : returns a list of all the sources for synonyms of this
                  VariationFeature
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

    if(!ref($vf) || !$vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
	 throw("Bio::EnsEMBL::Variation::VariationFeature argument expected");
    }
    
    if (!defined($vf->{'_variation_id'}) && !defined($vf->{'variation'})){
	warning("Not possible to get synonym sources for the VariationFeature: you need to attach a Variation first");
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

=head2 new_fake

  Arg [1]    : string $species
  Example    :
	$vfa = Bio::EnsEMBL::Variation::VariationFeatureAdaptor->new_fake('human');
  Description: Creates a VariationFeatureAdaptor with no underlying database
			   attached. Should be used only when getting consequence types for
			   species with no variation database available.
  Returntype : Bio::EnsEMBL::Variation::VariationFeatureAdaptor
  Exceptions : throw if no species given
  Caller     : called from webcode for species where no variation database present
  Status     : Stable

=cut

sub new_fake {
  my $class = shift;
  my $species = shift;
  
  throw("No species defined") unless defined $species;
  
  my $self = bless {}, $class;
  
  $self->{'species'} = $species;
  
  return $self;
}

sub _parse_hgvs_position {
    my $description = shift;
    
    my ($start,$start_offset,$end,$end_offset) = $description =~ m/^([\-\*]?\d+)((?:[\+\-]\d+)?)(?:_([\-\*]?\d+)((?:[\+\-]\d+)?))?/;
    
    #�Fill in missing values
    $start_offset = 0 unless (defined($start_offset) && length($start_offset));
    unless (defined($end) && length($end)) {
        $end = $start;
        $end_offset = $start_offset;
    }
    $end_offset = 0 unless (defined($end_offset) && length($end_offset));
    	
    #�Get rid of any '+' signs that may have been parsed into the variables
    $start_offset = substr($start_offset,1) if (substr($start_offset,0,1) eq '+');
    $end_offset = substr($end_offset,1) if (substr($end_offset,0,1) eq '+');
    
    return [$start,$end,$start_offset,$end_offset];
}

=head2 fetch_by_hgvs_notation

    Arg[1]      : String $hgvs
    Example     : my $hgvs = 'LRG_8t1:c.3891A>T';
		  $vf = $vf_adaptor->fetch_by_hgvs_notation($hgvs);
    Description : Parses an HGVS notation and tries to create a VariationFeature object
		  based on the notation. The object will have a Variation
		  and Alleles attached.
    ReturnType  : Bio::EnsEMBL::Variation::VariationFeature, undef on failure
    Exceptions  : thrown on error
    Caller      : general
    Status      : Stable

=cut

sub fetch_by_hgvs_notation {
    my $self = shift;
    my $hgvs = shift;
    my $user_slice_adaptor = shift;
    my $user_transcript_adaptor = shift;
    
    #�Split the HGVS notation into the reference, notation type and variation description
    my ($reference,$type,$description) = $hgvs =~ m/^([^\:]+)\:.*?([cgmrp]?)\.?([\*\-0-9]+.*)$/i;
    
    my $extra;
    
    if($description =~ m/\(.+\)/) {        
        ($description, $extra) = $description =~ /(.+?)(\(.+\))/;        
        throw ("Could not parse the HGVS notation $hgvs - can't interpret \'$extra\'") unless $extra eq '(p.=)';
    }
    
    #�If any of the fields are unknown, return undef
    throw ("Could not parse the HGVS notation $hgvs") unless (defined($reference) && defined($type) && defined($description));
    
    # strip version number from reference
    $reference =~ s/\.\d+//g if $reference =~ /^ENS/;
    
    my ($start,$end,$start_offset,$end_offset,$ref_allele,$alt_allele,$class);
    
    #�Parse differently depending on the type of the notation and the type of the variation
    if ($type =~ m/[gcrm]/i) {
	
	    # Get the positions
	    ($start,$end,$start_offset,$end_offset) = @{_parse_hgvs_position($description)};
	    
    	# A single nt substitution, reference and alternative alleles are required
    	if ($description =~ m/>/) {
    	    $class = 'snv';
    	    ($ref_allele,$alt_allele) = $description =~ m/([A-Z]+)>([A-Z]+)$/i;
    	}
    	#�A delins, the reference allele is optional
    	elsif ($description =~ m/del.*ins/i) {
    	    $class = 'delins';
    	    ($ref_allele,$alt_allele) = $description =~ m/del(.*?)ins([A-Z]+)$/i;    	    
    	}
    	# A deletion, the reference allele is optional
    	elsif ($description =~ m/del/i) {
    	    $class = 'del';
    	    ($ref_allele) = $description =~ m/del([A-Z]*)$/i; 
    	}
    	# A duplication, the reference allele is optional
    	elsif ($description =~ m/dup/i) {
    	    $class = 'dup';
    	    ($ref_allele) = $description =~ m/dup([A-Z]*)$/i;
    	}
    	# An inversion, the reference allele is optional
    	elsif ($description =~ m/inv/i) {
    	    $class = 'inv';
    	    ($ref_allele) = $description =~ m/inv([A-Z]*)$/i;
    	}
    	else {
    	    $class = 'unknown';
    	    die ("The variant class for HGVS notation '$hgvs' is unknown or could not be correctly recognized");
    	}
    	
    	# If the reference allele was omitted, set it to undef
    	$ref_allele = undef unless (defined($ref_allele) && length($ref_allele));
    	
    	#�If the position given was in a intronic or UTR position, it could be undefined what reference sequence the position actually refers to. Issue a warning that we will use the Ensembl reference sequence.
    	warn ("The position specified by HGVS notation '$hgvs' refers to a nucleotide that may not have a specific reference sequence. The current Ensembl genome reference sequence will be used.") if ($start_offset || $end_offset || substr($start,0,1) eq '*' || substr($end,0,1) eq '*' || $start < 0);
    	
    }
    # Else, it is protein notation
    else {
	throw ("Parsing of HGVS protein notation has not yet been implemented");
    }
    
    #�Get a slice representing this variation
    my $slice_adaptor = $user_slice_adaptor || $self->db()->dnadb()->get_SliceAdaptor();
    my $slice;
    my $strand;
    if ($type =~ m/c/i) {
	
    	#�A small fix in case the reference is a LRG and there is no underscore between name and transcript
    	$reference =~ s/^(LRG_[0-9]+)_?(t[0-9]+)$/$1\_$2/i;
    	
    	#�Get the Transcript object for this variation
    	my $transcript_adaptor = $user_transcript_adaptor || $self->db()->dnadb()->get_TranscriptAdaptor();
    	my $transcript = $transcript_adaptor->fetch_by_stable_id($reference) or throw ("Could not get a Transcript object for '$reference'");
    	
    	# Get the TranscriptMapper
    	my $tr_mapper = $transcript->get_TranscriptMapper();
    	
    	#�In case we have a position in the 3' UTR, we need to convert the coordinates by setting them to be the stop codon position + the UTR offset
    	$start = ($transcript->cdna_coding_end() - $transcript->cdna_coding_start() + 1) + int(substr($start,1)) if (substr($start,0,1) eq '*');
    	$end = ($transcript->cdna_coding_end() - $transcript->cdna_coding_start() + 1) + int(substr($end,1)) if (substr($end,0,1) eq '*');
    	
    	#�The mapper can only convert cDNA coordinates, but we have CDS (relative to the start codon), so we need to convert them
    	my ($cds_start,$cds_end) = (($start + $transcript->cdna_coding_start() - ($start > 0)),($end + $transcript->cdna_coding_start() - ($end > 0)));
    	
    	# Convert the cDNA coordinates to genomic coordinates.
    	my @coords = $tr_mapper->cdna2genomic($cds_start,$cds_end);
    	
    	#�Throw an error if we didn't get an unambiguous coordinate back
    	throw ("Unable to map the cDNA coordinates $start\-$end to genomic coordinates for Transcript $reference") if (scalar(@coords) != 1 || !$coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate'));
    	
    	#�Adjust any start and end offsets that resulted from e.g. intronic offsets or UTR positions
    	$strand = $coords[0]->strand();
    	$start = $coords[0]->start() + ($strand >= 0 ? $start_offset : $end_offset);
    	$end = $coords[0]->end() + ($strand < 0 ? $start_offset : $end_offset);
    	
    	#�Get a slice for this variation
    	$slice = $slice_adaptor->fetch_by_region($transcript->coord_system_name(),$transcript->seq_region_name());
    }
    
    #�Get the reference allele based on the coordinates
    my $refseq_allele = $slice->subseq($start,$end,$strand);
    
    # If the reference from the sequence does not correspond to the reference given in the HGVS notation, throw an exception 
    throw ("Reference allele extracted from $reference:$start-$end ($refseq_allele) does not match reference allele given by HGVS notation ($ref_allele)") if (defined($ref_allele) && $ref_allele ne $refseq_allele);
    
    # Use the reference allele from the sequence if none was specified in the notation
    $ref_allele ||= $refseq_allele;
    
    # If the variation type is an inversion, the alt allele is the reverse complement of the ref_allele
    if ($class eq 'inv') {
    	$alt_allele = $ref_allele;
    	reverse_comp(\$alt_allele);
    }
    #�Else, if we it is a duplication, set it to be a repeat of the reference allele
    elsif ($class eq 'dup') {
        my $repeat = 2;
        $alt_allele = ${ref_allele}x$repeat; 
    }
    #�Else if we have a deletion, the alt allele should be set to '-'
    elsif ($class eq 'del') {
        $alt_allele = '-';
    }
    
    #�Create Allele objects
    my @allele_objs;
    foreach my $allele ($ref_allele,$alt_allele) {
	push(@allele_objs,Bio::EnsEMBL::Variation::Allele->new('-adaptor' => $self, '-allele' => $allele));
    }
    
    #�Create a variation object. Use the HGVS string as its name
    my $variation = Bio::EnsEMBL::Variation::Variation->new(
	'-adaptor' => $self->db()->get_VariationAdaptor(),
	'-name' => $hgvs,
	'-source' => 'Parsed from HGVS notation',
	'-alleles' => \@allele_objs
    );
    
    #�Create a variation feature object
    my $variation_feature = Bio::EnsEMBL::Variation::VariationFeature->new(
	'-adaptor' => $self,
	'-start' => $start,
	'-end' => $end,
	'-strand' => $strand,
	'-slice' => $slice,
	'-map_weight' => 1,
	'-variation' => $variation,
	'-allele_string' => "$ref_allele/$alt_allele"
    );
    
    return $variation_feature;
}

1;

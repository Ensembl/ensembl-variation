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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor
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
By default, the 'fetch_all_by_...'-methods will not return structural variants
that have been flagged as failed in the Ensembl QC. This behaviour can be modified
by setting the include_failed_variations flag in Bio::EnsEMBL::Variation::DBSQL::DBAdaptor.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor;

use Bio::EnsEMBL::Variation::StructuralVariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Variation::Utils::Constants qw(%VARIATION_CLASSES);
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Iterator;

our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor', 'Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice 
               $slice the slice from which to obtain features
  Arg [2]    : int $include_supporting_evidence [optional]
  Example    : my $svfs = $svfa->fetch_all_by_Slice($slice);
  Description: Retrieves all germline structural variation features on the given Slice.
               If $include_supporting_evidence is set (i.e. $include_supporting_evidence=1), structural variation features from 
               both structural variation (SV) and their supporting structural variations (SSV) will be 
               returned. By default, it only returns features from structural variations (SV).
  Returntype : reference to list Bio::EnsEMBL::StructuralVariationFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my ($self, $slice, $include_supporting_evidence) = @_;
  return $self->fetch_all_by_Slice_constraint($slice, undef, $include_supporting_evidence);
}

=head2 fetch_all_by_Slice_size_range

  Arg [1]    : Bio::EnsEMBL::Slice 
               $slice the slice from which to obtain features
  Arg [2]    : int $size_min lower value for the size range
  Arg [3]    : int $size_max upper value for the size range 
  Arg [4]    : int $include_supporting_evidence [optional]
  Example    : my $svfs = $svfa->fetch_all_somatic_by_Slice($slice);
  Description: Retrieves all structural variation features on the given Slice for the given size range.
               If $include_supporting_evidence is set (i.e. $include_supporting_evidence=1), structural variation features from 
               both structural variation (SV) and their supporting structural variations (SSV) will be 
               returned. By default, it only returns features from structural variations (SV). 
  Returntype : reference to list Bio::EnsEMBL::StructuralVariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_size_range {
  my ($self, $slice, $size_min, $size_max, $include_supporting_evidence) = @_;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }

  my $constraint = qq{svf.seq_region_end - svf.seq_region_start >= $size_min};
     $constraint .= qq{ AND svf.seq_region_end - svf.seq_region_start < $size_max } if (defined $size_max);

  return $self->fetch_all_by_Slice_constraint($slice, $constraint, $include_supporting_evidence);
}

=head2 fetch_all_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice 
               $slice the slice from which to obtain features
  Arg [2]    : string $constraint [optional]
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [3]    : int $include_supporting_evidence [optional]
  Example    : my $svfs = $svfa->fetch_all_by_Slice($slice);
  Description: Retrieves all germline structural variation features on the given Slice, with the additional constraint(s).
               If $include_supporting_evidence is set (i.e. $include_supporting_evidence=1), structural variation features from 
               both structural variation (SV) and their supporting structural variations (SSV) will be 
               returned. By default, it only returns features from structural variations (SV).
  Returntype : reference to list Bio::EnsEMBL::StructuralVariationFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_constraint {
  my ($self, $slice, $constraint, $include_supporting_evidence) = @_;
  
  my $and_flag = (defined($constraint)) ? undef : 1;
  $constraint .= $self->_internal_exclude_failed_constraint('',$and_flag);
  $constraint .= ' AND ' if ($constraint);
  $constraint .= ' svf.somatic=0 ';
  $constraint .= ' AND '.$self->_internal_exclude_cnv_probe();
  $constraint .= ' AND svf.is_evidence=0 ' if (!$include_supporting_evidence);
  
  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}


=head2 fetch_all_somatic_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice 
               $slice the slice from which to obtain features
  Arg [2]    : int $include_supporting_evidence [optional]
  Example    : my $svfs = $svfa->fetch_all_somatic_by_Slice($slice);
  Description: Retrieves all somatic structural variation features on the given Slice.
               If $include_supporting_evidence is set (i.e. $include_supporting_evidence=1), structural variation features from 
               both structural variation (SV) and their supporting structural variations (SSV) will be 
               returned. By default, it only returns features from structural variations (SV). 
  Returntype : reference to list Bio::EnsEMBL::StructuralVariationFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_somatic_by_Slice {
  my ($self, $slice, $include_supporting_evidence) = @_;
  return $self->fetch_all_somatic_by_Slice_constraint($slice, undef, $include_supporting_evidence);
}

=head2 fetch_all_somatic_by_Slice_size_range

  Arg [1]    : Bio::EnsEMBL::Slice 
               $slice the slice from which to obtain features
  Arg [2]    : int $size_min lower value for the size range
  Arg [3]    : int $size_max upper value for the size range 
  Arg [4]    : int $include_supporting_evidence [optional]
  Example    : my $svfs = $svfa->fetch_all_somatic_by_Slice($slice);
  Description: Retrieves all somatic structural variation features on the given Slice for the given size range.
               If $include_supporting_evidence is set (i.e. $include_supporting_evidence=1), structural variation features from 
               both structural variation (SV) and their supporting structural variations (SSV) will be 
               returned. By default, it only returns features from structural variations (SV). 
  Returntype : reference to list Bio::EnsEMBL::StructuralVariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_somatic_by_Slice_size_range {
  my ($self, $slice, $size_min, $size_max, $include_supporting_evidence) = @_;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }

  my $constraint = qq{svf.seq_region_end - svf.seq_region_start >= $size_min};
     $constraint .= qq{ AND svf.seq_region_end - svf.seq_region_start < $size_max } if (defined $size_max);

  return $self->fetch_all_somatic_by_Slice_constraint($slice, $constraint, $include_supporting_evidence);
}

=head2 fetch_all_somatic_by_Slice_Source

  Arg [1]    : Bio::EnsEMBL::Slice 
               $slice the slice from which to obtain features
  Arg [2]    : Bio::EnsEMBL::Variation::Source $source only return somatic mutations for the given source 
  Arg [3]    : int $include_supporting_evidence [optional]
  Example    : my $svfs = $svfa->fetch_all_somatic_by_Slice($slice);
  Description: Retrieves all somatic structural variation features on the given Slice for the given source.
               If $include_supporting_evidence is set (i.e. $include_supporting_evidence=1), structural variation features from 
               both structural variation (SV) and their supporting structural variations (SSV) will be 
               returned. By default, it only returns features from structural variations (SV). 
  Returntype : reference to list Bio::EnsEMBL::StructuralVariationFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub fetch_all_somatic_by_Slice_Source {
  my $self  = shift;
  my $slice = shift;
  my $source = shift;
  my $include_supporting_evidence = shift;
  
  if (!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  if (!ref($source) || !$source->isa('Bio::EnsEMBL::Variation::Source')) {
    throw('Bio::EnsEMBL::Variation::Source arg expected');
  }
  
  my $constraint = $self->_internal_exclude_failed_constraint("svf.source_id = " . $source->dbID);
  
  my $svfs = $self->fetch_all_somatic_by_Slice_constraint($slice, $constraint, $include_supporting_evidence);
  return $svfs;
}


=head2 fetch_all_somatic_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice 
               $slice the slice from which to obtain features
  Arg [2]    : string $constraint [optional]
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [3]    : int $include_supporting_evidence [optional]
  Example    : my $svfs = $svfa->fetch_all_somatic_by_Slice_constraint($slice,$constraint);
  Description: Retrieves all somatic structural variation features on the given Slice, with the additional constraint(s).
               If $include_supporting_evidence is set (i.e. $include_supporting_evidence=1), structural variation features from 
               both structural variation (SV) and their supporting structural variations (SSV) will be 
               returned. By default, it only returns features from structural variations (SV). 
  Returntype : reference to list Bio::EnsEMBL::StructuralVariationFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_somatic_by_Slice_constraint {
  my ($self, $slice, $constraint, $include_supporting_evidence) = @_;
  
  my $and_flag = (defined($constraint)) ? undef : 1;
  $constraint .= $self->_internal_exclude_failed_constraint('',$and_flag);
  $constraint .= ' AND ' if ($constraint && $constraint !~/\s*AND\s*$/g);
  $constraint .= ' svf.somatic=1 ';
  $constraint .= ' AND svf.is_evidence=0 ' if (!$include_supporting_evidence);
  
  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}


=head2 fetch_all_by_StructuralVariation

  Arg [1]    : Bio::EnsEMBL:Variation::StructuralVariation or 
               Bio::EnsEMBL::Variation::SupportingStructuralVariation $var
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

  if(!ref($var) || (!$var->isa('Bio::EnsEMBL::Variation::StructuralVariation') &&
                    !$var->isa('Bio::EnsEMBL::Variation::SupportingStructuralVariation'))
  ) {
    throw('Bio::EnsEMBL::Variation::StructuralVariation or 
           Bio::EnsEMBL::Variation::SupportingStructuralVariation arg expected');
  }

  if(!defined($var->dbID())) {
    throw("StructuralVariation arg must have defined dbID");
  }
  
  my $constraint = $self->_internal_exclude_failed_constraint("svf.structural_variation_id = ".$var->dbID());
  
  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Slice_SO_term

  Arg [1]     : Bio::EnsEMBL::Slice
  Arg [2]    : SO term (string)
  Example    : $slice = $slice_adaptor->fetch_by_region("chromosome", 1, 100000, 200000);
               $SO_term = 'copy_number_variation';
               @svfs = @{$svf_adaptor->fetch_all_by_Slice_SO_term($slice,$SO_term)};
  Description: Retrieves all structural variation features in a slice with a variant type 
               (structural variation class) or an allele type (supporting structural variation class) 
               corresponding to the SO term.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationFeature
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
  
  my $constraint = $self->_internal_exclude_failed_constraint();
  
  my $from = 'structural_variation_feature svf';
  if (!$self->db->include_failed_variations()) {
    $from .= qq{ LEFT JOIN failed_structural_variation fsv 
                 ON (fsv.structural_variation_id=svf.structural_variation_id) };
  }

  my $sth = $self->prepare(qq{
    SELECT DISTINCT $cols
    FROM $from
    WHERE svf.seq_region_id = ?
      AND svf.seq_region_end > ?
      AND svf.seq_region_start < ?
      AND svf.class_attrib_id = ?
      $constraint
  });
  $sth->execute($slice->get_seq_region_id, $slice->start, $slice->end, $sv_class_id);
  
  my $result = $self->_objs_from_sth($sth);
  $sth->finish;
  
  return $result;
}



=head2 fetch_all_cnv_probe_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice 
               $slice the slice from which to obtain features
  Example    : my $svfs = $svfa->fetch_all_cnv_probe_by_Slice($slice);
  Description: Retrieves all CNV PROBE structural variation features on the given Slice.
  Returntype : reference to list Bio::EnsEMBL::StructuralVariationFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_cnv_probe_by_Slice {
  my ($self, $slice) = @_;
  my $SO_term = 'probe';
  
  return $self->fetch_all_by_Slice_SO_term($slice, $SO_term);
}


=head2 fetch_all_by_Slice_Study

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Arg [2]    : Bio::EnsEMBL:Variation::Study $study
  Arg [3]    : int $include_supporting_evidence [optional]
  Example    : my @vsfs = @{$svfa->fetch_all_by_Slice_Study($slice, $study)};
  Description: Retrieves all structural variation features in a slice that belong to a 
               given study.
               If $include_supporting_evidence is set (i.e. $include_supporting_evidence=1), structural variation features from 
               both structural variation (SV) and their supporting structural variations (SSV) will be 
               returned. By default, it only returns features from structural variations (SV). 
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_Study {

  my $self  = shift;
  my $slice = shift;
  my $study = shift;
  my $include_supporting_evidence = shift;
  
  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  if(!ref($study) || !$study->isa('Bio::EnsEMBL::Variation::Study')) {
    throw('Bio::EnsEMBL::Variation::Study arg expected');
  }
  
  # Add a constraint to only return StructuralVariationFeatures belonging to the given study
  my $constraint = $self->_internal_exclude_failed_constraint("svf.study_id = ".$study->dbID);
  
  # Include/exclude the supporting evidences
  if (!$include_supporting_evidence) {
    $constraint .= " AND " if (defined($constraint));
    $constraint .= " svf.is_evidence=0 ";
  }
  
  # Get the VariationFeatures by calling fetch_all_by_Slice_constraint
  my $svfs = $self->SUPER::fetch_all_by_Slice_constraint($slice,$constraint);

  return $svfs;
}

=head2 fetch_all_by_Slice_Source

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : Bio::EnsEMBL::Variation::Source $source
  Arg [3]    : int $include_supporting_evidence [optional]
  Example    : my @vsfs = @{$svfa->fetch_all_by_Slice_Source($slice, $source)};
  Description: Retrieves all structural variation features in a slice that belong to a 
               given source.
               If $include_supporting_evidence is set (i.e. $include_supporting_evidence=1), structural variation features from 
               both structural variation (SV) and their supporting structural variations (SSV) will be 
               returned. By default, it only returns features from structural variations (SV). 
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_Source {
  my $self  = shift;
  my $slice = shift;
  my $source = shift;
  my $include_supporting_evidence = shift;
  
  if (!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  if (!ref($source) || !$source->isa('Bio::EnsEMBL::Variation::Source')) {
    throw('Bio::EnsEMBL::Variation::Source arg expected');
  }
  
  my $constraint = $self->_internal_exclude_failed_constraint("svf.source_id = " . $source->dbID);
  
  # Include/exclude the supporting evidences
  if (!$include_supporting_evidence) {
    $constraint .= " AND " if (defined($constraint));
    $constraint .= " svf.is_evidence=0 ";
  }
  
  my $svfs = $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
  return $svfs;
}

=head2 fetch_all_by_Slice_VariationSet

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Arg [2]    : Bio::EnsEMBL:Variation::VariationSet $set
  Example    : my @vsfs = @{$svfa->fetch_all_by_Slice_VariationSet($slice, $set)};
  Description: Retrieves all structural variation features in a slice that belong to a 
               given variation set.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_VariationSet {

  my $self  = shift;
  my $slice = shift;
  my $set   = shift;
  
  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  if(!ref($set) || !$set->isa('Bio::EnsEMBL::Variation::VariationSet')) {
    throw('Bio::EnsEMBL::Variation::VariationSet arg expected');
  }
  
  # Get the bitvalue for this set and its subsets
  my $bitvalue = $set->_get_bitvalue();
  
  # Add a constraint to only return StructuralVariationFeatures having the 
  # primary keys of the supplied VariationSet or its subsets in the variation_set_id column
  my $constraint = " svf.variation_set_id & $bitvalue ";
  
  # Get the VariationFeatures by calling fetch_all_by_Slice_constraint
  my $svfs =  $self->SUPER::fetch_all_by_Slice_constraint($slice,$constraint);

  return $svfs;
}


# method used by superclass to construct SQL
sub _tables { 
  my $self = shift;
    
  my @tables = ( ['structural_variation_feature', 'svf'], [ 'source', 's'] );
  
  # If we are excluding failed_structural_variations, add that table
  push(@tables,['failed_structural_variation', 'fsv']) unless ($self->db->include_failed_variations());
    
  return @tables;
}

# Add a left join to the failed_structural_variation table
sub _left_join {
  my $self = shift;
  
  # If we are including failed structural variations, skip the left join
  return () if ($self->db->include_failed_variations());
  return (['failed_structural_variation', 'fsv.structural_variation_id=svf.structural_variation_id']);
}

sub _default_where_clause {
  my $self = shift;

  return 'svf.source_id = s.source_id';
}

sub _columns {
  return qw( svf.structural_variation_feature_id svf.seq_region_id svf.outer_start svf.seq_region_start 
             svf.inner_start svf.inner_end svf.seq_region_end svf.outer_end svf.seq_region_strand 
             svf.structural_variation_id svf.variation_name svf.source_id svf.study_id svf.class_attrib_id 
             svf.allele_string svf.somatic svf.breakpoint_order svf.length);
}

sub _objs_from_sth {
    my ($self, $sth, $mapper, $dest_slice) = @_;

    # 
    # This code is ugly because an attempt has been made to remove as many
    # function calls as possible for speed purposes.  Thus many caches and
    # a fair bit of gymnastics is used.
    #

    my $sa = $self->db()->dnadb()->get_SliceAdaptor();

    my $aa  = $self->db->get_AttributeAdaptor;

    my %slice_hash;
    my %sr_name_hash;
    my %sr_cs_hash;

    my ($structural_variation_feature_id, $seq_region_id, $outer_start, $seq_region_start, $inner_start, $inner_end, 
        $seq_region_end, $outer_end, $seq_region_strand, $structural_variation_id, $variation_name, $source_id, $study_id, $class_attrib_id, $allele_string, $is_somatic, $bp_order, $length, $last_svf_id);

    $sth->bind_columns(\$structural_variation_feature_id, \$seq_region_id, \$outer_start, \$seq_region_start, 
                       \$inner_start, \$inner_end, \$seq_region_end, \$outer_end, \$seq_region_strand, 
                       \$structural_variation_id, \$variation_name, \$source_id, \$study_id, 
                       \$class_attrib_id, \$allele_string, \$is_somatic, \$bp_order, \$length);

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
        
               {'outer_start'        => $outer_start,
                'start'              => $seq_region_start,
                'inner_start'        => $inner_start,
                'inner_end'          => $inner_end,
                'end'                => $seq_region_end,
                'outer_end'          => $outer_end,
                'strand'             => $seq_region_strand,
                'slice'              => $slice,
                'variation_name'     => $variation_name,
                'adaptor'            => $self,
                'dbID'               => $structural_variation_feature_id,
                '_source_id'         => $source_id,
                '_study_id'          => $study_id,
                '_structural_variation_id' => $structural_variation_id,
                'class_SO_term'      => $aa->attrib_value_for_id($class_attrib_id),
                'class_attrib_id'    => $class_attrib_id,
                'allele_string'      => $allele_string,
                'is_somatic'         => $is_somatic,
                'breakpoint_order'   => $bp_order,
                'length'             => $length
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


=head2 fetch_all_by_Study

  Arg [1]     : Bio::EnsEMBL::Variation::Study $study_id
  Arg [2]     : int $include_supporting_evidence [optional]
  Example     : my $study = $study_adaptor->fetch_by_name('estd1');
                foreach my $svf (@{$svf_adaptor->fetch_all_by_Study($study)}){
                   print $svf->variation_name,"\n";
                }
  Description : Retrieves all structural variation features from a specified study
                If $include_supporting_evidence is set (i.e. $include_supporting_evidence=1), structural variation features from 
                both structural variation (SV) and their supporting structural variations (SSV) will be 
                returned. By default, it only returns features from structural variations (SV).
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions  : throw if incorrect argument is passed
                warning if provided study does not have a dbID
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_by_Study {
  my $self = shift;
  my $study = shift;
  my $include_supporting_evidence = shift;

  if(!ref($study) || !$study->isa('Bio::EnsEMBL::Variation::Study')) {
    throw("Bio::EnsEMBL::Variation::Study arg expected");
  }
    
  if(!$study->dbID()) {
    warning("Study does not have dbID, cannot retrieve structural variants");
    return [];
  } 
  
  my $constraint = $self->_internal_exclude_failed_constraint('svf.study_id = '.$study->dbID);
  
  # Include/exclude the supporting evidences
  if (!$include_supporting_evidence) {
    $constraint .= " AND " if (defined($constraint));
    $constraint .= " svf.is_evidence=0 ";
  }
  
  my $result = $self->generic_fetch($constraint);

  return $result;
}


=head2 fetch_all_by_Source

  Arg [1]     : Bio::EnsEMBL::Variation::Source $source_id
  Arg [2]     : int $include_supporting_evidence [optional]
  Example     : my $source = $source_adaptor->fetch_by_name('DGVa');
                foreach my $svf (@{$svf_adaptor->fetch_all_by_Source($source)}){
                   print $svf->variation_name,"\n";
                }
  Description : Retrieves all structural variation features from a specified source
                If $include_supporting_evidence is set (i.e. $include_supporting_evidence=1), structural variation features from 
                both structural variation (SV) and their supporting structural variations (SSV) will be 
                returned. By default, it only returns features from structural variations (SV).
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions  : throw if incorrect argument is passed
                warning if provided study does not have a dbID
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_by_Source {
  my $self = shift;
  my $source = shift;
  my $include_supporting_evidence = shift;

  if(!ref($source) || !$source->isa('Bio::EnsEMBL::Variation::Source')) {
    throw("Bio::EnsEMBL::Variation::Source arg expected");
  }
    
  if(!$source->dbID()) {
    warning("Source does not have dbID, cannot retrieve structural variants");
    return [];
  } 
  
  my $constraint = $self->_internal_exclude_failed_constraint('svf.source_id = '.$source->dbID);
  
  # Include/exclude the supporting evidences
  if (!$include_supporting_evidence) {
    $constraint .= " AND " if (defined($constraint));
    $constraint .= " svf.is_evidence=0 ";
  }
  
  my $result = $self->generic_fetch($constraint);

  return $result;
}


# Exclude the constraint for failed structural variant
sub _internal_exclude_failed_constraint {
  my $self = shift;
  my $constraint = shift;
  my $no_and = shift;
  $constraint .= " AND " if (!$no_and or $constraint);
  $constraint .= $self->db->_exclude_failed_structural_variations_constraint();
  
  return $constraint;
}

# Retrieve the attribute ID in the database for the structural variation class "CNV PROBE"
sub _internal_exclude_cnv_probe {
  my $self = shift;
  my $cnv_probe_SO_term = 'probe'; # CNV_PROBE
  
  # Get the attrib_id
  my $at_adaptor = $self->db->get_AttributeAdaptor;
  my $attrib_id = $at_adaptor->attrib_id_for_type_value('SO_term','probe');
  
  return " class_attrib_id!=$attrib_id ";
}


=head2 new_fake

  Arg [1]    : string $species
  Example    :
  $vfa = Bio::EnsEMBL::Variation::StructuralVariationFeatureAdaptor->new_fake('human');
  Description: Creates a StructuralVariationFeatureAdaptor with no underlying database
         attached. Should be used only when getting consequence types for
         species with no variation database available.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariationFeatureAdaptor
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


sub store {
  my ($self, $svf) = @_;
    
  my $dbh = $self->dbc->db_handle;
    
  # look up source_id
  if(!defined($svf->{_source_id})) {
    my $sth = $dbh->prepare(q{
            SELECT source_id FROM source WHERE name = ?
    });
    $sth->execute($svf->source_name);
        
    my $source_id;
    $sth->bind_columns(\$source_id);
    $sth->fetch();
    $sth->finish();
    $svf->{_source_id} = $source_id;
  }
  throw("No source ID found for source name ", $svf->source_name) unless defined($svf->{_source_id});  
    
  # look up study_id
  if(!defined($svf->{_study_id}) && defined($svf->study)) {
    my $sth = $dbh->prepare(q{
       SELECT study_id FROM study WHERE name = ?
    });
    $sth->execute($svf->study->name);
      
    my $study_id;  
    $sth->bind_columns(\$study_id);
    $sth->fetch();
    $sth->finish();
    $svf->{_study_id} = $study_id;
  }
    
  # look up class_attrib_id
  my $class_attrib_id;
  if(defined($svf->{class_SO_term})) {
    my $sth = $dbh->prepare(q{
           SELECT attrib_id FROM attrib WHERE value = ?
    });
    $sth->execute($svf->{class_SO_term});
        
    $sth->bind_columns(\$class_attrib_id);
    $sth->fetch();
    $sth->finish();
  }
  throw("No class ID found for the class name ", $svf->{class_SO_term}) unless defined($class_attrib_id);
      
  my $sth = $dbh->prepare(q{
        INSERT INTO structural_variation_feature (
            seq_region_id,
            outer_start,
            seq_region_start,
            inner_start,
            inner_end,
            seq_region_end,
            outer_end,
            seq_region_strand,
            structural_variation_id,
            allele_string,
            variation_name,
            source_id,
            study_id,
            class_attrib_id,
            is_evidence,
            somatic,
            breakpoint_order,
            length
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    });
    
    $sth->execute(
        $svf->{slice} ? $svf->{slice}->get_seq_region_id : $svf->{seq_region_id},
        $svf->{outer_start} || undef,
        $svf->{slice} ? $svf->seq_region_start : $svf->{start},
        $svf->{inner_start} || undef,
        $svf->{inner_end} || undef,
        $svf->{slice} ? $svf->seq_region_end : $svf->{end},
        $svf->{outer_end} || undef,
        $svf->strand,
        $svf->structural_variation ? $svf->structural_variation->dbID : $svf->{_structural_variation_id},
        $svf->allele_string,
        $svf->variation_name,
        $svf->{_source_id},
        $svf->{_study_id} || undef,
        $class_attrib_id || 0,
        $svf->structural_variation ? $svf->structural_variation->is_evidence : 0,
        $svf->structural_variation ? $svf->structural_variation->is_somatic :  $svf->{is_somatic},
        $svf->{breakpoint_order} || undef,
        $svf->{length} || undef
    );
    
    $sth->finish;
    
    # get dbID
    my $dbID = $dbh->last_insert_id(undef, undef, 'structural_variation_feature', 'structural_variation_feature_id');
    $svf->{dbID}    = $dbID;
    $svf->{adaptor} = $self;
}


1;


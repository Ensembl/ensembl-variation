=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PhenotypeFeatureAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::PhenotypeFeatureAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $pfa = $reg->get_adaptor("human","variation","phenotypefeature");
  $va = $reg->get_adaptor("human","variation","variation");
  
  # Get a PhenotypeFeature by its internal identifier
  $pf = $pfa->fetch_by_dbID(45);

  # fetch all annotation for a particular variation
  $v = $va->fetch_by_name('rs56');

  foreach $pf (@{$pfa->fetch_all_by_Variation($v)}) {
    print $pf->phenotype->description(), $pf->source(), $pf->study_type(),"\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity between Variation and PhenotypeFeature objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::PhenotypeFeatureAdaptor;

use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::PhenotypeFeature qw(%TYPES);
use Bio::EnsEMBL::Variation::DBSQL::StudyAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');

# internal method
sub _fetch_all_by_object {
  my $self = shift;
  my $object = shift;
  my $type = shift;
  
  $type ||= (split '::', ref($object))[-1];
  throw("$type is not a valid object type, valid types are: ".(join ", ", sort %TYPES)) unless defined $type and defined($TYPES{$type});
  
  return $self->generic_fetch("pf.type = '".$type."' AND pf.object_id = '".$object->stable_id."'");
}


=head2 fetch_all_by_object_id

  Arg [1]    : string $id
  Arg [2]    : (optional) string $type
  Example    : my @pfs = @{$pfa->fetch_all_by_object_id('rs1333049')};
               my @pfs = @{$pfa->fetch_all_by_object_id('Activ1', 'QTL')};
  Description: Retrieves all phenotype features of the given object ID. Since a
               PhenotypeFeature may be associated with several object types, the
               object type may also be supplied to unambiguously identify the
               associated object.
  Returntype : reference to list Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_object_id {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  
  my $constraint = qq{ pf.object_id = $id };
  
  if(defined($type)) {
    throw("$type is not a valid object type, valid types are: ".(join ", ", sort %TYPES)) unless defined($TYPES{$type});
    $constraint .= qq{ AND pf.type = '$type' } if defined($type);
  }
  
  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Slice_type

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : string $type
  Example    : my @pfs = @{$pfa->fetch_all_by_Slice_type($slice, 'QTL')};
  Description: Retrieves all phenotype features of the given type on a slice.
  Returntype : reference to list Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_type {
  my $self = shift;
  my $slice = shift;
  my $type = shift;
  
  throw("No valid object type given, valid types are: ".(join ", ", sort %TYPES)) unless defined $type and defined($TYPES{$type});
  
  return $self->fetch_all_by_Slice_constraint($slice, "pf.type = '$type'");
}


=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation::Variation $var
  Example    : my @pfs = @{$pfa->fetch_all_by_Variation($var)};
  Description: Retrieves all PhenotypeFeatures for a given variation.
  Returntype : reference to list Bio::EnsEMBL::Variation::PhenotypeFeature
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
  
  return $self->_fetch_all_by_object($var, 'Variation');
}


=head2 fetch_all_by_Variation_list

  Arg [1]    : reference to a list of Bio::EnsEMBL::Variation::Variation objects
  Example    : my @pfs = @{$pfa->fetch_all_by_Variation_list($vars)};
  Description: Retrieves all PhenotypeFeatures for a given list of variations
  Returntype : reference to a list of Bio::EnsEMBL::Variation::PhenotypeFeature objects
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Variation_list {
  my $self = shift;
  my $vars  = shift;

  if(!ref($vars) || !$vars->[0]->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation arg expected');
  }

  if(!defined($vars->[0]->name())) {
    throw("Variation arg must have defined name");
  }
  
  my $in_str = join ',', map {"'".$_->name()."'"} @$vars;

  return $self->generic_fetch("pf.type = 'Variation' AND pf.object_id in (".$in_str.")");
}


=head2 fetch_all_by_StructuralVariation

  Arg [1]    : Bio::EnsEMBL::Variation::StructuralVariation $sv
  Example    : my @pfs = @{$pfa->fetch_all_by_StructuralVariation($sv)};
  Description: Retrieves all PhenotypeFeatures for a given structural variation.
  Returntype : reference to list Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_StructuralVariation {
  my $self = shift;
  my $sv   = shift;

  if(!ref($sv) || !$sv->isa('Bio::EnsEMBL::Variation::BaseStructuralVariation')) {
    throw('Bio::EnsEMBL::Variation::BaseStructuralVariation arg expected');
  }
  
  # don't provide type here as it could be a SupportingStructuralVariation
  return $self->_fetch_all_by_object($sv);
}


=head2 fetch_all_by_Gene

  Arg [1]    : Bio::EnsEMBL::Gene $gene
  Example    : my @pfs = @{$pfa->fetch_all_by_Gene($gene)};
  Description: Retrieves all PhenotypeFeatures for a given Gene.
  Returntype : reference to list Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Gene {
  my $self = shift;
  my $gene = shift;

  if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw('Bio::EnsEMBL::Gene arg expected');
  }
  
  return $self->_fetch_all_by_object($gene, 'Gene');
}


=head2 fetch_all_by_VariationFeature_list

  Arg [1]    : reference to a list of Bio::EnsEMBL::Variation::VariationFeature objects
  Example    : my @pfs = @{$pfa->fetch_all_by_VariationFeature_list($vfs)};
  Description: Retrieves all PhenotypeFeatures for a given list of variation features
  Returntype : reference to a list Bio::EnsEMBL::Variation::PhenotypeFeature objects
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_VariationFeature_list {
  my $self = shift;
  my $vfs  = shift;

  if(!ref($vfs) || !$vfs->[0]->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Listref of Bio::EnsEMBL::Variation::VariationFeature arg expected');
  }

  if(!defined($vfs->[0]->variation_name())) {
    throw("VariationFeatures in list must have defined names");
  }
  
  my $in_str = join ',', map {"'".$_->variation_name."'"} @$vfs;

  return $self->generic_fetch("pf.type = 'Variation' AND pf.object_id in (".$in_str.")");
}


=head2 fetch_all_by_Study

  Arg [1]    : Bio::EnsEMBL:Variation::Study $study
  Example    : my @studies = @{$studya->fetch_all_by_Study($study)};
  Description: Retrieves all PhenotypeFeatures for a given study.
  Returntype : reference to list Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Study {
  my $self   = shift;
  my $study  = shift;

  if(!ref($study) || !$study->isa('Bio::EnsEMBL::Variation::Study')) {
    throw('Bio::EnsEMBL::Variation::Study arg expected');
  }

  if(!defined($study->dbID())) {
    throw("Study arg must have defined dbID");
  }
  return $self->generic_fetch("pf.study_id = ".$study->dbID());
}


=head2 fetch_all_by_phenotype_name_source_name

  Arg [1]    : string $phenotype_name
  Arg [2]    : string $source_name (optional)
  Example    : $pf = $pf_adaptor->fetch_all_by_phenotype_source_name('BD','EGA');
  Description: Retrieves a PhenotypeFeature object via its phenotype/source name
  Returntype : list of ref of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw if phenotype name argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_phenotype_name_source_name {

  my $self = shift;
  my $phenotype_name  = shift;
  my $source_name = shift;

  throw('phenotype_name argument expected') if(!defined($phenotype_name));

  my $extra_sql = " p.name = $phenotype_name ";
  if (defined $source_name ) {
    $extra_sql .= qq( AND s.name = '$source_name' );
  }
  
  # Add the constraint for failed variations
  #$extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch("$extra_sql");
  
}

=head2 fetch_all_by_phenotype_description_source_name

  Arg [1]    : string $phenotype_description
  Arg [2]    : string $source_name (optional)
  Example    : $pf = $pf_adaptor->fetch_all_by_phenotype_description_source_name('diabetes','EGA');
  Description: Retrieves a PhenotypeFeature object via its phenotype description/source name
  Returntype : list of ref of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw if phenotype name argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_phenotype_description_source_name {

  my $self = shift;
  my $phenotype_description  = shift;
  my $source_name = shift;

  throw('phenotype_description argument expected') if(!defined($phenotype_description));

  my $extra_sql = qq( p.description like '%$phenotype_description%' );
  if (defined $source_name ) {
    $extra_sql .= qq( AND s.name = '$source_name' );
  }
  
  # Add the constraint for failed variations
  #$extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch("$extra_sql");
  
}

=head2 fetch_all_by_phenotype_id_source_name

  Arg [1]    : integer $phenotype_id
  Arg [2]    : string $source_name (optional)
  Example    : $pf = $pf_adaptor->fetch_all_by_phenotype_id_source_name(999,'EGA');
  Description: Retrieves a PhenotypeFeature object via its phenotype id/source name
  Returntype : list of ref of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw if phenotype id argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_phenotype_id_source_name {

  my $self = shift;
  my $phenotype_id  = shift;
  my $source_name = shift;

  throw('phenotype_id argument expected') if(!defined($phenotype_id));

  my $extra_sql = sprintf('p.phenotype_id = %s', $self->dbc->db_handle->quote( $phenotype_id, SQL_INTEGER ) );

  if (defined $source_name ) {
    $extra_sql .= sprintf(" AND s.name = '%s'", $self->dbc->db_handle->quote( $source_name, SQL_VARCHAR ) );
  }
  
  # Add the constraint for failed variations
  #$extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch("$extra_sql");
  
}


=head2 fetch_all_by_associated_gene

  Arg [1]    : string $gene_name
  Example    : $pf = $pf_adaptor->fetch_all_by_associated_gene('CAV3');
  Description: Retrieves the PhenotypeFeature objects via which are associated with the gene.
  Returntype : list of ref of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw if the gene_name argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_associated_gene {

  my $self = shift;
  my $gene_name  = shift;

  throw('gene_name argument expected') if(!defined($gene_name));
  
  $self->_include_attrib(1);

	my $extra_sql = " at.code = 'associated_gene' and pfa.value REGEXP '^(.+,)?[. .]*$gene_name(,.+)?\$'";
  
  # Add the constraint for failed variations
  #$extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();

  my $result = $self->generic_fetch($extra_sql);
  
  $self->_include_attrib(0);
  
  return $result;
}

=head2 count_all_by_associated_gene

  Description: Retrieves count of variation_annotation objects associated with a
               given gene
  Returntype : integer
  Exceptions : none
  Caller     : general

=cut

sub count_all_by_associated_gene {

  my $self = shift;
  my $gene_name  = shift;

  throw('gene_name argument expected') if(!defined($gene_name));
  
  $self->_include_attrib(1);

  my $extra_sql = " at.code = 'associated_gene' and pfa.value REGEXP '^(.+,)?[. .]*$gene_name(,.+)?\$'";
  
  # Add the constraint for failed variations
  #$extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();

  my $result = $self->generic_count($extra_sql);
  
  $self->_include_attrib(0);
  
  return $result;
}

sub count_all_by_Phenotype {
  my $self = shift;
  return $self->count_all_by_phenotype_id($_[0]->dbID);
}

sub count_all_by_Gene {
  my $self = shift;
  my $gene = shift;
  
  return $self->generic_count("pf.object_id = '".$gene->stable_id."' AND pf.type = 'Gene'");
}

sub count_all_by_phenotype_id {
  my $self = shift;
  my $id = shift;
  return $self->generic_count("pf.phenotype_id = $id");
}


=head2 fetch_all

  Description: Retrieves all available PhenotypeFeature objects
  Returntype : list of ref of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : none
  Caller     : general

=cut

sub fetch_all {

  my $self = shift;
  
  # Add the constraint for failed variations
  #my $extra_sql = $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch();#"$extra_sql");
  
}

# stub method used by web
sub _check_gene_by_HGNC {
  my $self = shift;
  my $hgnc = shift;
  
  my $sth = $self->dbc->prepare(qq{
    SELECT count(*)
    FROM phenotype_feature_attrib pfa, attrib_type at
    WHERE pfa.attrib_type_id = at.attrib_type_id
    AND at.code = 'associated_gene'
    AND pfa.value = ?
  });
  
  $sth->bind_param(1, $hgnc, SQL_VARCHAR);
  $sth->execute();
  
  my $count;
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish();
  
  return $count;
}

# fetches attributes
sub _fetch_attribs_by_dbID {
  my $self = shift;
  my $id = shift;
  
  throw("Cannot fetch attributes without dbID") unless defined($id);
  
  my $attribs = {};
  
  my $sth = $self->dbc->prepare(qq{
    SELECT at.code, a.value
    FROM phenotype_feature_attrib a, attrib_type at
    WHERE a.attrib_type_id = at.attrib_type_id
    AND a.phenotype_feature_id = ?
  });
  
  $sth->bind_param(1,$id,SQL_INTEGER);
  $sth->execute();
  
  my ($key, $value);
  
  $sth->bind_columns(\$key, \$value);
  $attribs->{$key} = $value while $sth->fetch;
  $sth->finish;
  
  return $attribs;
}

sub _include_attrib {
  my $self = shift;
  
  # default on first call
  $self->{_include_attrib} = 0 unless defined($self->{_include_attrib});
  
  # update if given
  $self->{_include_attrib} = shift if @_;
  
  # return
  return $self->{_include_attrib};
}

# method used by superclass to construct SQL
sub _tables {
  my $self = shift;
  
  my @tables = (
    [ 'phenotype_feature', 'pf' ],
    [ 'phenotype', 'p' ],
    [ 'source', 's' ]
  );
  
  # include attrib tables?
  push @tables, (
    [ 'phenotype_feature_attrib', 'pfa', ],
    [ 'attrib_type', 'at' ]
  ) if $self->_include_attrib;
  
  return @tables; 
}

# Add left join to the source table
sub _left_join {
  my $self = shift;
  
  my @lj = (
    [ 'source', 'pf.source_id = s.source_id' ],
    [ 'phenotype', 'pf.phenotype_id = p.phenotype_id' ]
  );
  
  push @lj, (
    [ 'phenotype_feature_attrib', 'pf.phenotype_feature_id = pfa.phenotype_feature_id' ],
    [ 'attrib_type', 'pfa.attrib_type_id = at.attrib_type_id' ]
  ) if $self->_include_attrib;
  
  return @lj;
}

#sub _default_where_clause {
#  my $self = shift;
#
#  return 'pf.phenotype_id = p.phenotype_id';
#}

sub _columns {
  return qw(
    pf.phenotype_feature_id pf.object_id pf.type pf.is_significant
    pf.seq_region_id pf.seq_region_start pf.seq_region_end pf.seq_region_strand
    pf.phenotype_id s.name pf.study_id
  );
}



sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;
  
  # get required adaptors
  my $sa = $self->db()->dnadb()->get_SliceAdaptor();
  my $aa = $self->db->get_AttributeAdaptor;

  my (
    @features, %slice_hash, %sr_name_hash, %sr_cs_hash,
    $asm_cs, $cmp_cs, $asm_cs_vers, $asm_cs_name, $cmp_cs_vers, $cmp_cs_name,
    $dest_slice_start, $dest_slice_end, $dest_slice_strand, $dest_slice_length
  );
  
  if($mapper) {
    $asm_cs = $mapper->assembled_CoordSystem();
    $cmp_cs = $mapper->component_CoordSystem();
    $asm_cs_name = $asm_cs->name();
    $asm_cs_vers = $asm_cs->version();
    $cmp_cs_name = $cmp_cs->name();
    $cmp_cs_vers = $cmp_cs->version();
  }
  
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
  }
  
  my (
    $phenotype_feature_id, $object_id, $object_type, $is_significant,
    $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand,
    $phenotype_id, $source_name, $study_id
  );
  
  $sth->bind_columns(
    \$phenotype_feature_id, \$object_id, \$object_type, \$is_significant,
    \$seq_region_id, \$seq_region_start, \$seq_region_end, \$seq_region_strand,
    \$phenotype_id, \$source_name, \$study_id
  );
  
  my $sta = $self->db()->get_StudyAdaptor();
  
  FEATURE: while($sth->fetch()) {
    
    # remap
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
      }
      else {
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
        }
        else {
          my $tmp_seq_region_start = $seq_region_start;
          $seq_region_start = $dest_slice_end - $seq_region_end + 1;
          $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
          #$seq_region_strand *= -1;
        }
        
        #throw away features off the end of the requested slice
        if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
          next FEATURE;
        }
      }
      $slice = $dest_slice;
    }
    
    push @features, $self->_create_feature_fast(
      'Bio::EnsEMBL::Variation::PhenotypeFeature', {
        'dbID'           => $phenotype_feature_id,
        'start'          => $seq_region_start,
        'end'            => $seq_region_end,
        'strand'         => $seq_region_strand,
        'slice'          => $slice,
        '_object_id'     => $object_id,
        'type'           => $object_type,
        '_phenotype_id'  => $phenotype_id,
        'adaptor'        => $self,
        '_study_id'      => $study_id,
        'source'         => $source_name,
        'is_significant' => $is_significant,
      }
    );
  }

  return \@features;

}

1;

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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PhenotypeFeatureAdaptor
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
    print $pf->phenotype->description(), $pf->source(),"\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity between Variation and PhenotypeFeature objects.

By default, only the phenotype features with associations labelled as "significant" will be fetched.
To fetch all the phenotype features, set the DBAdaptor variable "include_non_significant_phenotypes", 
e.g. $pfa->db->include_non_significant_phenotype_associations(1);
See Bio::EnsEMBL::Variation::DBSQL::DBAdaptor


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
=head2 _is_significant_constraint

  Arg [1]    : string $constraint
  Example    : $self->_is_significant_constraint($constraint)
  Description: Internal method to add a constraint on the column "is_significant".
               By default, the constraint is on the significant value (is_significant=1).
               To remove this constraint, set the DBAdaptor variable "include_non_significant_phenotypes", using the
               method "include_non_significant_phenotype_associations":  
               e.g. $pfa->db->include_non_significant_phenotype_associations(1);
               See the module Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::include_non_significant_phenotype_associations
  Returntype : string
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _is_significant_constraint {
  my $self = shift;
  my $constraint = shift;
    
  # If we should include non significant objects, no extra condition is needed
  return $constraint if ($self->db->include_non_significant_phenotype_associations());

  # Otherwise, add a constraint on the phenotype_feature table
  my $ns_constraint = qq{ pf.is_significant=1 };
  $constraint  .= (defined($constraint)) ? " AND$ns_constraint" : $ns_constraint;
    
  return $constraint;
}


# internal method
sub _fetch_all_by_object {
  my $self   = shift;
  my $object = shift;
  my $type   = shift;
  
  $type ||= (split '::', ref($object))[-1];
  throw("$type is not a valid object type, valid types are: ".(join ", ", sort keys %TYPES)) unless defined $type and defined($TYPES{$type});
  
  my $constraint = "pf.type = '".$type."' AND pf.object_id = '".$object->stable_id."'";
  
  # Add the constraint for significant data
  $constraint = $self->_is_significant_constraint($constraint);
     
  return $self->generic_fetch($constraint);
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
  
  my $constraint = qq{ pf.object_id = '$id' };
  
  if(defined($type)) {
    throw("$type is not a valid object type, valid types are: ".(join ", ", sort keys %TYPES)) unless defined($TYPES{$type});
    $constraint .= qq{ AND pf.type = '$type' } if defined($type);
  }
  
  # Add the constraint for significant data
  $constraint = $self->_is_significant_constraint($constraint);
  
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
  my $self  = shift;
  my $slice = shift;
  my $type  = shift;
  
  throw("$type is not a valid object type, valid types are: ".(join ", ", sort keys %TYPES)) unless defined $type and defined($TYPES{$type});
  
  my $constraint = qq{pf.type = '$type'};
  
  # Add the constraint for significant data
  $constraint = $self->_is_significant_constraint($constraint);
  
  return $self->fetch_all_by_Slice_constraint($slice, $constraint);
}


=head2 fetch_all_by_Slice_Study

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Arg [2]    : Bio::EnsEMBL:Variation::Study $study
  Example    : my @pfs = @{$pfa->fetch_all_by_Slice_Study($slice, $study)};
  Description: Retrieves all phenotype features in a slice that belong to a 
               given study.
  Returntype : reference to list Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_Study {

  my $self  = shift;
  my $slice = shift;
  my $study = shift;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Bio::EnsEMBL::Slice arg expected');
  }
  if(!ref($study) || !$study->isa('Bio::EnsEMBL::Variation::Study')) {
    throw('Bio::EnsEMBL::Variation::Study arg expected');
  }
  if(!defined($study->dbID())) {
    throw("Study arg must have defined dbID");
  }

  # Add a constraint to only return PhenotypeFeatures belonging to the given study, within the given slice
  my $constraint = "pf.study_id = ".$study->dbID;

  # Get the VariationFeatures by calling fetch_all_by_Slice_constraint
  my $pfs = $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);

  return $pfs;
}

=head2 fetch_all_by_Slice_with_ontology_accession

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : (optional) string $type
  Example    : my @pfs = @{$pfa->fetch_all_by_Slice_with_ontology_accession($slice, 'Variation')};
  Description: Retrieves all phenotype features of the given type (e.g. 'Variation','StructuralVariation','Gene','QTL')
               on a slice with its phenotype ontology accession data.
  Returntype : reference to list Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_with_ontology_accession {
  my $self  = shift;
  my $slice = shift;
  my $type  = shift;

  my $constraint;
  if ($type) {
    throw("$type is not a valid object type, valid types are: ".(join ", ", sort keys %TYPES)) unless defined $type and defined($TYPES{$type});
    $constraint .= qq{pf.type = '$type'};
  }

  $self->_include_ontology(1);

  # Add the constraint for significant data
  $constraint = $self->_is_significant_constraint($constraint);

  my $result = $self->fetch_all_by_Slice_constraint($slice, $constraint);

  ## reset flags
  $self->_include_ontology(0);

  return $result;
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
  my $vars = shift;

  if(!ref($vars) || !$vars->[0]->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation arg expected');
  }

  if(!defined($vars->[0]->name())) {
    throw("Variation arg must have defined name");
  }
  
  my $in_str = join ',', map {"'".$_->name()."'"} @$vars;

  my $constraint = qq{pf.type = 'Variation' AND pf.object_id in ($in_str)};
  
  # Add the constraint for significant data
  $constraint = $self->_is_significant_constraint($constraint);

  return $self->generic_fetch($constraint);
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

  my $constraint = qq{pf.type = 'Variation' AND pf.object_id in ($in_str)};
  
  # Add the constraint for significant data
  $constraint = $self->_is_significant_constraint($constraint);
     
  return $self->generic_fetch($constraint);
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
  
  my $constraint = "pf.study_id = ".$study->dbID();
  
  # Add the constraint for significant data
  $constraint = $self->_is_significant_constraint($constraint);
  
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

  my $extra_sql = " p.name = '$phenotype_name' ";
  if (defined $source_name ) {
    $extra_sql .= qq( AND s.name = '$source_name' );
  }
  
  # Add the constraint for significant data
  $extra_sql = $self->_is_significant_constraint($extra_sql);
  
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

  my $extra_sql = qq( p.description like "%$phenotype_description%" );
  if (defined $source_name ) {
    $extra_sql .= qq( AND s.name = '$source_name' );
  }
  
  # Add the constraint for significant data
  $extra_sql = $self->_is_significant_constraint($extra_sql);
  
  # Add the constraint for failed variations
  #$extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch("$extra_sql");
  
}

=head2 fetch_all_by_Phenotype

  Arg [1]    : Bio::EnsEMBL::Variation::Phenotype
  Example    : $pf = $pf_adaptor->fetch_all_by_Phenotype($phenotype_object);
  Description: Retrieves a PhenotypeFeature object via its phenotype id/source name
  Returntype : list of ref of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw if phenotype argument is not supplied
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Phenotype {

  my $self = shift;
  my $phenotype  = shift;

  if(!ref($phenotype) || !$phenotype->isa('Bio::EnsEMBL::Variation::Phenotype')) {
    throw('Bio::EnsEMBL::Variation::Phenotype arg expected');
  }

  return $self->fetch_all_by_phenotype_id_source_name($phenotype->dbID());
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
    $extra_sql .= sprintf(" AND s.name = %s", $self->dbc->db_handle->quote( $source_name, SQL_VARCHAR ) );
  }
  
  # Add the constraint for significant data
  $extra_sql = $self->_is_significant_constraint($extra_sql);
  
  # Add the constraint for failed variations
  #$extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch("$extra_sql");
  
}

=head2 fetch_all_by_phenotype_id_feature_type

  Arg [1]    : integer $phenotype_id
  Arg [2]    : string  feature type
  Example    : $pf = $pf_adaptor->fetch_all_by_phenotype_id_feature_type(999,'Gene');
  Description: Retrieves a PhenotypeFeature object via its phenotype id and feature type
  Returntype : list of ref of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw if phenotype id argument is not defined or type not supported
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_phenotype_id_feature_type {

  my $self = shift;
  my $phenotype_id  = shift;
  my $feature_type  = shift;

  throw('phenotype_id argument expected') if(!defined($phenotype_id));

  throw('feature_type not recognised') unless $TYPES{$feature_type};

  my $extra_sql = sprintf('p.phenotype_id = %s ', $self->dbc->db_handle->quote( $phenotype_id, SQL_INTEGER ) );
  $extra_sql .= sprintf(' AND pf.type = %s ', $self->dbc->db_handle->quote( $feature_type, SQL_VARCHAR ) );


  # Add the constraint for significant data
  $extra_sql = $self->_is_significant_constraint($extra_sql);

  return $self->generic_fetch("$extra_sql");

}

=head2 fetch_all_by_phenotype_accession_source

  Arg [1]    : string phenotype ontology_accession
  Arg [2]    : string source name (optional)
  Example    : $pf = $pf_adaptor->fetch_all_by_phenotype_accession_source('EFO:0004330','ClinVar');
  Description: Retrieves a PhenotypeFeature object via an ontology accession and optional source
  Returntype : list of ref of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw if phenotype ontology accession argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_phenotype_accession_source {

  my $self      = shift;
  my $accession = shift;
  my $source    = shift;

  throw('phenotype ontology accession argument expected') if(!defined($accession));
 
  $self->_include_ontology(1);
  $self->_include_attrib(1);

  my $extra_sql = ' poa.accession = "' . $accession . '" ';
  $extra_sql   .= ' and s.name = "' . $source . '" ' if $source;


  my $result = $self->generic_fetch($extra_sql);

  ## reset flags  
  $self->_include_ontology(0);
  $self->_include_attrib(0);

  return $result;
}

=head2 fetch_all_by_phenotype_accession_type_source

  Arg [1]    : string phenotype ontology_accession
  Arg [2]    : mapping type - default 'is', option 'involves'
  Arg [3]    : string source name (optional)
  Example    : $pf = $pf_adaptor->fetch_all_by_phenotype_accession_source('EFO:0004330','is', ClinVar');
  Description: Retrieves a PhenotypeFeature object via an ontology accession mapping type and optional source
  Returntype : list of ref of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw if phenotype ontology accession argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_phenotype_accession_type_source {

  my $self      = shift;
  my $accession = shift;
  my $map_type  = shift;
  my $source    = shift;

  throw('phenotype ontology accession argument expected') if(!defined($accession));
  $map_type ||= 'is';

  $self->_include_ontology(1);
  $self->_include_attrib(1);

  my $extra_sql = ' poa.accession = "' . $accession . '" ';
  $extra_sql   .= ' and poa.mapping_type = "' . $map_type . '" ';
  $extra_sql   .= ' and s.name = "' . $source . '" ' if $source;


  my $result = $self->generic_fetch($extra_sql);

  ## reset flags  
  $self->_include_ontology(0);
  $self->_include_attrib(0);

  return $result;
}

    
=head2 fetch_all_by_associated_gene_phenotype_description

  Arg [1]    : string $gene_name
  Arg [2]    : string $phenotype
  Example    : $pf = $pf_adaptor->fetch_all_by_associated_gene_phenotype_description('HFE','Blood pressure');
  Description: Retrieves the PhenotypeFeature objects via which are associated with the gene, for a given phenotype.
  Returntype : list of ref of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : throw if the gene_name and the phenotype arguments are not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_associated_gene_phenotype_description {
  my $self = shift;
  my $gene_name = shift;
  my $phenotype = shift;
  
  throw('gene_name argument expected') if(!defined($gene_name));
  throw('phenotype argument expected') if(!defined($phenotype));
  
  my $constraint = qq( p.description="$phenotype" );
  
  return $self->fetch_all_by_associated_gene($gene_name,$constraint);
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
  my $constraint = shift;

  throw('gene_name argument expected') if(!defined($gene_name));

  $self->_include_attrib(1);

  my $extra_sql  = " at.code = 'associated_gene' and (pfa.value = '$gene_name' OR pfa.value like '%,$gene_name' OR pfa.value like '$gene_name,%' OR pfa.value like '%,$gene_name,%')";
     $extra_sql .= " and $constraint" if ($constraint);
  # Add the constraint for significant data
  $extra_sql = $self->_is_significant_constraint($extra_sql);
  
  # Add the constraint for failed variations
  #$extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();

  my $result = $self->generic_fetch($extra_sql);
  
  $self->_include_attrib(0);
  $result = [grep {$_->type ne 'Gene'} @$result];  
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

  my $extra_sql  = " at.code = 'associated_gene' and pfa.value REGEXP '^(.+,)?[. .]*$gene_name(,.+)?\$'";
     $extra_sql .= " and pf.type!='Gene'";
  # Add the constraint for significant data
  $extra_sql = $self->_is_significant_constraint($extra_sql);
  
  # Add the constraint for failed variations
  #$extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();

  my $result = $self->generic_count($extra_sql);
  
  $self->_include_attrib(0);

  return $result;
}



=head2 count_all_by_Phenotype

  Arg [1]    : Bio::EnsEMBL:Variation::Phenotype object
  Example    : $count = $pf_adaptor->count_all_by_Phenotype($phe_object);
  Description: Retrieves count of the phenotype_feature objects associated with a
               given phenotype
  Returntype : integer
  Exceptions : none
  Caller     : general

=cut

sub count_all_by_Phenotype {
  my $self = shift;
  return $self->count_all_by_phenotype_id($_[0]->dbID);
}

=head2 count_all_by_Gene

  Arg [1]    : Bio::EnsEMBL:Gene object
  Example    : $count = $pf_adaptor->count_all_by_Gene($gene_object);
  Description: Retrieves count of the phenotype_feature objects associated with a
               given gene
  Returntype : integer
  Exceptions : none
  Caller     : general

=cut
sub count_all_by_Gene {
  my $self = shift;
  my $gene = shift;
  
  my $constraint = "pf.object_id = '".$gene->stable_id."' AND pf.type = 'Gene'";
  
  # Add the constraint for significant data
  $constraint = $self->_is_significant_constraint($constraint);
  
  return $self->generic_count($constraint);
}


=head2 count_all_by_phenotype_id

  Description: Retrieves phenotype_feature counts for a given phenotype
  Returntype : the phenotype_feature counts
  Exceptions : none
  Caller     : web

=cut
sub count_all_by_phenotype_id {
  my $self = shift;
  my $id = shift;
  
  my $constraint = qq{pf.phenotype_id = $id};
  
  # Add the constraint for significant data
  $constraint = $self->_is_significant_constraint($constraint);
  
  return $self->generic_count($constraint);
}

=head2 count_all_type_by_phenotype_id

  Description: Retrieves phenotype_feature counts by type
               (e.g. Variation, StructuralVariation, Gene, QTL)
  Returntype : a hash ref type => phenotype_feature counts
  Exceptions : none
  Caller     : web

=cut
sub count_all_type_by_phenotype_id {
  my $self = shift;
  my $id = shift;

  my %count_by_type;

  my $constraint = qq{pf.phenotype_id = ?};

  # Add the constraint for significant data
  $constraint = $self->_is_significant_constraint($constraint);

  my $sth = $self->dbc->prepare(qq{
    SELECT pf.type, count(*)
    FROM phenotype_feature pf
    WHERE $constraint
    GROUP BY pf.type
  });

  $sth->execute($id);
  my $counts = $sth->fetchall_arrayref();

  foreach my $c (@{$counts}){
    $count_by_type{$c->[0]} = $c->[1];
  }

  return \%count_by_type;
}

=head2 count_all_with_source_by_Phenotype

  Description: Retrieves phenotype_feature counts by source name
  Returntype : a hash ref source name => phenotype_feature counts
  Exceptions : none
  Caller     : web

=cut
sub count_all_with_source_by_Phenotype {
  my $self     = shift;
  my $phenotype = shift;

  my %count_by_source;

  my $sth = $self->dbc->prepare(qq{
    SELECT s.name, count(*)
    FROM phenotype_feature pf, source s
    where pf.phenotype_id = ?
    and pf.source_id = s.source_id 
    group by s.name
   });

  $sth->execute($phenotype->dbID);
  my $counts = $sth->fetchall_arrayref();

  foreach my $c (@{$counts}){
    $count_by_source{$c->[0]} = $c->[1];
  }

  return \%count_by_source;

}

=head2 fetch_all

  Description: Retrieves all available PhenotypeFeature objects
  Returntype : list of ref of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : none
  Caller     : general

=cut

sub fetch_all {

  my $self = shift;
  
  # Add the constraint for significant data
  my $extra_sql = $self->_is_significant_constraint();
  
  # Add the constraint for failed variations
  #my $extra_sql = $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch("$extra_sql");
  
}

# stub method used by web
sub _check_gene_by_HGNC {
  my $self = shift;
  my $hgnc = shift;

  my $extra_sql = $self->_is_significant_constraint();
  
  my $sth = $self->dbc->prepare(qq{
    SELECT count(*)
    FROM phenotype_feature pf, phenotype_feature_attrib pfa, attrib_type at, source s
    WHERE pf.phenotype_feature_id = pfa.phenotype_feature_id
    AND pfa.attrib_type_id = at.attrib_type_id
    AND pf.source_id = s.source_id
    AND s.name != 'COSMIC'
    AND at.code = 'associated_gene'
    AND $extra_sql
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
  while ($sth->fetch){

    if($key =~/submitter_id/){
      ## expand the attrib ids stored in phenotype_feature_attrib to the full submitter names
      foreach my $id(split/\,/,$value){
        push @{ $attribs->{submitter_names}}, $self->_get_submitter_name($id);
      }
    }
    else{
      $attribs->{$key} = $value;
    }
  }
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

sub _include_ontology {
  my $self = shift;
  
  # default on first call
  $self->{_include_ontology} = 0 unless defined($self->{_include_ontology});
  
  # update if given
  $self->{_include_ontology} = shift if @_;
  
  # return
  return $self->{_include_ontology};
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
    [ 'phenotype_feature_attrib', 'pfa' ],
    [ 'attrib_type', 'at' ]
  ) if $self->_include_attrib;

  # include ontology tables for search?
  push @tables, (
    [ 'phenotype_ontology_accession', 'poa' ]
  ) if $self->_include_ontology;
 
  return @tables; 
}

# Add left join to the source table
sub _left_join {
  my $self = shift;
  
  my @lj = ();
  
  push @lj, (
    [ 'phenotype_feature_attrib', 'pf.phenotype_feature_id = pfa.phenotype_feature_id' ],
    [ 'attrib_type', 'pfa.attrib_type_id = at.attrib_type_id' ]
  ) if $self->_include_attrib;

  push @lj, (
    [ 'phenotype_ontology_accession', 'pf.phenotype_id = poa.phenotype_id' ]
  ) if $self->_include_ontology;
  
  return @lj;
}

sub _default_where_clause {
  my $self = shift;

  return 'pf.phenotype_id = p.phenotype_id and pf.source_id = s.source_id';
}

sub _columns {
  my $self = shift;

  my @cols = qw(
    pf.phenotype_feature_id pf.object_id pf.type pf.is_significant
    pf.seq_region_id pf.seq_region_start pf.seq_region_end pf.seq_region_strand
    pf.phenotype_id pf.source_id s.name pf.study_id p.description
  ) ;
  push @cols, ('pfa.value','at.code') if $self->_include_attrib;
  push @cols, 'poa.accession'         if $self->_include_ontology;

  return @cols;
}



sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;
   
  my %row;
  # Create the row hash using column names as keys
  $sth->bind_columns( \( @row{ @{$sth->{NAME_lc} } } ));

  while ($sth->fetch) {

      # we don't actually store the returned object because
      # the _obj_from_row method stores them in a temporary
      # hash _temp_objs in $self 
      $self->_obj_from_row(\%row, $mapper, $dest_slice);
  }

  # Get the created objects from the temporary hash
  my @objs = values %{ $self->{_temp_objs} };
  delete $self->{_temp_objs};
 
  # Return the created objects 
  return \@objs;
}


sub _obj_from_row {

  my ($self, $row, $mapper, $dest_slice) = @_;


  my $obj = $self->{_temp_objs}{$row->{phenotype_feature_id}}; 
    
  unless (defined($obj)) {
 
    # get required adaptors
    my $sa = $self->db()->dnadb()->get_SliceAdaptor();

    my (
      @features, %slice_hash, %sr_name_hash, %sr_cs_hash,
      $asm_cs, $cmp_cs, $asm_cs_vers, $asm_cs_name, $cmp_cs_vers, $cmp_cs_name,
      $dest_slice_start, $dest_slice_end, $dest_slice_strand, $dest_slice_length,
      $sr_name
    );

    my $seq_region_start   = $row->{seq_region_start};
    my $seq_region_end     = $row->{seq_region_end};
    my $seq_region_strand  = $row->{seq_region_strand};

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
  
    my $sta = $self->db()->get_StudyAdaptor();

    # remap
    my $slice = $slice_hash{"ID:".$row->{seq_region_id}};
    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($row->{seq_region_id}, undef, undef, undef, 1);
      next unless $slice;
      $slice_hash{"ID:".$row->{seq_region_id}} = $slice;
      $sr_name_hash{$row->{seq_region_id}} = $slice->seq_region_name();
      $sr_cs_hash{$row->{seq_region_id}} = $slice->coord_system();
    }
    
    # remap the feature coordinates to another coord system
    # if a mapper was provided
    if($mapper) {
      my $sr_name = $sr_name_hash{$row->{seq_region_id}};
      my $sr_cs   = $sr_cs_hash{$row->{seq_region_id}};
      
      ($sr_name,$seq_region_start,$seq_region_end,$seq_region_strand) =
      $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
      $seq_region_strand, $sr_cs);
      
      #skip features that map to gaps or coord system boundaries
      return if(!defined($sr_name));
      
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
          return;
        }
      }
      $slice = $dest_slice;
    }

    $obj = $self->_create_feature_fast(
      'Bio::EnsEMBL::Variation::PhenotypeFeature', {
        'dbID'           => $row->{phenotype_feature_id},
        'start'          => $seq_region_start,
        'end'            => $seq_region_end,
        'strand'         => $seq_region_strand,
        'slice'          => $slice,
        '_object_id'     => $row->{object_id},
        'type'           => $row->{type},
        '_phenotype_id'  => $row->{phenotype_id},
        '_phenotype_description'  => $row->{description},
        'adaptor'        => $self,
        '_study_id'      => $row->{study_id},
        '_source_id'     => $row->{source_id},
        '_source_name'   => $row->{name}, 
        'is_significant' => $row->{is_significant},
      }
    );

    $self->{_temp_objs}{$row->{phenotype_feature_id}} = $obj;
  }

  ## add attribs if extracted
  $obj->{attribs}->{$row->{code}} = $row->{value} if $row->{code};

  ## add ontology accession if extracted
  # Gets only the unique ontology accessions in a hash
  $obj->{ontology_accessions_hash}->{$row->{accession}} = 1 if $row->{accession};
  # Then we transform the hash into an array
  my @accessions = keys(%{$obj->{ontology_accessions_hash}});
  $obj->{ontology_accessions} = \@accessions;
}

sub store{
   my ($self, $pf) = @_;

   my $dbh = $self->dbc->db_handle;
   
   # if only object supplied check type matches declared type
   if(! $pf->{_object_id} && ! $pf->object->isa("Bio::EnsEMBL::Variation::$pf->{type}")) {
       throw("Bio::EnsEMBL::Variation::$pf->{type} object expected");
   }

    # look up source_id
    if(defined $pf->source_name && !defined($pf->{source_id})) {
        my $sth = $dbh->prepare(q{
            SELECT source_id FROM source WHERE name = ?
        });
        $sth->execute($pf->source_name);

        my $source_id;
        $sth->bind_columns(\$source_id);
        $sth->fetch();
        $sth->finish();
        $pf->{source_id} = $source_id;
    }

    throw("No source ID found for source name ", $pf->source_name)
        unless (defined($pf->{source_id}) || defined ( $pf->source->dbID));


    # look up submitter from projects such as ClinVar if supplied
    if(defined $pf->{attribs}->{submitter_name} && !defined($pf->{attribs}->{submitter_id})) {
        $pf->{attribs}->{submitter_id} = $self->_get_submitter_id($pf);
    }

    my $sth = $dbh->prepare(q{
        INSERT INTO phenotype_feature (
            phenotype_id,
            source_id,
            study_id,
            type,
            object_id,
            is_significant,
            seq_region_id,
            seq_region_start,
            seq_region_end,
            seq_region_strand            
        ) VALUES (?,?,?,?,?,?,?,?,?,?)
    });

    $sth->execute(
        $pf->phenotype->dbID(),        
        $pf->source ? $pf->source->dbID : $pf->{source_id},
        defined($pf->study)? $pf->study->dbID() : undef,
        $pf->{type},
        defined($pf->{_object_id})? $pf->{_object_id} :  $pf->object->stable_id(),
        defined($pf->{is_significant})? $pf->{is_significant} : 0,
        defined($pf->{slice}) ? $pf->slice()->get_seq_region_id() : undef,
        defined($pf->{start}) ? $pf->{start} :undef,
        defined($pf->{end})   ? $pf->{end} : undef,
        defined($pf->{strand})? $pf->{strand} : undef         
    );
  
   $sth->finish;

   # get dbID
   my $dbID = $dbh->last_insert_id(undef, undef, 'phenotype_feature', 'phenotype_feature_id');
   $pf->{dbID}    = $dbID;
   $pf->{adaptor} = $self;
   
   # add phenotype_feature_attributes
   my $aa = $self->db->get_AttributeAdaptor;
   my $pfa_sth = $dbh->prepare(q{
        INSERT INTO phenotype_feature_attrib (
            phenotype_feature_id,
            attrib_type_id,
            value                  
        ) VALUES (?,?,?)
    });
   
   foreach my $attrib_type( keys %{$pf->{attribs}} ){

       ## we don't store these directly - the id is smaller
       if ($attrib_type eq 'submitter_names'){
         my @submitter_ids;
         my $attrib_type_id = $aa->attrib_id_for_type_code("submitter_id");
         foreach my $submittter_name(@{$pf->{attribs}->{$attrib_type}}){
           my $submitter_id = $self->_get_submitter_id($submittter_name);

           throw("No attrib type ID found for attrib_type submitter_id") unless defined  $attrib_type_id;
            push @submitter_ids, $submitter_id;
         }
         my $sub_ids = join(",", @submitter_ids);
         $pfa_sth->execute( $pf->{dbID},  $attrib_type_id, $sub_ids );

         next;
       }
       my $attrib_type_id = $aa->attrib_id_for_type_code($attrib_type);
       throw("No attrib type ID found for attrib_type $attrib_type") unless defined  $attrib_type_id;
       $pf->{attribs}->{$attrib_type} =~ s/\s+$//;
       $pfa_sth->execute( $pf->{dbID}, $attrib_type_id, $pf->{attribs}->{$attrib_type} );
   }
   $pfa_sth->finish;
}

## Return the id of a submitter from a project such as ClinVar, given the name
## Updates the database if the submitter name is new.
sub _get_submitter_id{

  my $self = shift;
  my $submitter_name = shift;

  my $dbh = $self->dbc->db_handle;

  my $sth = $dbh->prepare(q{
      SELECT submitter_id FROM submitter WHERE description = ?
  });

  $sth->execute($submitter_name);

  my $submitter_id;
  $sth->bind_columns(\$submitter_id);
  $sth->fetch();
  $sth->finish();

  if (! defined $submitter_id){
    ## if it is a new submitter description, enter it.

    my $sth = $dbh->prepare(q{
      INSERT INTO submitter (description) values (?)
    });
    $sth->execute($submitter_name);
    $submitter_id = $dbh->last_insert_id(undef, undef, 'submitter', 'submitter_id');
  }

  return $submitter_id;
}

## Return a submitter name from project such as ClinVar, given the submitter id
## Submitter name s and ids are cached on first access to speed up queries
sub _get_submitter_name{
  my $self         = shift;
  my $submitter_id = shift;

  #load all for speed in extracting a large number
  unless( $self->{submitter_names_lookup}){

    my $dbh = $self->dbc->db_handle;
    my $sth = $dbh->prepare(qq{ SELECT submitter_id, description FROM submitter });
    $sth->execute()||die;
    my $names = $sth->fetchall_arrayref();
    foreach my $n(@{$names}){
      $self->{submitter_names_lookup}->{$n->[0]} = $n->[1];
    }
  }
  return $self->{submitter_names_lookup}->{$submitter_id};

}

1;

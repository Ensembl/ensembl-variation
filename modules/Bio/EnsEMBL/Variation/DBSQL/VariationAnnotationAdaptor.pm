
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VariationAnnotationAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::VariationAnnotationAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $vaa = $reg->get_adaptor("human","variation","variationannotation");
  $va = $reg->get_adaptor("human","variation","variation");
  
  # Get a VariationAnotation by its internal identifier
  $va = $vaa->fetch_by_dbID(45);

  # fetch all annotation for a particular variation
  $v = $va->fetch_by_name('rs56');

  foreach $va (@{$vaa->fetch_all_by_Variation($v)}) {
    print $va->phenotype_name(), $va->phenotype_description(), $va->source_name(), $va->study_type(), $va->local_stable_id(),"\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity between Variation and VariationAnnotation objects.

=head1 AUTHOR - Yuan Chen

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationAnnotationAdaptor;

use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::VariationAnnotation;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;


our @ISA = ('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');

=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL:Variation::Variation $var
  Example    : my @vas = @{$vaa->fetch_all_by_Variation($var)};
  Description: Retrieves all variation annotations for a given variation.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationAnnotation
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

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

  return $self->generic_fetch("va.variation_id = ".$var->dbID());
}

=head2 fetch_all_by_Variation_list

  Arg [1]    : reference to a list of Bio::EnsEMBL:Variation::Variation objects
  Example    : my @vas = @{$vaa->fetch_all_by_Variation_list($vars)};
  Description: Retrieves all variation annotations for a given list of variations
  Returntype : reference to a list of Bio::EnsEMBL::Variation::VariationAnnotation objects
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut


sub fetch_all_by_Variation_list {
  my $self = shift;
  my $vars  = shift;

  if(!ref($vars) || !$vars->[0]->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation arg expected');
  }

  if(!defined($vars->[0]->dbID())) {
    throw("Variation arg must have defined dbID");
  }
  
  my $in_str = join ',', map {$_->dbID()} @$vars;

  return $self->generic_fetch("va.variation_id in (".$in_str.")");
}

=head2 fetch_all_by_VariationFeature_list

  Arg [1]    : reference to a list of Bio::EnsEMBL::Variation::VariationFeature objects
  Example    : my @vas = @{$vaa->fetch_all_by_VariationFeature_list($vfs)};
  Description: Retrieves all variation annotations for a given list of variation features
  Returntype : reference to a list Bio::EnsEMBL::Variation::VariationAnnotation objects
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_VariationFeature_list {
  my $self = shift;
  my $vfs  = shift;

  if(!ref($vfs) || !$vfs->[0]->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
    throw('Listref of Bio::EnsEMBL::Variation::VariationFeature arg expected');
  }

  if(!defined($vfs->[0]->dbID())) {
    throw("VariationFeatures in list must have defined dbIDs");
  }
  
  my $in_str = join ',', map {$_->{'_variation_id'}} @$vfs;

  return $self->generic_fetch("va.variation_id in (".$in_str.")");
}

=head2 fetch_all_by_phenotype_source_name

  Arg [1]    : string $phenotype_name
  Arg [2]    : string $source_name (optional)
  Example    : $vaa = $va_adaptor->fetch_by_phenotype_source_name('BD','EGA');
  Description: Retrieves a variation annotation object via its phenotype/source name
  Returntype : list of ref of Bio::EnsEMBL::Variation::VariationAnnotation
  Exceptions : throw if phenotype name argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_phenotype_name_source_name {

  my $self = shift;
  my $phenotype_name  = shift;
  my $source_name = shift;

  throw('phenotype_name argument expected') if(!defined($phenotype_name));

  my $extra_sql = " p.name = ? ";
  if (defined $source_name ) {
    $extra_sql .= qq( AND s.name = ?);
  }
  
  return $self->generic_fetch("$extra_sql");
  
}

=head2 fetch_all_by_phenotype_description_source_name

  Arg [1]    : string $phenotype_description
  Arg [2]    : string $source_name (optional)
  Example    : $vaa = $va_adaptor->fetch_by_phenotype_description_source_name('diabetes','EGA');
  Description: Retrieves a variation annotation object via its phenotype description/source name
  Returntype : list of ref of Bio::EnsEMBL::Variation::VariationAnnotation
  Exceptions : throw if phenotype name argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_phenotype_description_source_name {

  my $self = shift;
  my $phenotype_description  = shift;
  my $source_name = shift;

  throw('phenotype_description argument expected') if(!defined($phenotype_description));

  my $extra_sql = qq( p.description like '%?%' );
  if (defined $source_name ) {
    $extra_sql .= qq( AND s.name = ?);
  }
  
  return $self->generic_fetch("$extra_sql");
  
}

=head2 fetch_all_by_phenotype_id_source_name

  Arg [1]    : integer $phenotype_id
  Arg [2]    : string $source_name (optional)
  Example    : $vaa = $va_adaptor->fetch_by_phenotype_id_source_name(999,'EGA');
  Description: Retrieves a variation annotation object via its phenotype id/source name
  Returntype : list of ref of Bio::EnsEMBL::Variation::VariationAnnotation
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
    $extra_sql .= sprintf(' AND s.name = %s', $self->dbc->db_handle->quote( $source_name, SQL_VARCHAR ) );
  }
  
  return $self->generic_fetch("$extra_sql");
  
}

# method used by superclass to construct SQL
sub _tables { return (['variation_annotation', 'va'],
                       [ 'phenotype', 'p'],
                       [ 'source', 's']); 
}

sub _default_where_clause {
  my $self = shift;

  return 'va.phenotype_id = p.phenotype_id AND va.source_id = s.source_id';
}

sub _columns {
  return qw( va.variation_annotation_id va.variation_id p.phenotype_id p.name p.description
             s.name va.study va.study_type va.local_stable_id
             va.associated_gene va.associated_variant_risk_allele
	     va.variation_names va.risk_allele_freq_in_controls va.p_value
             );
}



sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @features;

  my ($variation_annotation_id,$var_id,$phenotype_id,$phenotype_name,$phenotype_description,$source_name,$study,$study_type,$local_stable_id,$associated_gene,$associated_variant_risk_allele,$variation_names,$risk_allele_freq_in_controls,$p_value);
  $sth->bind_columns(\$variation_annotation_id,\$var_id,\$phenotype_id,\$phenotype_name,\$phenotype_description,\$source_name,\$study,\$study_type,\$local_stable_id,\$associated_gene,\$associated_variant_risk_allele,\$variation_names,\$risk_allele_freq_in_controls, \$p_value);

  while($sth->fetch()) {
    push @features, $self->_create_feature_fast('Bio::EnsEMBL::Variation::VariationAnnotation',

    {'dbID' => $variation_annotation_id,
     '_variation_id'         => $var_id,
	 '_phenotype_id'		 => $phenotype_id,
     'phenotype_name'        => $phenotype_name,
     'phenotype_description' => $phenotype_description,
     'source_name'           => $source_name,
     'study'                 => $study,
     'study_type'            => $study_type,
     'local_stable_id'       => $local_stable_id,
     'associated_gene'       => $associated_gene,
     'associated_variant_risk_allele' => $associated_variant_risk_allele,
     'variation_names'       => $variation_names,
     'risk_allele_freq_in_controls'   => $risk_allele_freq_in_controls,
     'p_value'               => $p_value,
     'adaptor'  => $self,
    });
  }

  return \@features;


}


1;

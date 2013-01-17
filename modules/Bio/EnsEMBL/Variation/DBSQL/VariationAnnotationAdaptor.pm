=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut


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
    print $va->phenotype_name(), $va->phenotype_description(), $va->source_name(), $va->study_type(),"\n";
  }

=head1 DESCRIPTION

This adaptor provides database connectivity between Variation and VariationAnnotation objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationAnnotationAdaptor;

use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::VariationAnnotation;
use Bio::EnsEMBL::Variation::DBSQL::StudyAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');

=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation::Variation $var
  Example    : my @vas = @{$vaa->fetch_all_by_Variation($var)};
  Description: Retrieves all variation annotations for a given variation.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationAnnotation
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

  return $self->generic_fetch("va.variation_id = ".$var->dbID());
}

=head2 fetch_all_by_Variation_list

  Arg [1]    : reference to a list of Bio::EnsEMBL::Variation::Variation objects
  Example    : my @vas = @{$vaa->fetch_all_by_Variation_list($vars)};
  Description: Retrieves all variation annotations for a given list of variations
  Returntype : reference to a list of Bio::EnsEMBL::Variation::VariationAnnotation objects
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
  Status     : Stable

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


=head2 fetch_all_by_Study

  Arg [1]    : Bio::EnsEMBL:Variation::Study $study
  Example    : my @studies = @{$studya->fetch_all_by_Study($study)};
  Description: Retrieves all variation annotations for a given study.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationAnnotation
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
  return $self->generic_fetch("va.study_id = ".$study->dbID());
}


=head2 fetch_all_by_phenotype_name_source_name

  Arg [1]    : string $phenotype_name
  Arg [2]    : string $source_name (optional)
  Example    : $va = $va_adaptor->fetch_all_by_phenotype_source_name('BD','EGA');
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

  my $extra_sql = " p.name = $phenotype_name ";
  if (defined $source_name ) {
    $extra_sql .= qq( AND s.name = '$source_name' );
  }
  
  # Add the constraint for failed variations
  $extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch("$extra_sql");
  
}

=head2 fetch_all_by_phenotype_description_source_name

  Arg [1]    : string $phenotype_description
  Arg [2]    : string $source_name (optional)
  Example    : $va = $va_adaptor->fetch_all_by_phenotype_description_source_name('diabetes','EGA');
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

  my $extra_sql = qq( p.description like '%$phenotype_description%' );
  if (defined $source_name ) {
    $extra_sql .= qq( AND s.name = '$source_name' );
  }
  
  # Add the constraint for failed variations
  $extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch("$extra_sql");
  
}

=head2 fetch_all_by_phenotype_id_source_name

  Arg [1]    : integer $phenotype_id
  Arg [2]    : string $source_name (optional)
  Example    : $va = $va_adaptor->fetch_all_by_phenotype_id_source_name(999,'EGA');
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
    $extra_sql .= sprintf(" AND s.name = '%s'", $self->dbc->db_handle->quote( $source_name, SQL_VARCHAR ) );
  }
  
  # Add the constraint for failed variations
  $extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch("$extra_sql");
  
}


=head2 fetch_all_by_associated_gene

  Arg [1]    : string $gene_name
  Example    : $va = $va_adaptor->fetch_all_by_associated_gene('CAV3');
  Description: Retrieves the variation annotation objects via which are associated with the gene.
  Returntype : list of ref of Bio::EnsEMBL::Variation::VariationAnnotation
  Exceptions : throw if the gene_name argument is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_associated_gene {

  my $self = shift;
  my $gene_name  = shift;

  throw('gene_name argument expected') if(!defined($gene_name));

	my $extra_sql = " va.associated_gene REGEXP '^(.+,)?[. .]*$gene_name(,.+)?\$'";
  
  # Add the constraint for failed variations
  $extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();

  return $self->generic_fetch($extra_sql);
  
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

  my $extra_sql = " va.associated_gene REGEXP '^(.+,)?[. .]*$gene_name(,.+)?\$'";
  
  # Add the constraint for failed variations
  $extra_sql .= " AND " . $self->db->_exclude_failed_variations_constraint();

  return $self->generic_count($extra_sql);
}


=head2 fetch_all

  Description: Retrieves all available variation annotation objects
  Returntype : list of ref of Bio::EnsEMBL::Variation::VariationAnnotation
  Exceptions : none
  Caller     : general

=cut

sub fetch_all {

  my $self = shift;
  
  # Add the constraint for failed variations
  my $extra_sql = $self->db->_exclude_failed_variations_constraint();
  
  return $self->generic_fetch("$extra_sql");
  
}

# method used by superclass to construct SQL
sub _tables { return (['variation_annotation', 'va'],
                      [ 'failed_variation', 'fv'],
                      [ 'phenotype', 'p'],
											[ 'study', 'st'],
                      [ 'source', 's']); 
}

# Add a left join to the failed_variation table
sub _left_join { return ([ 'failed_variation', 'fv.variation_id = va.variation_id']); }

sub _default_where_clause {
  my $self = shift;

  return 'va.phenotype_id = p.phenotype_id AND st.source_id = s.source_id AND va.study_id=st.study_id';
}

sub _columns {
  return qw( va.variation_annotation_id va.variation_id p.phenotype_id p.name p.description
             va.study_id va.associated_gene va.associated_variant_risk_allele
	         	 va.variation_names va.risk_allele_freq_in_controls va.p_value
           );
}



sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @features;

  my ($variation_annotation_id,$var_id,$phenotype_id,$phenotype_name,$phenotype_description,
      $study_id,$associated_gene,$associated_variant_risk_allele,$variation_names,
	  	$risk_allele_freq_in_controls,$p_value,$last_va_id,$study);
  $sth->bind_columns(\$variation_annotation_id,\$var_id,\$phenotype_id,\$phenotype_name,
                     \$phenotype_description,\$study_id,
					 					 \$associated_gene,\$associated_variant_risk_allele,\$variation_names,
					 					 \$risk_allele_freq_in_controls,\$p_value);
	
	my $sta = $self->db()->get_StudyAdaptor();
	
  while($sth->fetch()) {
    
    next if (defined($last_va_id) && $last_va_id == $variation_annotation_id);
    $last_va_id = $variation_annotation_id;
    
		$study = $sta->fetch_by_dbID($study_id);
		
    push @features, $self->_create_feature_fast('Bio::EnsEMBL::Variation::VariationAnnotation',

    {'dbID' => $variation_annotation_id,
     '_variation_id'         => $var_id,
	   '_phenotype_id'		     => $phenotype_id,
     'phenotype_name'        => $phenotype_name,
     'phenotype_description' => $phenotype_description,
     'associated_gene'       => $associated_gene,
     'associated_variant_risk_allele' => $associated_variant_risk_allele,
     'variation_names'       => $variation_names,
     'risk_allele_freq_in_controls'   => $risk_allele_freq_in_controls,
     'p_value'               => $p_value,
     'adaptor'  => $self,
	   'study' => $study,
    });
  }

  return \@features;

}


=head2 fetch_phenotype_description_by_id
  
  Arg [1]    : int $phenotype_id
  Example    : $phenotype = $va_adaptor->fetch_phenotype_description(10);
  Description: Retrieves the phenotype description from the phenotype ID
  Returntype : string
  Exceptions : throw if the phenotype_id argument is not defined
  Caller     : general

=cut

sub fetch_phenotype_description_by_id {
	my $self = shift;
  my $phenotype_id  = shift;
	
	throw('phenotype_id argument expected') if(!defined($phenotype_id));
	
	my $sth = $self->prepare(qq{SELECT description FROM phenotype WHERE phenotype_id = ?});
	$sth->bind_param(1,$phenotype_id,SQL_INTEGER);
	$sth->execute();
	
	return ($sth->fetchrow_array)[0];
}


=head2 fetch_annotation_number_by_phenotype_id
  
  Arg [1]    : int $phenotype_id
  Example    : $phenotype = $va_adaptor->fetch_annotation_number_by_phenotype_id(10);
  Description: Retrieves the number of variation annotation with the given phenotype ID
  Returntype : integer
  Exceptions : throw if the phenotype_id argument is not defined
  Caller     : general

=cut

sub fetch_annotation_number_by_phenotype_id {
	my $self = shift;
  my $phenotype_id  = shift;
	
	throw('phenotype_id argument expected') if(!defined($phenotype_id));
	
	my $sth = $self->prepare(qq{SELECT count(variation_annotation_id) FROM variation_annotation WHERE phenotype_id = ?});
	$sth->bind_param(1,$phenotype_id,SQL_INTEGER);
	$sth->execute();
	
	return ($sth->fetchrow_array)[0];
}
1;

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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAnnotationAdaptor
#
# Copyright (c) 2011 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAnnotationAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $svaa = $reg->get_adaptor("human","variation","structuralvariationannotation");
	$ssva = $reg->get_adaptor("human","variation","supportingstructuralvariation");
	
	my $ssv = $ssva->fetch_by_name('nssv706165');
	my $sva = $svaa->fetch_all_by_StructuralVariation($ssv);
 

=head1 DESCRIPTION

This adaptor provides database connectivity between StructuralVariation/SupportingStructuralVariation 
and StructuralVariationAnnotation objects.
By default, the 'fetch_all_by_...'-methods will not return structural variants
that have been flagged as failed in the Ensembl QC. This behaviour can be modified
by setting the include_failed_variations flag in Bio::EnsEMBL::Variation::DBSQL::DBAdaptor.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAnnotationAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Variation::BaseStructuralVariation;
use Bio::EnsEMBL::Variation::StructuralVariationAnnotation;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_all_by_StructuralVariation

  Arg [1]    : Bio::EnsEMBL::Variation::StructuralVariation or 
	             Bio::EnsEMBL::Variation::SupportingStructuralVariation $svar
  Example    : my @vas = @{$vaa->fetch_all_by_StructuralVariation($svar)};
  Description: Retrieves all variation annotations for a given variation.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationAnnotation
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_StructuralVariation {
  my $self = shift;
  my $svar = shift;

  if(!ref($svar) || (!$svar->isa('Bio::EnsEMBL::Variation::StructuralVariation') && 
	                   !$svar->isa('Bio::EnsEMBL::Variation::SupportingStructuralVariation'))
	) {
    throw('Bio::EnsEMBL::Variation::StructuralVariation or Bio::EnsEMBL::Variation::SupportingStructuralVariation arg expected');
  }

  if(!defined($svar->dbID())) {
    throw("StructuralVariation arg must have defined dbID");
  }
	
	my $constraint = $self->_internal_exclude_failed_constraint("sva.structural_variation_id = ".$svar->dbID());

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_StructuralVariation_list

  Arg [1]    : reference to a list of Bio::EnsEMBL::Variation::StructuralVariation or 
	             Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Example    : my @svas = @{$svaa->fetch_all_by_StructuralVariation_list($svars)};
  Description: Retrieves all variation annotations for a given list of structural variants
  Returntype : reference to a list of Bio::EnsEMBL::Variation::StructuralVariationAnnotation objects
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_StructuralVariation_list {
  my $self  = shift;
  my $svars = shift;

  if(!ref($svars) || (!$svars->[0]->isa('Bio::EnsEMBL::Variation::StructuralVariation') &&
	                    !$svars->[0]->isa('Bio::EnsEMBL::Variation::SupportingStructuralVariation'))
	) {
		throw('Bio::EnsEMBL::Variation::StructuralVariation or Bio::EnsEMBL::Variation::SupportingStructuralVariation arg expected');
  }

  if(!defined($svars->[0]->dbID())) {
    throw("StructuralVariation arg must have defined dbID");
  }
  
  my $in_str = join ',', map {$_->dbID()} @$svars;

	my $constraint = $self->_internal_exclude_failed_constraint("sva.structural_variation_id in (".$in_str.")");

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_StructuralVariationFeature_list

  Arg [1]    : reference to a list of Bio::EnsEMBL::Variation::StructuralVariationFeature objects
  Example    : my @svas = @{$svaa->fetch_all_by_StructuralVariationFeature_list($svfs)};
  Description: Retrieves all variation annotations for a given list of structural variation features
  Returntype : reference to a list Bio::EnsEMBL::Variation::StructuralVariationAnnotation objects
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_StructuralVariationFeature_list {
  my $self = shift;
  my $svfs = shift;

  if(!ref($svfs) || !$svfs->[0]->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {
    throw('Listref of Bio::EnsEMBL::Variation::StructuralVariationFeature arg expected');
  }

  if(!defined($svfs->[0]->dbID())) {
    throw("VariationFeatures in list must have defined dbIDs");
  }
  
  my $in_str = join ',', map {$_->{'structural_variation_id'}} @$svfs;
	
	my $constraint = $self->_internal_exclude_failed_constraint("sva.structural_variation_id in (".$in_str.")");

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Study

  Arg [1]    : Bio::EnsEMBL:Variation::Study $study
  Example    : my @studies = @{$studya->fetch_all_by_Study($study)};
  Description: Retrieves all structural variation annotations for a given study.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationAnnotation
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
	
	my $constraint = $self->_internal_exclude_failed_constraint('sv.study_id = '.$study->dbID());
  
	return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_sample_name

  Arg [1]    : string $sample_name
  Example    : my $sva = @{$svaa->fetch_all_by_sample_name($sample_name)};
  Description: Retrieves all structural variation annotations for a given sample.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationAnnotation
  Exceptions : throw if sample_name argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_sample_name {
  my $self   = shift;
  my $sample_name  = shift;

	throw('sample_name argument expected') if(!defined($sample_name));
 	
	my $constraint = $self->_internal_exclude_failed_constraint("s1.name = '$sample_name'");
 
	return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_strain_name

  Arg [1]    : string $strain_name
  Example    : my $sva = @{$svaa->fetch_all_by_strain_name($strain_name)};
  Description: Retrieves all structural variation annotations for a given strain.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationAnnotation
  Exceptions : throw if strain_name argument is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_strain_name {
  my $self   = shift;
  my $strain_name  = shift;

	throw('strain_name argument expected') if(!defined($strain_name));
 
 	my $constraint = $self->_internal_exclude_failed_constraint("s2.name = '$strain_name'");
 
	return $self->generic_fetch($constraint);
}


# method used by superclass to construct SQL
sub _tables { 
  my $self = shift;
  my @tables = ([ 'structural_variation_annotation', 'sva'],
	              [ 'study', 'st'],
								[ 'phenotype', 'p'],
								[ 'sample s1', ''],
								[ 'sample s2', ''],
								[ 'structural_variation', 'sv']
							 ); 
							 
	# If we are excluding failed_structural_variations, add that table
  push(@tables,['failed_structural_variation', 'fsv']) unless ($self->db->include_failed_variations());
	
	return @tables;
}

# Add a left join to the failed_variation table
sub _left_join {
  my $self = shift;
  my @tables = ([ 'phenotype', 'p.phenotype_id = sva.phenotype_id'],
								[ 'sample s1', 's1.sample_id = sva.sample_id'],
								[ 'sample s2', 's2.sample_id = sva.strain_id']
							 );
	
	# If we are excluding failed_structural_variations, add that table
  push(@tables,['failed_structural_variation', 'fsv.structural_variation_id=sv.structural_variation_id']) unless ($self->db->include_failed_variations());
  
	return @tables;
}

sub _default_where_clause {
  my $self = shift;

  return 'sva.structural_variation_id = sv.structural_variation_id AND sv.study_id=st.study_id';
}

sub _columns {
  return qw( sva.structural_variation_annotation_id sva.structural_variation_id p.phenotype_id p.description
             sv.study_id sva.clinical_attrib_id s1.name s2.name
           );
}


sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @svas;
  my ($structural_variation_annotation_id,$svar_id,$phenotype_id,$phenotype_description,
      $study_id,$clinical_attrib_id,$sample_name,$strain_name,$study);
  $sth->bind_columns(\$structural_variation_annotation_id,\$svar_id,\$phenotype_id,\$phenotype_description,
                     \$study_id,\$clinical_attrib_id,\$sample_name,\$strain_name);
										 
	my $aa = $self->db()->get_AttributeAdaptor();
	my $sta = $self->db()->get_StudyAdaptor();
	
  while($sth->fetch()) {
    
		$study = $sta->fetch_by_dbID($study_id);
		
		push @svas, Bio::EnsEMBL::Variation::StructuralVariationAnnotation->new(
      -dbID                     => $structural_variation_annotation_id,
      -_STRUCTURAL_VARIATION_ID => $svar_id,
	    -_PHENOTYPE_ID		        => $phenotype_id,
      -PHENOTYPE_DESCRIPTION    => $phenotype_description,
      -SAMPLE_NAME              => $sample_name,
      -STRAIN_NAME              => $strain_name,
      -CLINICAL_SIGNIFICANCE    => $aa->attrib_value_for_id($clinical_attrib_id),
      -ADAPTOR                  => $self,
	    -STUDY                    => $study,
    );
  }

  return \@svas;
}


# Exclude the constraint for failed structural variant
sub _internal_exclude_failed_constraint {
	my $self = shift;
	my $constraint = shift;
	$constraint .= " AND " . $self->db->_exclude_failed_structural_variations_constraint();
	
	return $constraint;
}

1;

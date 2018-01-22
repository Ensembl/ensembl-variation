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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::StructuralVariationSampleAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::StructuralVariationSampleAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $svsa = $reg->get_adaptor("human","variation","structuralvariationsample");
	$ssva = $reg->get_adaptor("human","variation","supportingstructuralvariation");
	
	my $ssv = $ssva->fetch_by_name('nssv706165');
	my $svs = $svsa->fetch_all_by_StructuralVariation($ssv);
 

=head1 DESCRIPTION

This adaptor provides database connectivity between StructuralVariation/SupportingStructuralVariation 
and StructuralVariationSample objects.
By default, the 'fetch_all_by_...'-methods will not return structural variants
that have been flagged as failed in the Ensembl QC. This behaviour can be modified
by setting the include_failed_variations flag in Bio::EnsEMBL::Variation::DBSQL::DBAdaptor.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::StructuralVariationSampleAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Variation::BaseStructuralVariation;
use Bio::EnsEMBL::Variation::StructuralVariationSample;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_all_by_StructuralVariation

  Arg [1]    : Bio::EnsEMBL::Variation::StructuralVariation or 
	             Bio::EnsEMBL::Variation::SupportingStructuralVariation $svar
  Example    : my @vas = @{$vaa->fetch_all_by_StructuralVariation($svar)};
  Description: Retrieves all variation samples for a given variation.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationSample
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
	
	my $constraint = $self->_internal_exclude_failed_constraint("svs.structural_variation_id = ".$svar->dbID());

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_StructuralVariation_list

  Arg [1]    : reference to a list of Bio::EnsEMBL::Variation::StructuralVariation or 
	             Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Example    : my @svss = @{$svsa->fetch_all_by_StructuralVariation_list($svars)};
  Description: Retrieves all variation samples for a given list of structural variants
  Returntype : reference to a list of Bio::EnsEMBL::Variation::StructuralVariationSample objects
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

	my $constraint = $self->_internal_exclude_failed_constraint("svs.structural_variation_id in (".$in_str.")");

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_StructuralVariationFeature_list

  Arg [1]    : reference to a list of Bio::EnsEMBL::Variation::StructuralVariationFeature objects
  Example    : my @svss = @{$svsa->fetch_all_by_StructuralVariationFeature_list($svfs)};
  Description: Retrieves all variation samples for a given list of structural variation features
  Returntype : reference to a list Bio::EnsEMBL::Variation::StructuralVariationSample objects
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
  
  my $in_str = join ',', map {$_->{'_structural_variation_id'}} @$svfs;
	
	my $constraint = $self->_internal_exclude_failed_constraint("svs.structural_variation_id in (".$in_str.")");

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Study

  Arg [1]    : Bio::EnsEMBL:Variation::Study $study
  Example    : my @studies = @{$studya->fetch_all_by_Study($study)};
  Description: Retrieves all structural variation samples for a given study.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationSample
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

=head2 fetch_all_by_Individual

  Arg [1]    : Bio::EnsEMBL:Variation::Individual $individual
  Example    : my $svs = @{$svsa->fetch_all_by_Individual($individual)};
  Description: Retrieves all structural variation samples from a given individual.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationSample
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Individual {
  my $self       = shift;
  my $individual = shift;
  
	if(!ref($individual) || !$individual->isa('Bio::EnsEMBL::Variation::Individual')) {
    throw('Bio::EnsEMBL::Variation::Individual arg expected');
  }
 	
	my $constraint = $self->_internal_exclude_failed_constraint('s.individual_id = '.$individual->dbID());
 
	return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Sample

  Arg [1]    : Bio::EnsEMBL:Variation::Sample $sample
  Example    : my $svs = @{$svsa->fetch_all_by_Sample($sample)};
  Description: Retrieves all structural variation samples from a given sample.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationSample
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Sample {
  my $self   = shift;
  my $sample = shift;
  
	if(!ref($sample) || !$sample->isa('Bio::EnsEMBL::Variation::Sample')) {
    throw('Bio::EnsEMBL::Variation::Sample arg expected');
  }
 	
	my $constraint = $self->_internal_exclude_failed_constraint('s.sample_id = '.$sample->dbID());
 
	return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_strain

  Arg [1]    : Bio::EnsEMBL:Variation::Individual $strain
  Example    : my $svs = @{$svsa->fetch_all_by_strain($strain)};
  Description: Retrieves all structural variation samples from a given strain.
               The strains are stored as an Individual object.
  Returntype : reference to list Bio::EnsEMBL::Variation::StructuralVariationSample
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_strain {
  my $self   = shift;
 
	return $self->fetch_all_by_Individual(@_);
}


# method used by superclass to construct SQL
sub _tables { 
  my $self = shift;
  my @tables = ([ 'structural_variation_sample', 'svs'],
								[ 'structural_variation', 'sv'],
								[ 'sample', 's']
							 ); 
							 
	# If we are excluding failed_structural_variations, add that table
  push(@tables,['failed_structural_variation', 'fsv']) unless ($self->db->include_failed_variations());
	
	return @tables;
}

# Add a left join to the failed_variation table
sub _left_join {
  my $self = shift;
  my @tables = ( ['sample s', 's.sample_id = svs.sample_id'] );
	
	# If we are excluding failed_structural_variations, add that table
  push(@tables,['failed_structural_variation', 'fsv.structural_variation_id=sv.structural_variation_id']) unless ($self->db->include_failed_variations());
  
	return @tables;
}

sub _default_where_clause {
  my $self = shift;

  return 'svs.structural_variation_id = sv.structural_variation_id';
}

sub _columns {
  return qw( svs.structural_variation_sample_id svs.structural_variation_id 
             sv.study_id s.sample_id s.individual_id svs.zygosity
           );
}


sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @svas;
  my ($structural_variation_sample_id,$svar_id,$study_id,$sample_id,$strain_id,$study,$sample,$strain,$zygosity);
  $sth->bind_columns(\$structural_variation_sample_id,\$svar_id,\$study_id,\$sample_id,\$strain_id,\$zygosity);
										 
	my $sample_adapt = $self->db()->get_SampleAdaptor();
	my $ind_adapt    = $self->db()->get_IndividualAdaptor();
	
  while($sth->fetch()) {
    
		$sample = $sample_adapt->fetch_by_dbID($sample_id) if (defined($sample_id));
		
		push @svas, Bio::EnsEMBL::Variation::StructuralVariationSample->new(
      -dbID                     => $structural_variation_sample_id,
      -_STRUCTURAL_VARIATION_ID => $svar_id,
      -SAMPLE                   => $sample,
      -_STRAIN_ID               => $strain_id,
      -ADAPTOR                  => $self,
      -_STUDY_ID                => $study_id,
      -ZYGOSITY                 => $zygosity
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

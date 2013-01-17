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

# Ensembl module for Bio::EnsEMBL::Variation::StructuralVariationAnnotation
#
# Copyright (c) 2011 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::StructuralVariationAnnotation - Annotations for a structural variant (sample and phenotype annotations).

=head1 SYNOPSIS
  
  # A study object
  $study = $study_adaptor->fetch_by_name('nstd37');

  $sva = Bio::EnsEMBL::Variation::StructuralVariationAnnotation->new
        (-sample_name => 'ISCA_ID_5554',
         -clinical_significance => 'Not tested',
         -study => $study);
  ...
  
  $sva->structural_variation->variation_name(),":", $sva->sample_name();      

=head1 DESCRIPTION

This is a class representing the annotation of a structural variant
from the ensembl-variation database. The actual structural variant information is
represented by an associated Bio::EnsEMBL::Variation::StructuralVariation object. 

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::StructuralVariationAnnotation;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Variation::BaseStructuralVariation;
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);

=head2 new

  Arg [-dbID] :
    int - unique internal identifier for variation_annotation
  
  Arg [-ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAnnotationAdaptor
    Adaptor which provides database connectivity for this StructuralVariationAnnotation object
  
  Arg [-_PHENOTYPE_ID] :
    int _ the internal id of the phenotype  
  
  Arg [-PHENOTYPE_DESCRIPTION] :
    string - description of the phenotype
  
  Arg [-SAMPLE_NAME] :
    string - name of the associated sample
  
  Arg [-STRAIN_NAME] :
    string - name of the associated strain
  
  Arg [-CLINICAL_SIGNIFICANCE] :
    string - clinical annotation for this structural variant.
  
  Arg [-STUDY] :
    object ref - the study object describing where the annotated variation comes from.  
  
  Arg [_STRUCTURAL_VARIATION_ID] :
    int _ the internal id of the structural variant object associated with this
    identifier. TUsing this identifier the structural variant may be lazy-loaded from 
    the database on demand.
  
  Example    :
    $study = $study_adaptor->fetch_by_name('nstd37');

    $sva = Bio::EnsEMBL::Variation::StructuralVariationAnnotation->new
          (-sample_name => 'ISCA_ID_5554',
           -strain_name => 'ISCA',
           -clinical_significance => 'Not tested',
           -study => $study);

  Description: Constructor. Instantiates a new StructuralVariationAnnotation object.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariationAnnotation
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($dbID,$adaptor,$phenotype_id,$phenotype_description,$structural_variation_id,$sample_name,
      $strain_name,$clinical_significance,$study) =
    rearrange([qw(dbID ADAPTOR _PHENOTYPE_ID PHENOTYPE_DESCRIPTION _STRUCTURAL_VARIATION_ID 
                  SAMPLE_NAME STRAIN_NAME CLINICAL_SIGNIFICANCE STUDY)],@_); 

  $self->{'dbID'}                     = $dbID;
  $self->{'adaptor'}                  = $adaptor;
  $self->{'_phenotype_id'}            = $phenotype_id;
  $self->{'phenotype_description'}    = $phenotype_description;
  $self->{'_structural_variation_id'} = $structural_variation_id;
  $self->{'sample_name'}              = $sample_name;
  $self->{'strain_name'}              = $strain_name;
  $self->{'clinical_significance'}    = $clinical_significance;
  $self->{'study'}                    = $study;
  return $self;
}



sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 structural_variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::StructuralVariation or 
               Bio::EnsEMBL::Variation::SupportingStructuralVariation $structural_variation
  Example    : $sv = $svf->structural_variation();
  Description: Getter/Setter for the structural variant associated with this feature.
               If not set, and this StructuralVariationFeature has an associated adaptor
               an attempt will be made to lazy-load the structural variation from the
               database.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation or Bio::EnsEMBL::Variation::SupportingStructuralVariation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub structural_variation {
  my $self = shift;

  if(@_) {
    if(!ref($_[0]) || (!$_[0]->isa('Bio::EnsEMBL::Variation::StructuralVariation') &&
                       !$_[0]->isa('Bio::EnsEMBL::Variation::SupportingStructuralVariation')
    )) {
      throw("Bio::EnsEMBL::Variation::StructuralVariation or Bio::EnsEMBL::Variation::SupportingStructuralVariation argument expected");
    }
    $self->{'_structural_variation_id'} = shift;
  }
  elsif(!defined($self->{'structural_variation'}) && $self->{'adaptor'} &&
         defined($self->{'_structural_variation_id'})) {
    # lazy-load from database on demand
    my $sva = $self->{'adaptor'}->db()->get_StructuralVariationAdaptor();
    $self->{'structural_variation'} = $sva->fetch_by_dbID($self->{'_structural_variation_id'});
    if (!defined($self->{'structural_variation'})) {
      $sva = $self->{'adaptor'}->db()->get_SupportingStructuralVariationAdaptor();
      $self->{'structural_variation'} = $sva->fetch_by_dbID($self->{'_structural_variation_id'});
    }
  }

  return $self->{'structural_variation'};
}


=head2 study

  Arg [1]    : Bio::EnsEMBL::Variation::Study (optional)
  Example    : $study = $sv->study()
  Description: Getter/Setter for the study object
  Returntype : Bio::EnsEMBL::Variation::Study
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub study {
  my $self = shift;
  return $self->{'study'} = shift if(@_);
  return $self->{'study'};
}


=head2 study_type  

  Arg [1]    : string study_type (optional)               
               The new value to set the study_type attribute to  
  Example    : $study_type = $obj->study_type()  
  Description: Getter/Setter for the study_type attribute.  
  Returntype : string  
  Exceptions : none  
  Caller     : general  
  Status     : At Risk

=cut

sub study_type{  
  my $self = shift;  
  return $self->{'study'}->type = shift if(@_);  
  return $self->{'study'}->type;
}


=head2 study_name

  Arg [1]    : string $study_name (optional)
               The new value to set the study_name attribute to
  Example    : $study = $sva->study_name()
  Description: Getter/Setter for the study_name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub study_name{
  my $self = shift;
  return $self->{'study'}->name = shift if(@_);
  return $self->{'study'}->name;
}


=head2 study_description

  Arg [1]    : string $study_description (optional)
               The new value to set the study_description attribute to
  Example    : $study_description = $obj->study_description()
  Description: Getter/Setter for the study_description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub study_description{
  my $self = shift;
  return $self->{'study'}->description = shift if(@_);
  return $self->{'study'}->description;
}


=head2 external_reference

  Arg [1]    : string $newval (optional)
               The new value to set the external reference attribute to
  Example    : $external_reference = $obj->external_reference()
  Description: Getter/Setter for the external reference attribute.  This is the
               pubmed/id or project name associated with this study.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub external_reference{
  my $self = shift;
  return $self->{'study'}->external_reference = shift if(@_);
  return $self->{'study'}->external_reference;
}


=head2 study_url

  Arg [1]    : string $newval (optional)
               The new value to set the study_url attribute to
  Example    : $url = $obj->study_url()
  Description: Getter/Setter for the study_url attribute. This is the link to the website where the data are stored.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub study_url{
  my $self = shift;
  return $self->{'study'}->url = shift if(@_);
  return $self->{'study'}->url;
}


=head2 sample_name

  Arg [1]    : string sample_name (optional)
               The new value to set the sample attribute to
  Example    : $sample_name = $obj->sample_name()
  Description: Getter/Setter for the sample attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub sample_name {
  my $self = shift;
  return $self->{'sample_name'} = shift if(@_);
  return $self->{'sample_name'};
}


=head2 strain_name

  Arg [1]    : string strain_name (optional)
               The new value to set the strain attribute to
  Example    : $strain_name = $obj->strain_name()
  Description: Getter/Setter for the strain attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub strain_name {
  my $self = shift;
  return $self->{'strain_name'} = shift if(@_);
  return $self->{'strain_name'};
}


=head2 phenotype_description

  Arg [1]    : string phenotype_description (optional)
               The new value to set the phenotype_description attribute to
  Example    : $phenotype_description = $obj->phenotype_description()
  Description: Getter/Setter for the phenotype_description attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub phenotype_description{
  my $self = shift;
  return $self->{'phenotype_description'} = shift if(@_);
	
  ## Hack to hide wrong phenotypes for COSMIC data ##
	return undef if ($self->{'study'}->description =~ /COSMIC/ && $self->{'phenotype_description'} !~ /^COSMIC/);
	
	return $self->{'phenotype_description'};
}


=head2 clinical_significance

  Arg [1]    : string clinical_significance (optional)
               The new value to set the clinical significance attribute to
  Example    : $clinical_significance = $obj->clinical_significance()
  Description: Getter/Setter for the clinical significance attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub clinical_significance {
  my $self = shift;
  return $self->{'clinical_significance'} = shift if(@_);
  return $self->{'clinical_significance'};
}

1;

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

# Ensembl module for Bio::EnsEMBL::Variation::BaseStructuralVariation
#
#


=head1 NAME

Bio::EnsEMBL::Variation::BaseStructuralVariation - Ensembl representation of a structural variant.

=head1 SYNOPSIS
    # A Study object
    $study = $study_adaptor->fetch_by_name('estd59');
    
    # Structural variation representing a CNV
    $sv = Bio::EnsEMBL::Variation::StructuralVariation->new
       (-variation_name => 'esv234231',
        -class_so_term => ''copy_number_variation',
        -is_evidence => 0,
        -is_somatic => 0,
        -source => $source,
        -study => $study
       );

    ...

    print $sv->variation_name(), ":", $sv->var_class();

=head1 DESCRIPTION

This is a class representing a structural variation from the
ensembl-variation database. A structural variant may have a copy number variation, a tandem duplication, 
an inversion of the sequence or others structural variations. 

The position of a StructuralVariation object on the Genome is represented
by the B<Bio::EnsEMBL::Variation::StructuralVariationFeature> class.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::BaseStructuralVariation;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Variation::Utils::Constants qw(%VARIATION_CLASSES); 
use Bio::EnsEMBL::Variation::Failable;
our @ISA = ('Bio::EnsEMBL::Storable','Bio::EnsEMBL::Variation::Failable');

=head2 new

  Arg [-dbID] :
    see superclass constructor

  Arg [-ADAPTOR] :
    see superclass constructor

  Arg [-VARIATION_NAME] :
    string - the name of the structural variant.
    
  Arg [-CLASS_SO_TERM] :
    string - the sequence ontology term defining the class of the structural variant.
    
  Arg [-SOURCE] :
    object ref - the source object describing where the structural variant comes from.

  Arg [-TYPE] :
     string - the class of structural variant e.g. 'copy_number_variation'
  
  Arg [-STUDY] :
    object ref - the study object describing where the structural variant comes from.
  
  Arg [-VALIDATION_STATUS] :
    string - the status of the structural variant (e.g. validated, not validated, high quality, ...)
  
  Arg [-IS_EVIDENCE] :
    int - flag to inform whether the structural variant is a supporting evidence (1) or not (0).
    
  Arg [-IS_SOMATIC] :
    int - flag to inform whether the structural variant is a somatic (1) or germline (0).
    
  Arg [-ALIAS] :
    string - other name given to the structural variant.
    
  Arg [-CLINICAL_SIGNIFICANCE] :
    reference to list of strings - clinical significance(s) associated with the structural variant, e.g. 'pathogenic'.
    
  Arg [-COPY_NUMBER] :
    int - Number of sequence copies for the supporting evidence of a structural variant classified as "copy number variant" (CNV) when available.

  Example for a structural variation:
    $sv = Bio::EnsEMBL::Variation::StructuralVariation->new
       (-variation_name => 'esv25480',
        -class_so_term => ''copy_number_variation',
        -is_evidence => 0,
        -is_somatic => 0,
        -source => $source,
        -study => $study
       );
    
  Description: Constructor. Instantiates a new structural variant object.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation or 
               Bio::EnsEMBL::Variation::SupportingStructuralVariation
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my (
    $dbID,
    $adaptor,
    $var_name,
    $source_id,
    $source,
    $class_so_term,
    $study_id,
    $study,
    $validation_status,
    $is_evidence,
    $is_somatic,
    $alias,
    $clinical_significance,
    $copy_number
  ) = rearrange([qw(
          dbID
          ADAPTOR
          VARIATION_NAME
          _SOURCE_ID
          SOURCE
          CLASS_SO_TERM
          _STUDY_ID
          STUDY
          VALIDATION_STATUS
          IS_EVIDENCE
          IS_SOMATIC
          ALIAS
          CLINICAL_SIGNIFICANCE
          COPY_NUMBER
    )], @_);
    
  my $self = bless {
    'dbID'                  => $dbID,
    'adaptor'               => $adaptor,
    'variation_name'        => $var_name,
    '_source_id'            => $source_id,
    'source'                => $source,
    'class_SO_term'         => $class_so_term,
    '_study_id'             => $study_id,
    'study'                 => $study,
    'validation_status'     => $validation_status,
    'is_evidence'           => $is_evidence || 0,
    'is_somatic'            => $is_somatic || 0,
    'alias'                 => $alias,
    'clinical_significance' => $clinical_significance,
    'copy_number'           => $copy_number
  };
  return $self;
}



sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 display_id

  Arg [1]    : none
  Example    : print $sv->display_id(), "\n";
  Description: Returns the 'display' identifier for this structural variant.
  Returntype : string
  Exceptions : none
  Caller     : webcode
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  return $self->{'variation_name'} || '';
}



=head2 variation_name

  Arg [1]    : string $newval (optional)
               The new value to set the variation_name attribute to
  Example    : $variation_name = $obj->variation_name()
  Description: Getter/Setter for the variation_name attribute.  This is the
               name of the variation associated with this feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub variation_name{
  my $self = shift;
  return $self->{'variation_name'} = shift if(@_);
  return $self->{'variation_name'};
}

=head2 var_class

    Args        : None
    Example     : my $sv_class = $sv->var_class()
    Description : Getter/setter for the class of structural variant
    ReturnType  : String
    Exceptions  : none
    Caller      : General
    Status      : Stable

=cut

sub var_class {
  my $self = shift;
    
  unless ($self->{class_display_term}) {
    my $display_term = $VARIATION_CLASSES{$self->{class_SO_term}}->{display_term};

    warn "No display term for SO term: ".$self->{class_SO_term} unless $display_term;

    $self->{class_display_term} = $display_term || $self->{class_SO_term};
  }

  return $self->{class_display_term};
}


=head2 class_SO_term

    Args        : None
    Example     : my $sv_so_term = $svf->class_SO_term()
    Description : Getter for the class of structural variant, returning the SO term
    ReturnType  : String
    Exceptions  : none
    Caller      : General
    Status      : Stable

=cut

sub class_SO_term {
  my $self = shift;

  return $self->{class_SO_term};
}


=head2 class_SO_accession

    Args        : None
    Example     : my $sv_so_accession = $svf->class_SO_accession()
    Description : Returns the SO accession corresponding to the structural variant class
    ReturnType  : String
    Exceptions  : none
    Caller      : General
    Status      : Stable

=cut

sub class_SO_accession {
  my $self = shift;

  return $VARIATION_CLASSES{$self->class_SO_term}->{SO_accession};
}


=head2 source

  Arg [1]    : Bio::EnsEMBL::Variation::Source $src (optional)
               The new value to set the source attribute to
  Example    : $source = $sv->source()
  Description: Getter/Setter for the source object attribute
  Returntype : Bio::EnsEMBL::Variation::Source
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source{
  my $self = shift;
  
  # set
  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Source')) {
      throw("Bio::EnsEMBL::Variation::Source argument expected");
    }
    $self->{'source'} = shift;
  }
  # get
  elsif(!defined($self->{'source'}) && $self->adaptor() && defined($self->{'_source_id'})) {
    # lazy-load from database on demand
    my $sa = $self->adaptor->db()->get_SourceAdaptor();
    $self->{'source'} = $sa->fetch_by_dbID($self->{'_source_id'});
  }
  
  return $self->{'source'};
}


=head2 source_name

  Arg [1]    : string $source_name (optional)
               The new value to set the source name attribute to
  Example    : $source_name = $sv->source_name()
  Description: Getter/Setter for the source name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source_name{
  my $self = shift;
  my $source = $self->source;
  return unless defined $source;
  
  $source->name(@_) if(@_);
  return $source->name;
}


=head2 source_version

  Arg [1]    : string $source_version (optional)
               The new value to set the source version attribute to
  Example    : $source_version = $sv->source_version()
  Description: Getter/Setter for the source version attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source_version{
  my $self = shift;
  my $source = $self->source;
  return unless defined $source;
  
  $source->version(@_) if(@_);
  return $source->version;
}


=head2 source_description

  Arg [1]    : string $source_description (optional)
               The new value to set the source description attribute to
  Example    : $source_description = $sv->source_description()
  Description: Getter/Setter for the source description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source_description{
  my $self = shift;
  my $source = $self->source;
  return unless defined $source;
  
  $source->description(@_) if(@_);
  return $source->description;
}


=head2 alias

  Arg [1]    : string $alias (optional)
  Example    : $alias = $sv->alias()
  Description: Getter/Setter for the alias name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub alias {
  my $self = shift;
  return $self->{'alias'} = shift if(@_);
  return $self->{'alias'};
}

=head2 get_all_clinical_significance_states

  Arg [1]    : none
  Example    : my @csstates = @{$sv->get_all_clinical_significance_states()};
  Description: Retrieves all clinical_significance states associated with this structural variant, as reported by dbVar/ClinVar,
               e.g. 'pathogenic','benign', 'drug response',...
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_clinical_significance_states {
    my $self = shift;
    
    return $self->{'clinical_significance'};
}

=head2 validation_status

  Arg [1]    : none
  Example    : my $status = $sv->validation_status();
  Description: Getter/Setter of the validation status for the structural variant. Current
               possible validation statuses are 'validated','not validated',
               'high quality'
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub validation_status {
  my $self = shift;
  return $self->{'validation_status'} = shift if(@_);
  return $self->{'validation_status'};
}


=head2 is_evidence

  Arg [1]    : int $flag (optional)
  Example    : $is_evidence = $obj->is_evidence()
  Description: Getter/Setter of a flag to inform whether the structural variant is a 
               supporting evidence (1) or not (0).
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_evidence{
  my $self = shift;
  return $self->{'is_evidence'} = shift if(@_);
  return $self->{'is_evidence'};
}

=head2 is_somatic

  Arg [1]    : int $flag (optional)
  Example    : $is_somatic = $obj->is_somatic()
  Description: Getter/Setter of a flag to inform whether the structural variant is somatic (1) or germline (0).
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_somatic{
  my $self = shift;
  return $self->{'is_somatic'} = shift if(@_);
  return $self->{'is_somatic'};
}


=head2 copy_number

  Arg [1]    : int $copy (optional)
  Example    : $copy_number = $sv->copy_number()
  Description: Getter/Setter of the number of copies from a structural variant
               classified as 'copy number variant' (CNV) when available.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub copy_number {
  my $self = shift;
  return $self->{'copy_number'} = shift if(@_);
  return $self->{'copy_number'};
}


=head2 study

  Arg [1]    : Bio::EnsEMBL::Variation::Study (optional)
  Example    : $study = $sv->study()
  Description: Getter/Setter for the study object
  Returntype : Bio::EnsEMBL::Variation::Study
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub study {
  my $self = shift;
  
  # set
 if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Study')) {
      throw("Bio::EnsEMBL::Variation::Study argument expected");
    }
    $self->{'study'} = shift;
  }
  # get
  elsif(!defined($self->{'study'}) && $self->adaptor() && defined($self->{'_study_id'})) {
    # lazy-load from database on demand
    my $sa = $self->adaptor->db()->get_StudyAdaptor();
    $self->{'study'} = $sa->fetch_by_dbID($self->{'_study_id'});
  }
  
  return $self->{'study'};
}


sub stable_id {
  my $self = shift;
  return $self->variation_name(@_);
}


=head2 get_all_StructuralVariationFeatures

  Args        : None
  Example     : $svfs = $sv->get_all_StructuralVariationFeatures();
  Description : Retrieves all StructuralVariationFeature for this structural variant
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions  : None
  Caller      : general
  Status      : Stable

=cut

sub get_all_StructuralVariationFeatures{
  my $self = shift;
  
  if(defined $self->{'adaptor'}) {
  
    # get structural variation feature adaptor
    my $svf_adaptor = $self->{'adaptor'}->db()->get_StructuralVariationFeatureAdaptor();
  
    return $svf_adaptor->fetch_all_by_StructuralVariation($self);
  }
  else {
    warn("No variation database attached");
    return [];
  }
}


=head2 get_all_PhenotypeFeatures

  Args        : None
  Example     : $pfs = $sv->get_all_PhenotypeFeatures();
  Description : Retrieves all PhenotypeFeatures for this structural variant
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions  : None
  Caller      : general
  Status      : Stable

=cut

sub get_all_PhenotypeFeatures {
  my $self = shift;
  
  if(defined $self->{'adaptor'}) {
  
    # get phenotype feature adaptor
    my $pfa_adaptor = $self->{'adaptor'}->db()->get_PhenotypeFeatureAdaptor();
  
    return $pfa_adaptor->_fetch_all_by_object($self);
  }
  else {
    warn("No variation database attached");
    return [];
  }
}


=head2 get_all_StructuralVariationSamples

  Args        : None
  Example     : $svss = $sv->get_all_StructuralVariationSamples();
  Description : Retrieves all StructuralVariationSamples for this structural variant
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::StructuralVariationSample
  Exceptions  : None
  Caller      : general
  Status      : Stable

=cut

sub get_all_StructuralVariationSamples {
  my $self = shift;
  
  if(defined $self->{'adaptor'}) {
  
    # get structural variation sample adaptor
    my $svsa_adaptor = $self->{'adaptor'}->db()->get_StructuralVariationSampleAdaptor();
  
    return $svsa_adaptor->fetch_all_by_StructuralVariation($self);
  }
  else {
    warn("No variation database attached");
    return [];
  }
}


=head2 summary_as_hash

  Example       : $sv_summary = $sv->summary_as_hash();
  Description   : Retrieves a textual summary of this StructuralVariation object.
  Returns       : hashref of descriptive strings

=cut

sub summary_as_hash {
  my $self = shift;
  my %summary;
  $summary{'display_id'} = $self->display_id;
  $summary{'study_name'} = $self->study->name;
  $summary{'study_description'} = $self->study->description;
  $summary{'class'} = $self->var_class;
  return \%summary;

}

1;

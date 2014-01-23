=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
        -source => 'DGVa',
        -source_description => 'Database of Genomic Variants Archive',
        -is_evidence => 0,
        -is_somatic => 0);

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
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
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
    string - the name of the source where the structural variant comes from
  
  Arg [-SOURCE_DESCRIPTION] :
    string - description of the source

  Arg [-TYPE] :
     string - the class of structural variant e.g. 'copy_number_variation'
  
  Arg [-STUDY] :
    object ref - the study object describing where the structural variant comes from.
  
  Arg [-VALIDATION_STATUS] :
    string - the status of the structural variant (e.g. validated, not validated, ...)
  
  Arg [-IS_EVIDENCE] :
    int - flag to inform whether the structural variant is a supporting evidence (1) or not (0).
    
  Arg [-IS_SOMATIC] :
    int - flag to inform whether the structural variant is a somatic (1) or germline (0).
    
  Arg [-ALIAS] :
    string - other name given to the structural variant.
    
  Arg [-CLINICAL_SIGNIFICANCE] :
    string - clinical significance associated with the structural variant, e.g. 'Pathogenic'.
    
  Example for a structural variation:
    $sv = Bio::EnsEMBL::Variation::StructuralVariation->new
       (-variation_name => 'esv25480',
        -class_so_term => ''copy_number_variation',
        -source => 'DGVa',
        -source_description => 'Database of Genomic Variants Archive',
        -is_evidence => 0,
        -is_somatic => 0);
    
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
    $source, 
    $source_version, 
    $source_description, 
    $class_so_term,
    $study,
    $validation_status,
    $is_evidence,
    $is_somatic,
    $alias,
    $clinical_significance
  ) = rearrange([qw(
          dbID
          ADAPTOR
          VARIATION_NAME
          SOURCE 
          SOURCE_VERSION
          SOURCE_DESCRIPTION 
          CLASS_SO_TERM
          STUDY
          VALIDATION_STATES
          IS_EVIDENCE
          IS_SOMATIC
          ALIAS
          CLINICAL_SIGNIFICANCE
    )], @_);
    
  my $self = bless {
    'dbID'                  => $dbID,
    'adaptor'               => $adaptor,
    'variation_name'        => $var_name,
    'source'                => $source,
    'source_version'        => $source_version,
    'source_description'    => $source_description,
    'class_SO_term'         => $class_so_term,
    'study'                 => $study,
    'validation_status'     => $validation_status,
    'is_evidence'           => $is_evidence || 0,
    'is_somatic'            => $is_somatic || 0,
    'alias'                 => $alias,
    'clinical_significance' => $clinical_significance
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

    Args         : None
    Example      : my $sv_class = $sv->var_class()
    Description  : Getter/setter for the class of structural variant
    ReturnType   : String
    Exceptions   : none
    Caller       : General
    Status       : Stable

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

    Args         : None
    Example      : my $sv_so_term = $svf->class_SO_term()
    Description  : Getter for the class of structural variant, returning the SO term
    ReturnType   : String
    Exceptions   : none
    Caller       : General
    Status       : Stable

=cut

sub class_SO_term {
  my $self = shift;

  return $self->{class_SO_term};
}


=head2 class_SO_accession

    Args         : None
    Example      : my $sv_so_accession = $svf->class_SO_accession()
    Description  : Returns the SO accession corresponding to the structural variant class
    ReturnType   : String
    Exceptions   : none
    Caller       : General
    Status       : Stable

=cut

sub class_SO_accession {
  my $self = shift;

  return $VARIATION_CLASSES{$self->class_SO_term}->{SO_accession};
}


=head2 source

  Arg [1]    : string $source (optional)
               The new value to set the source attribute to
  Example    : $source = $svf->source()
  Description: Getter/Setter for the source attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source {
  my $self = shift;
  return $self->{'source'} = shift if(@_);
  return $self->{'source'};
}


=head2 source_version

  Arg [1]    : string $source_version (optional)
               The new value to set the source_version attribute to
  Example    : $source_version = $svf->source_version()
  Description: Getter/Setter for the source_version attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub source_version {
  my $self = shift;
  return $self->{'source_version'} = shift if(@_);
  return $self->{'source_version'};
}


=head2 source_description

  Arg [1]    : string $source_description (optional)
               The new value to set the source_description attribute to
  Example    : $source_description = $svf->source_description()
  Description: Getter/Setter for the source_description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub source_description {
  my $self = shift;
  return $self->{'source_description'} = shift if(@_);
  return $self->{'source_description'};
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


=head2 clinical_significance

  Arg [1]    : string $clinical_significance (optional)
  Example    : $clinical_significance = $sv->clinical_significance()
  Description: Getter/Setter for the clinical significance associated with the
               structural variant
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


=head2 get_all_validation_states

  Arg [1]    : none
  Example    : my @vstates = @{$v->get_all_validation_states()};
  Description: Retrieves all validation states for this structural variation. Current
               possible validation statuses are 'validated','not validated',
               'high quality'
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_validation_states {
  my $self = shift;

  return $self->{'validation_status'} || [];
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
  return $self->{'study'} = shift if(@_);
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
  $summary{'study_name'} = $self->study_name;
  $summary{'study_description'} = $self->study_description;
  $summary{'class'} = $self->var_class;
  return \%summary;

}

1;

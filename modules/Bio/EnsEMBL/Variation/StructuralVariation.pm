=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

# Ensembl module for Bio::EnsEMBL::Variation::StructuralVariation
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::StructuralVariation - Ensembl representation of a structural variation.

=head1 SYNOPSIS

    # Structural variation representing a CNV
    $sv = Bio::EnsEMBL::Variation::StructuralVariation->new
       (-variation_name => 'esv25480',
				-class_so_term => 'structural_variant',
				-source => 'DGVa',
				-source_description => 'Database of Genomic Variants Archive',
				-study_name => 'estd20',
				-study_description => 'Conrad 2009 "Origins and functional impact of copy number variation in the human genome." PMID:19812545 [remapped from build NCBI36]',
				-study_url => 'ftp://ftp.ebi.ac.uk/pub/databases/dgva/estd20_Conrad_et_al_2009',
				-external_reference => 'pubmed/19812545');

    ...

    print $sv->name(), ":", $sv->var_class();

=head1 DESCRIPTION

This is a class representing a structural variation from the
ensembl-variation database. A structural variant may have a copy number variation, a tandem duplication, 
an inversion of the sequence or others structural variations. 

The position of a StrucutralVariation object on the Genome is represented
by the B<Bio::EnsEMBL::Variation::StructuralVariationFeature> class.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::StructuralVariation;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Variation::Utils::Constants qw(%VARIATION_CLASSES); 

our @ISA = ('Bio::EnsEMBL::Storable');

=head2 new

  Arg [-dbID] :
    see superclass constructor

  Arg [-ADAPTOR] :
    see superclass constructor

  Arg [-VARIATION_NAME] :
    string - the name of the variation this feature is for (denormalisation
    from Variation object).

	Arg [-CLASS_SO_TERM] :
		string - the sequence ontology term defining the class of the structural variation.
		
  Arg [-SOURCE] :
    string - the name of the source where the variation comes from
	
  Arg [-SOURCE_DESCRIPTION] :
	  string - description of the source

  Arg [-TYPE] :
     string - the class of structural variation e.g. 'copy_number_variation'
	
	Arg [-STUDY] :
    object ref - the study object describing where the annotated variation comes from.
	
	Arg [-VALIDATION_STATUS] :
	  string - the status of the structural variation (e.g. validated, not validated, ...)
	
	
  Example    :
    $sv = Bio::EnsEMBL::Variation::StructuralVariation->new
       (-variation_name => 'esv25480',
		    -class_so_term => 'copy_number_variation',
		    -source => 'DGVa',
		    -source_description => 'Database of Genomic Variants Archive',
		
  Description: Constructor. Instantiates a new StructuralVariation object.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : none
  Caller     : general
  Status     : At Risk

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
		$validation_status
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
    )], @_);
		
	my $self = bless {
		'dbID'               => $dbID,
		'adaptor'            => $adaptor,
  	'variation_name'     => $var_name,
  	'source'             => $source,
  	'source_version'     => $source_version,
  	'source_description' => $source_description,
  	'class_SO_term'      => $class_so_term,
  	'study'              => $study,
		'validation_status'  => $validation_status,
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
  Description: Returns the 'display' identifier for this structural variation. For
               StructuralVariations this is simply the name of the variation
               it is associated with.
  Returntype : string
  Exceptions : none
  Caller     : webcode
  Status     : At Risk

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


=head2 get_all_SupportingStructuralVariants

  Example     : $sv->get_all_SupportingStructuralVariants();
  Description : Retrieves all SupportingStructuralVariation associated with this structural variation.
                Return empty list if there are none.
  Returntype  : reference to list of Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Exceptions  : None
  Caller      : general
  Status      : At Risk

=cut

sub get_all_SupportingStructuralVariants {
	my $self = shift;
	
	if (defined ($self->{'adaptor'})){
		my $ssv_adaptor = $self->{'adaptor'}->db()->get_SupportingStructuralVariationAdaptor();
		return $ssv_adaptor->fetch_all_by_StructuralVariation($self);
    }
    return [];
}


=head2 var_class

    Args         : None
    Example      : my $sv_class = $sv->var_class()
    Description  : Getter/setter for the class of structural variation
    ReturnType   : String
    Exceptions   : none
    Caller       : General
    Status       : At Risk

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
    Description  : Getter for the class of structural variation, returning the SO term
    ReturnType   : String
    Exceptions   : none
    Caller       : General
    Status       : At Risk

=cut

sub class_SO_term {
	my $self = shift;

	return $self->{class_SO_term};
}



=head2 source

  Arg [1]    : string $source (optional)
               The new value to set the source attribute to
  Example    : $source = $svf->source()
  Description: Getter/Setter for the source attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub source{
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

sub source_description{
  my $self = shift;
  return $self->{'source_description'} = shift if(@_);
  return $self->{'source_description'};
}


=head2 get_all_validation_states

  Arg [1]    : none
  Example    : my @vstates = @{$v->get_all_validation_states()};
  Description: Retrieves all validation states for this structural variation.  Current
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


=head2 study_name

  Arg [1]    : string $study (optional)
               The new value to set the study attribute to
  Example    : $study = $sv->study()
  Description: Getter/Setter for the study attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Deprecated

=cut

sub study_name{
  my $self = shift;
	deprecate('Use the method "study" instead (returns a Bio::EnsEMBL::Variation::Study object).');
	return undef if (!$self->study);
  return $self->study->name = shift if(@_);
  return $self->study->name;
}



=head2 study_description

  Arg [1]    : string $study_description (optional)
               The new value to set the study_description attribute to
  Example    : $study_description = $sv->study_description()
  Description: Getter/Setter for the study_description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Deprecated

=cut

sub study_description{
  my $self = shift;
	deprecate('Use the method "study" instead (returns a Bio::EnsEMBL::Variation::Study object).');
	return undef if (!$self->study);
  return $self->study->description = shift if(@_);
  return $self->study->description;
}

=head2 study_url

  Arg [1]    : string $newval (optional)
               The new value to set the study_url attribute to
  Example    : $paper = $obj->study_url()
  Description: Getter/Setter for the study_url attribute.This is the link to the website where the data are stored.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Deprecated

=cut

sub study_url{
  my $self = shift;
	deprecate('Use the method "study" instead (returns a Bio::EnsEMBL::Variation::Study object).');
	return undef if (!$self->study);
  return $self->study->url = shift if(@_);
  return $self->study->url;
}


=head2 external_reference

  Arg [1]    : string $newval (optional)
               The new value to set the external reference attribute to
  Example    : $paper = $obj->external_reference()
  Description: Getter/Setter for the external reference attribute. This is the
               pubmed/id or project name associated with this study.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Deprecated

=cut

sub external_reference{
  my $self = shift;
	deprecate('Use the method "study" instead (returns a Bio::EnsEMBL::Variation::Study object).');
	return undef if (!$self->study);
  return $self->study->external_reference = shift if(@_);
  return $self->study->external_reference;
}


=head2 is_supporting_structural_variation
  Example    : $sv = $obj->is_supporting_structural_variation()
  Description: Getter of the structural variation object for which this structural variant 
	             is a supporting evidence. 
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub is_supporting_structural_variation{
  my $self = shift;

	my $ssva = $self->adaptor->db()->get_SupportingStructuralVariationAdaptor();
	my $ssv  = $ssva->fetch_by_name($self->{'variation_name'});
	if (defined($ssv)) {
		return $ssv->get_StructuralVariation;
	}
	else { return undef; }
}


=head2 get_all_StructuralVariationFeatures

  Args        : None
  Example     : $svfs = $sv->get_all_StructuralVariationFeatures();
  Description : Retrieves all StructuralVariationFeatures for this StructuralVariation
  ReturnType  : reference to list of Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions  : None
  Caller      : general
  Status      : At Risk

=cut

sub get_all_StructuralVariationFeatures{
  my $self = shift;
  
  if(defined $self->{'adaptor'}) {
  
  	# get variation feature adaptor
  	my $svf_adaptor = $self->{'adaptor'}->db()->get_StructuralVariationFeatureAdaptor();
  
  	return $svf_adaptor->fetch_all_by_StructuralVariation($self);
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

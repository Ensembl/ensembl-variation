=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

# Ensembl module for Bio::EnsEMBL::Variation::PhenotypeFeature
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::PhenotypeFeature - a genomic locus associated with a
phenotype.

=head1 SYNOPSIS


=head1 DESCRIPTION

This is a class representing the location of a genomic locus associated with a
phenotype. This region can represent a variant (short or structural), gene, QTL.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::PhenotypeFeature;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Feature;

use Exporter;
use vars qw(@EXPORT_OK @ISA);

our @ISA = ('Bio::EnsEMBL::Feature', 'Exporter');
@EXPORT_OK = qw(%TYPES);

# define valid object types
# this must correspond to the types defined in the type column of
# the phenotype_feature table
our %TYPES = (
  'Variation'                     => 1,
  'StructuralVariation'           => 1,
  'SupportingStructuralVariation' => 1,
  'QTL'                           => 1,
  'Gene'                          => 1,
  'RegulatoryFeature'             => 1,
);

=head2 new

  Arg [-dbID] :
    int - unique internal identifier for variation_annotation
  Arg [-ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::PhenotypeFeatureAdaptor
  Arg [-START] :
    see superclass constructor
  Arg [-END] :
    see superclass constructor
  Arg [-STRAND] :
    see superclass constructor
  Arg [-SLICE] :
    see superclass constructor
  Arg [-PHENOTYPE] :
    Bio::EnsEMBL::Variation::Phenotype
  Arg [-TYPE]
    string - associated object type (e.g. Variation, Gene)
  Arg [-OBJECT]
    object - associated object
  Arg [-SOURCE] :
    string - source of the phenotype association
  Arg [-STUDY_NAME] :
    string - name of study reporting the association
  Arg [-STUDY_DESCRIPTION] :
    string - description of study reporting the association
  Arg [-ATTRIBS] :
	hashref - contains key-value pairs of additional data e.g. p-value, risk
	allele, associated gene
	
  Example    :
  my $pf = Bio::EnsEMBL::Variation::PhenotypeFeature->new(
		-slice     => $slice,
		-start     => 100,
		-end       => 100,
		-phenotype => $phenotype,
		-type      => 'Variation',
		-object    => $variation,
		-source    => 'OMIM',
		-attribs   => {
			p_value => 0.0000023,
		},
  );

  Description: Constructor. Instantiates a new PhenotypeFeature object.
  Returntype : Bio::EnsEMBL::Variation::PhenotypeFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($dbID,$adaptor,$phenotype_id,$phenotype,$type,$object,$object_id,$source,$study,$study_id,$is_significant,$attribs) =
    rearrange([qw(
			dbID ADAPTOR _PHENOTYPE_ID PHENOTYPE
			TYPE OBJECT _OBJECT_ID
			SOURCE STUDY _STUDY_ID
			IS_SIGNIFICANT
			ATTRIBS
		)], @_);

  $self->{'dbID'} = $dbID;
  $self->{'adaptor'} = $adaptor;
  
  # can get phenotype or phenotype ID
  if(defined($phenotype)) {
		$self->{phenotype} = $phenotype;
  }
  elsif(defined($phenotype_id)) {
		$self->{_phenotype_id} = $phenotype_id;
  }
  
  # can get object or object ID
  if(defined($object)) {
		$self->{object} = $object;
  }
  elsif(defined($object_id)) {
		$self->{_object_id} = $object_id;
  }
  
  # can get study or study ID
  if(defined($study)) {
		$self->{study} = $study;
  }
  elsif(defined($study_id)) {
		$self->{_study_id} = $study_id;
  }
  
  $self->{type}           = $type;
  $self->{source}         = $source;
  $self->{is_significant} = $is_significant;
  $self->{attribs}        = $attribs || {};
  
  return $self;
}

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 phenotype

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Phenotype $ph
  Example    : $ph = $pf->phenotype();
  Description: Getter/Setter for the phenotype associated with this annotation.
               If not set, and this PhenotypeFeature has an associated adaptor
               an attempt will be made to lazy-load the phenotype from the
               database.
  Returntype : Bio::EnsEMBL::Variation::Phenotype
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub phenotype {
  my ($self, $phenotype) = @_;
  
  # set
  $self->{phenotype} = $phenotype if defined($phenotype);
  
  # get
  if(!defined($self->{phenotype}) && defined($self->{_phenotype_id})) {
		my $pa = $self->adaptor->db->get_PhenotypeAdaptor();
		
		$self->{phenotype} = $pa->fetch_by_dbID($self->{_phenotype_id});
  }
  
  return $self->{phenotype};
}


=head2 object

  Arg [1]    : (optional) Bio::EnsEMBL::* object $ph
  Example    : $object = $pf->object();
  Description: Getter/Setter for the object associated with this annotation.
	             PhenotypeFeatures may be associated with several Ensembl object
							 types e.g. Variation, StruturalVariation, Gene.
               If not set, and this PhenotypeFeature has an associated adaptor
               an attempt will be made to lazy-load the object from the
               database.
  Returntype : Bio::EnsEMBL::*
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub object {
  my ($self, $object) = @_;
  
  # set
	if(defined($object)) {
		my $type = (split '::', ref($object))[-1];
		throw("$type is not a valid object type, valid types are: ".(join ", ", sort %TYPES)) unless defined $type and defined($TYPES{$type});
		
		$self->{object} = $object;
		
		# update type
		$self->type($type);
	}
  
  # get
  if(!defined($self->{object})) {
		throw("No object or internal identifier found for PhenotypeFeature") unless defined($self->{_object_id});
		
		# get object type and correct adaptor
		my $type = $self->type;
		my $method = 'get_'.$type.'Adaptor';
		my $adaptor = $self->adaptor->db->$method || $self->adaptor->db->dnadb->$method;
		
		# fetch the object
		$self->{object} = $adaptor->fetch_by_stable_id($self->{_object_id});
  }
  
  return $self->{object};
}


=head2 object_id

  Arg [1]    : (optional) string $object_id
  Example    : $object_id = $pf->object_id();
  Description: Getter/Setter for the ID of the object associated with this
	             annotation.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub object_id {
	my ($self, $object_id) = @_;
	
	$self->{_object_id} = $object_id if defined($object_id);
	return $self->{_object_id};
}


=head2 variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Variation $variation
  Example    : $v = $pf->variation();
  Description: Getter/Setter for the variation associated with this annotation.
               If not set, and this PhenotypeFeature has an associated adaptor
               an attempt will be made to lazy-load the variation from the
               database. Can only be called when $pf->type() is 'Variation'; for
							 PhenotypeFeatures with other object types, use $pf->object();
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub variation {
  my $self = shift;

  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Variation')) {
      throw("Bio::EnsEMBL::Variation::Variation argument expected");
    }
    $self->{'variation'} = shift;
  }
  elsif(!defined($self->{'variation'}) && $self->{'adaptor'} &&
        defined($self->{'_variation_id'})) {
    # lazy-load from database on demand
    my $va = $self->{'adaptor'}->db()->get_VariationAdaptor();
    $self->{'variation'} = $va->fetch_by_dbID($self->{'_variation_id'});
  }

  return $self->{'variation'};
}


=head2 type

  Arg [1]    : string $type (optional)
               The new value to set the type attribute to
  Example    : $type = $obj->type()
  Description: Getter/Setter for the object type of the PhenotypeFeature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub type {
  my $self = shift;
	my $type = shift;
	
	if(defined($type)) {
		throw("$type is not a valid object type, valid types are: ".(join ", ", sort %TYPES)) unless defined($TYPES{$type});
		$self->{'type'} = $type;
	}
  
  return $self->{'type'};
}


=head2 is_significant

  Arg [1]    : boolean $is_significant(optional)
               The new value to set the is_significant attribute to
  Example    : $is_significant = $obj->is_significant()
  Description: Getter/Setter for the is_significant attribute - identifies
	             whether this phenotype association should be considered as
							 significant.
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_significant {
  my $self = shift;
  return $self->{'is_significant'} = shift if(@_);
  return $self->{'is_significant'};
}


=head2 source

  Arg [1]    : string source (optional)
               The new value to set the source attribute to
  Example    : $source = $obj->source()
  Description: Getter/Setter for the source attribute.
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


=head2 study

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Study $st
  Example    : $st = $pf->study();
  Description: Getter/Setter for the study associated with this annotation.
               If not set, and this PhenotypeFeature has an associated adaptor
               an attempt will be made to lazy-load the study from the
               database.
  Returntype : Bio::EnsEMBL::Variation::Study
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub study {
  my ($self, $study) = @_;
  
  # set
  $self->{study} = $study if defined($study);
  
  # get
  if(!exists($self->{study})) {
		if(!defined($self->{_study_id})) {
      $self->{study} = undef;
    }
		
    else {
      my $pa = $self->adaptor->db->get_StudyAdaptor();
      
      $self->{study} = $pa->fetch_by_dbID($self->{_study_id});
    }
  }
  
  return $self->{study};
}


=head2 study_name

  Arg [1]    : string $study_name (optional)
               The new value to set the study_name attribute to
  Example    : $study = $sva->study_name()
  Description: Getter/Setter for the study_name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub study_name {
  my $self = shift;
	
	my $study = $self->study;
	return unless defined $study;
	
  $study->name(@_) if(@_);
  return $study->name;
}


=head2 study_description

  Arg [1]    : string $study_description (optional)
               The new value to set the study_description attribute to
  Example    : $study_description = $obj->study_description()
  Description: Getter/Setter for the study_description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub study_description {
  my $self = shift;
	
	my $study = $self->study;
	return unless defined $study;
	
  $study->description(@_) if(@_);
  return $study->description;
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

sub external_reference {
  my $self = shift;
	
	my $study = $self->study;
	return unless defined $study;
	
  $study->external_reference(@_) if(@_);
  return $study->external_reference;
}


=head2 study_url

  Arg [1]    : string $newval (optional)
               The new value to set the study_url attribute to
  Example    : $url = $obj->study_url()
  Description: Getter/Setter for the study_url attribute. This is the link to
	             the website where the data are stored.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub study_url {
  my $self = shift;
	
	my $study = $self->study;
	return unless defined $study;
	
  $study->url(@_) if(@_);
  return $study->url;
}


=head2 associated_studies
  Example    : $name = $obj->associate_studies()
  Description: Getter/Setter for the associated_studies attribute 
	            (e.g. EGA studies can be associated to NHGRI studies). 
  Returntype : reference to list of Bio::EnsEMBL::Variation::Study
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub associated_studies {
  my $self = shift;
	return $self->study ? $self->study->associated_studies : undef;
}



=head2 get_all_attributes

  Example    : %attribs = %{$obj->get_all_attributes}
  Description: Retrieves attributes of this PhenotypeFeature as a hash reference
               containing key-value pairs e.g. "p_value" => 0.0000012
  Returntype : hashref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_attributes {
  my $self = shift;
  
  if(!defined($self->{attribs})) {
	$self->{attribs} = $self->adaptor->_fetch_attribs_by_dbID($self->dbID);
  }
  
  return $self->{attribs};
}


sub _set_attribute {
	my $self  = shift;
	my $key   = shift;
	my $value = shift;
	
	$self->get_all_attributes;
	$self->{attribs}->{$key} = $value;
}


=head2 variation_names

  Arg [1]    : string $newval (optional)
               The new value to set the variation_names attribute to
  Example    : $variation_names = $obj->variation_names()
  Description: Getter/Setter for the variation_names attribute.  This is the
               variation name(s) linked with this phenotype association.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub variation_names {
  my $self = shift;
	my $new  = shift;
	
	$self->_set_attribute('variation_names', $new) if defined($new);
	
  return defined($self->get_all_attributes->{'variation_names'}) ? $self->get_all_attributes->{'variation_names'} : undef;
}


=head2 associated_gene

  Arg [1]    : string $newval (optional)
               The new value to set the associated_gene attribute to
  Example    : $associated_gene = $obj->associated_gene()
  Description: Getter/Setter for the associated_gene attribute.  This is the
               gene name(s) linked with this phenotype association.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub associated_gene {
  my $self = shift;
	my $new  = shift;
	
	$self->_set_attribute('associated_gene', $new) if defined($new);
	
  return defined($self->get_all_attributes->{'associated_gene'}) ? $self->get_all_attributes->{'associated_gene'} : undef;
}


=head2 risk_allele

  Example    : $risk_allele = $obj->risk_allele()
  Description: Getter/setter for the risk_allele attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub risk_allele {
  my $self = shift;
	my $new  = shift;
	
	$self->_set_attribute('risk_allele', $new) if defined($new);
	
  return defined($self->get_all_attributes->{'risk_allele'}) ? $self->get_all_attributes->{'risk_allele'} : undef;
}


=head2 p_value

  Example    : $p_value = $obj->p_value()
  Description: Getter/setter for the p_value attribute.
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub p_value {
  my $self = shift;
	my $new  = shift;
	
	$self->_set_attribute('p_value', $new) if defined($new);
	
  return defined($self->get_all_attributes->{'p_value'}) ? $self->get_all_attributes->{'p_value'} : undef;
}

=head2 clinical_significance

  Example    : $clinical_significance = $obj->clinical_significance()
  Description: Getter/setter for the clinical_significance attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub clinical_significance {
	my $self = shift;
	my $new  = shift;
	
	$self->_set_attribute('dgva_clin_sig', $new) if defined($new);
	
	return defined($self->get_all_attributes->{'dgva_clin_sig'}) ? $self->get_all_attributes->{'dgva_clin_sig'} : undef;
}


=head2 external_id

  Example    : $external_id = $obj->external_id()
  Description: Getter/setter for the external_id attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub external_id {
	my $self = shift;
	my $new  = shift;
	
	$self->_set_attribute('external_id', $new) if defined($new);
	
	return defined($self->get_all_attributes->{'external_id'}) ? $self->get_all_attributes->{'external_id'} : undef;
}


=head2 sample_name

  Example    : $sample_name = $obj->sample_name()
  Description: Getter for the sample_name attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub sample_name {
	my $self = shift;
	
	if(!exists($self->{sample_name})) {
		if(my $sample_id = $self->get_all_attributes->{'sample_id'}) {
			$self->{sample_name} = $self->_generic_sample_name_fetch($sample_id);
		}
		else {
			$self->{sample_name} = undef;
		}
	}
	
	return $self->{sample_name};
}


=head2 strain_name

  Example    : $strain_name = $obj->strain_name()
  Description: Getter for the strain_name attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub strain_name {
	my $self = shift;
	
	if(!exists($self->{strain_name})) {
		if(my $strain_id = $self->get_all_attributes->{'strain_id'}) {
			$self->{strain_name} = $self->_generic_sample_name_fetch($strain_id);
		}
		else {
			$self->{strain_name} = undef;
		}
	}
	
	return $self->{strain_name};
}

# internal method used by strain_name and sample_name
sub _generic_sample_name_fetch {
	my $self = shift;
	my $id = shift;
	
	# method we want is on SampleAdaptor, but that is a base adaptor so we can't
	# use it directly
	my $ia = $self->adaptor->db->get_IndividualAdaptor;
	return $ia->_get_sample_name_by_dbID($id);
}


1;

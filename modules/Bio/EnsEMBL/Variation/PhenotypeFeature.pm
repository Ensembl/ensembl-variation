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

# Ensembl module for Bio::EnsEMBL::Variation::PhenotypeFeature
#
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

use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
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
  Arg [-SOURCE_VERSION] :
    string - version of the source of the phenotype association
  Arg [-SOURCE_OBJECT] :
    object  - source of the phenotype association
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

  my ($dbID,$adaptor,$phenotype_id,$phenotype,$type,$object,$object_id,$source_name,$source_id,$source,$study,$study_id,$is_significant,$attribs, $ontology_accessions) =
    rearrange([qw(
      dbID ADAPTOR _PHENOTYPE_ID PHENOTYPE
      TYPE OBJECT _OBJECT_ID
      SOURCE_NAME _SOURCE_ID SOURCE STUDY _STUDY_ID
      IS_SIGNIFICANT
      ATTRIBS ONTOLOGY_ACCESSIONS
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
  
  # can get source or source ID
  if(defined($source)) {
    $self->{source} = $source;
  }
  elsif(defined($source_id)) {
    $self->{_source_id} = $source_id;
  }
   elsif(defined($source_name)) {
    $self->{source_name} = $source_name;
  }

  # can get study or study ID
  if(defined($study)) {
    $self->{study} = $study;
  }
  elsif(defined($study_id)) {
    $self->{_study_id} = $study_id;
  }


  $self->{type}                = $type;
  $self->{is_significant}      = $is_significant;
  $self->{attribs}             = $attribs || {};
  $self->{ontology_accessions} = $ontology_accessions || undef;

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


=head2 phenotype_description

  Example    : $desc = $pf->phenotype_description();
  Description: Convenience method to get the phenotype description
               associated with this annotation.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : experimental

=cut

sub phenotype_description {
  my $self = shift;

  return $self->{_phenotype_description} if $self->{_phenotype_description};

  if(!defined($self->{phenotype}) && defined($self->{_phenotype_id})) {
    my $pa = $self->adaptor->db->get_PhenotypeAdaptor();

    $self->{phenotype} = $pa->fetch_by_dbID($self->{_phenotype_id});
  }
  $self->{_phenotype_description} =  $self->{phenotype}->description();

  return $self->{_phenotype_description};
}


=head2 phenotype_id

  Example    : $id = $pf->phenotype_id();
  Description: Convenience method to get the phenotype internal ID
               associated with this annotation.
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : experimental

=cut

sub phenotype_id {
  my $self = shift;

  return $self->{_phenotype_id} if $self->{_phenotype_id};
  
  if (defined($self->{phenotype})) {
    if ($self->{phenotype}->dbID) {
      $self->{_phenotype_id} = $self->{phenotype}->dbID;
      return $self->{phenotype}->dbID;
    }
  }
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
    my $adaptor;
    if ($type eq 'Gene') {
      $adaptor = $self->adaptor->db->dnadb->$method;
    } else {
      $adaptor = $self->adaptor->db->$method;
    }
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
        defined($self->{'_object_id'}) && $self->{'type'} eq 'Variation') {
    # lazy-load from database on demand
    my $va = $self->{'adaptor'}->db()->get_VariationAdaptor();
    $self->{'variation'} = $va->fetch_by_name($self->{'_object_id'});
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

  Arg [1]    : Bio::EnsEMBL::Variation::Source $src (optional)
               The new value to set the source attribute to
  Example    : $source = $pf->source()
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
  Example    : $source_name = $pf->source_name()
  Description: Getter/Setter for the source name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source_name{
  my $self = shift;

  return $self->{'_source_name'} if $self->{'_source_name'};

  my $source = $self->source;
  return unless defined $source;
  
  $source->name(@_) if(@_);
  return $source->name;
}


=head2 source_version

  Arg [1]    : string $source_version (optional)
               The new value to set the source version attribute to
  Example    : $source_version = $pf->source_version()
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
               Here is a list of the keys used: associated_gene, beta_coef, clinvar_clin_sig, external_id,
               inheritance_type, odds_ratio, p_value, review_status, risk_allele, variation_names
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
  ## standardise clinical significance terms to lowercase
  $self->{attribs}->{'clinvar_clin_sig'} = "\L$self->{attribs}->{'clinvar_clin_sig'}"
     if defined $self->{attribs}->{'clinvar_clin_sig'};

  return $self->{attribs};
}


=head2 get_all_ontology_accessions

  Example    : @ontology_acc = @{$obj->get_all_ontology_accessions}
  Description: Retrieves all ontology accessions linked to the phenotype entry of this PhenotypeFeature
  Returntype : listref of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_ontology_accessions {
  my $self = shift;

  if(!defined($self->{ontology_accessions})) {
    $self->{ontology_accessions} = $self->phenotype->ontology_accessions;
  }
  return $self->{ontology_accessions};
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

  $self->_set_attribute('clinvar_clin_sig', $new) if defined($new);
  
  return defined($self->get_all_attributes->{'clinvar_clin_sig'}) ? $self->get_all_attributes->{'clinvar_clin_sig'} : undef;
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

=head2 allele_symbol

  Example    : $allele_symbol = $obj->allele_symbol()
  Description: Getter/setter for allele_symbol attribute. This is only stored for mouse phenotype data.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub allele_symbol {
  my $self = shift;
  my $new  = shift;
  
  $self->_set_attribute('allele_symbol', $new) if defined($new);
  
  return defined($self->get_all_attributes->{'allele_symbol'}) ? $self->get_all_attributes->{'allele_symbol'} : undef;
}

=head2 allele_accession_id

  Example    : $allele_accession_id = $obj->allele_accession_id()
  Description: Getter/setter for allele_accession_id attribute. This is only stored for mouse phenotype data.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub allele_accession_id {
  my $self = shift;
  my $new  = shift;
  
  $self->_set_attribute('allele_accession_id', $new) if defined($new);
  
  return defined($self->get_all_attributes->{'allele_accession_id'}) ? $self->get_all_attributes->{'allele_accession_id'} : undef;
}

=head2 marker_accession_id

  Example    : $marker_accession_id = $obj->marker_accession_id()
  Description: Getter/setter for marker_accession_id attribute. This is only stored for mouse phenotype data.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub marker_accession_id {
  my $self = shift;
  my $new  = shift;
  
  $self->_set_attribute('marker_accession_id', $new) if defined($new);
  
  return defined($self->get_all_attributes->{'marker_accession_id'}) ? $self->get_all_attributes->{'marker_accession_id'} : undef;
}

=head2 pipeline_name

  Example    : $pipeline_name = $obj->pipeline_name()
  Description: Getter/setter for pipeline_name attribute. This is only stored for mouse phenotype data.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub pipeline_name {
  my $self = shift;
  my $new  = shift;
  
  $self->_set_attribute('pipeline_name', $new) if defined($new);
  
  return defined($self->get_all_attributes->{'pipeline_name'}) ? $self->get_all_attributes->{'pipeline_name'} : undef;
}

=head2 procedure_name

  Example    : $procedure_name = $obj->procedure_name()
  Description: Getter/setter for procedure_name attribute. This is only stored for mouse phenotype data.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub procedure_name {
  my $self = shift;
  my $new  = shift;
  
  $self->_set_attribute('procedure_name', $new) if defined($new);
  
  return defined($self->get_all_attributes->{'procedure_name'}) ? $self->get_all_attributes->{'procedure_name'} : undef;
}

=head2 parameter_name

  Example    : $parameter_name = $obj->parameter_name()
  Description: Getter/setter for parameter_name attribute. This is only stored for mouse phenotype data.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub parameter_name {
  my $self = shift;
  my $new  = shift;
  
  $self->_set_attribute('parameter_name', $new) if defined($new);
  
  return defined($self->get_all_attributes->{'parameter_name'}) ? $self->get_all_attributes->{'parameter_name'} : undef;
}

=head2 project_name

  Example    : $project_name = $obj->project_name()
  Description: Getter/setter for project_name attribute. This is only stored for mouse phenotype data.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub project_name {
  my $self = shift;
  my $new  = shift;
  
  $self->_set_attribute('project_name', $new) if defined($new);
  
  return defined($self->get_all_attributes->{'project_name'}) ? $self->get_all_attributes->{'project_name'} : undef;
}

=head2 project_fullname

  Example    : $project_fullname = $obj->project_fullname()
  Description: Getter/setter for project_fullname attribute. This is only stored for mouse phenotype data.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub project_fullname {
  my $self = shift;
  my $new  = shift;
  
  $self->_set_attribute('project_fullname', $new) if defined($new);
  
  return defined($self->get_all_attributes->{'project_fullname'}) ? $self->get_all_attributes->{'project_fullname'} : undef;
}

=head2 strain

  Example    : $strain = $pf->strain();
  Description: Getter/Setter for the strain associated with this annotation.
               If not set, and this PhenotypeFeature has an associated adaptor
               an attempt will be made to lazy-load the variation from the
               database. 
  Returntype : Bio::EnsEMBL::Variation::Individual
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub strain {
  my $self = shift;
  if (@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Individual')) {
      throw("Bio::EnsEMBL::Variation::Variation argument expected");
    } 
    $self->{'strain'} = shift;
  } elsif ( !defined($self->{'strain'}) && $self->{'adaptor'} ) {
    $self->{'strain_id'} = $self->get_all_attributes->{'strain_id'};  
  # lazy-load from database on demand
    if (!$self->{'strain_id'}) {
      my $sample_adaptor = $self->{'adaptor'}->db()->get_SampleAdaptor();
      my $reference_strain = $sample_adaptor->fetch_reference_strain;
      $self->{'strain_id'} = $reference_strain->individual->dbID;
    }
    my $ia = $self->{'adaptor'}->db()->get_IndividualAdaptor();
    $self->{'strain'} = $ia->fetch_by_dbID($self->{'strain_id'});
  } else {
    throw("Adaptor is not defined.");
  }

  return $self->{'strain'};
}

=head2 odds_ratio

  Example    : $odds_ratio = $obj->odds_ratio()
  Description: Getter/setter for odds_ratio attribute. This is only stored for human phenotype data.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub odds_ratio {
  my $self = shift;
  my $new  = shift;
  
  $self->_set_attribute('odds_ratio', $new) if defined($new);
  
  return defined($self->get_all_attributes->{'odds_ratio'}) ? $self->get_all_attributes->{'odds_ratio'} : undef;
}

=head2 beta_coefficient

  Example    : $beta_coef = $obj->beta_coefficient()
  Description: Getter/setter for beta_coef attribute. This is only stored for human phenotype data.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub beta_coefficient {
  my $self = shift;
  my $new  = shift;
  
  $self->_set_attribute('beta_coef', $new) if defined($new);
  
  return defined($self->get_all_attributes->{'beta_coef'}) ? $self->get_all_attributes->{'beta_coef'} : undef;
}

=head2 submitter_names

  Example    : $names = $obj->submitter_names()
  Description: Get all submitter_names. Eg submitters to ClinVar.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub submitter_names {
  my $self = shift;
  return defined($self->get_all_attributes->{'submitter_names'}) ? $self->get_all_attributes->{'submitter_names'} : undef;
}

=head2 date_last_evaluated

  Example    : $date = $obj->date_last_evaluated()
  Description: Get the date evidence for the assertion was last evaluated. Eg for ClinVar.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub date_last_evaluated{
  my $self = shift;

  return defined($self->get_all_attributes->{'DateLastEvaluated'}) ? $self->get_all_attributes->{'DateLastEvaluated'} : undef;
}


1;

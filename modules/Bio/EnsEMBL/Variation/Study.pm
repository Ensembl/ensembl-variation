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

# Ensembl module for Bio::EnsEMBL::Variation::Study
#
# Copyright (c) 2011 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::Study - Ensembl representation of a study.

=head1 SYNOPSIS

    # Structural variation feature representing a CNV
    $ssv = Bio::EnsEMBL::Variation::StructuralVariation->new
       (-supporting_structural_evidence => 'ssv001'
        -structural_variation   => $structural_variation);

    ...


=head1 DESCRIPTION

This is a class representing the supporting evidence of a structural variation
from the ensembl-variation database.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Study;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);

our @ISA = ('Bio::EnsEMBL::Storable');


=head2 new

  Arg [-dbID] :
    see superclass constructor

  Arg [-ADAPTOR] :
    see superclass constructor

  Arg [-NAME] :
    see superclass constructor
    
  Arg [-STRUCTURAL_VARIATION_ID] :
    see superclass constructor

  Example    :
		
    $svv = Bio::EnsEMBL::Variation::SupportingStructuralVariation->new
       (-name => 'esv25480',
        -structural_variation_id   => $structural_variation->dbID);

  Description: Constructor. Instantiates a new StructuralVariation object.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
	my ($dbID,$adaptor,$study_name,$study_description,$study_url,$external_reference,
	    $study_type,$source_name,$associate) = 
			rearrange([qw(dbID ADAPTOR NAME DESCRIPTION URL EXTERNAL_REFERENCE TYPE SOURCE ASSOCIATE)], @_);

  $self = {
			'dbID' => $dbID,
			'adaptor' => $adaptor,
  		'name' => $study_name,
			'description' => $study_description,
			'url' => $study_url,
			'external_reference' => $external_reference,
			'type' => $study_type,
			'source' => $source_name,
			'associate' => $associate
	};
	
	return bless $self, $class;
}


=head2 name

  Arg [1]    : string $newval (optional)
               The new value to set the name attribute to
  Example    : $name = $obj->name()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub name{
  my $self = shift;
  return $self->{'name'} = shift if(@_);
  return $self->{'name'};
}


=head2 description

  Arg [1]    : string $newval (optional)
               The new value to set the description attribute to
  Example    : $name = $obj->description()
  Description: Getter/Setter for the description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub description{
  my $self = shift;
  return $self->{'description'} = shift if(@_);
  return $self->{'description'};
}


=head2 url

  Arg [1]    : string $newval (optional)
               The new value to set the url attribute to
  Example    : $name = $obj->url()
  Description: Getter/Setter for the url attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub url{
  my $self = shift;
  return $self->{'url'} = shift if(@_);
  return $self->{'url'};
}


=head2 external_reference

  Arg [1]    : string $newval (optional)
               The new value to set the external_reference attribute to
  Example    : $name = $obj->external_reference()
  Description: Getter/Setter for the external_reference attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub external_reference{
  my $self = shift;
  return $self->{'external_reference'} = shift if(@_);
  return $self->{'external_reference'};
}


=head2 type

  Arg [1]    : string $newval (optional)
               The new value to set the type attribute to
  Example    : $name = $obj->type()
  Description: Getter/Setter for the type attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub type{
  my $self = shift;
  return $self->{'type'} = shift if(@_);
  return $self->{'type'};
}


=head2 source

  Arg [1]    : string $newval (optional)
               The new value to set the source attribute to
  Example    : $name = $obj->source()
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


=head2 associated_studies
  Example    : $name = $obj->associate_studies()
  Description: Getter/Setter for the associated_studies attribute
  Returntype : reference to list of Bio::EnsEMBL::Variation::Study
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub associated_studies{
  my $self = shift;
	
	my $results;
	
  if (defined($self->{'associate'}) && defined($self->{'adaptor'})) {
		my $studya = $self->{'adaptor'}->db()->get_StudyAdaptor();
		return $studya->fetch_all_by_dbID_list($self->{'associate'});
	}
	else {
  	return [];
	}
}

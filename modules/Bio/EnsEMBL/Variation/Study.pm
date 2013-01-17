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

# Ensembl module for Bio::EnsEMBL::Variation::Study
#
# Copyright (c) 2011 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::Study - Ensembl representation of a study.

=head1 SYNOPSIS

    # Study
		$study = Bio::EnsEMBL::Variation::Study->new
       (-name => 'EGAS00000000001',
        -external_reference   => 'pubmed/17554300',
				-url => 'http://www.ebi.ac.uk/ega/page.php?page=study&study=EGAS00000000001&cat=www.wtccc.studies.xml.ega&subcat=BD'
				);
    ...


=head1 DESCRIPTION

This is a class representing a study from the ensembl-variation database.

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
    name of the study 
  Arg [-DESCRIPTION] :
    study description
	Arg [-URL] :
    string - url of the database/file where the data are stored
	Arg [-EXTERNAL_REFERENCE] :
    string - the pubmed/ids or project/study names
	Arg [-TYPE] :
    string - type of the study (e.g. GWAS)
	Arg [-SOURCE] :
    string - name of the source
	Arg [-ASSOCIATE] :
    array ref - list of the study objects associated with the current study
			
  Example    :
		
    $study = Bio::EnsEMBL::Variation::Study->new
       (-name => 'EGAS00000000001',
        -external_reference   => 'pubmed/17554300',
				-url => 'http://www.ebi.ac.uk/ega/page.php?page=study&study=EGAS00000000001&cat=www.wtccc.studies.xml.ega&subcat=BD'
				);

  Description: Constructor. Instantiates a new Study object.
  Returntype : Bio::EnsEMBL::Variation::Study
  Exceptions : none
  Caller     : general
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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

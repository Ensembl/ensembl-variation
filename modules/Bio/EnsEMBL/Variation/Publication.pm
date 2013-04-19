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


=head1 NAME

Bio::EnsEMBL::Variation::Publication - Ensembl representation of a publication.

=head1 SYNOPSIS

    # Publication
		$publication = Bio::EnsEMBL::Variation::Publication->new
                                       (-title   => 'A brief history of Alzheimer\'s disease gene discovery',
                                        -authors => 'Tanzi RE',
		 		        -pmid    => 23042215
				);
    ...


=head1 DESCRIPTION

This is a class representing a publication from the ensembl-variation database.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Publication;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);

our @ISA = ('Bio::EnsEMBL::Storable');


=head2 new

  Arg [-dbID] :
    see superclass constructor
  Arg [-ADAPTOR] :
    see superclass constructor
  Arg [-TITLE] :
    title of the publication 
  Arg [-AUTHORS] :
    publication authors
  Arg [-PMID] :
    string - the PubMed id for this publication
  Arg [-PMCID :
    string - the PubMed Central id for this publication
			
  Example    : 	$publication = Bio::EnsEMBL::Variation::Publication->new
                                       (-title   => 'A brief history of Alzheimer's disease gene discovery',
                                        -authors => 'Tanzi RE',
		 		        -pmid    => 23042215
				);
		    

  Description: Constructor. Instantiates a new Publication object.
  Returntype : Bio::EnsEMBL::Variation::Publication
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
	my ($dbID, $adaptor, $title, $authors, $pmid, $pmcid, $variants ) = 
			rearrange([qw(dbID ADAPTOR TITLE AUTHORS PMID PMCID VARIANTS)], @_);

  $self = {
      'dbID'     => $dbID,
      'adaptor'  => $adaptor,
      'title'    => $title,
      'authors'  => $authors,
      'pmid'     => $pmid,
      'pmcid'    => $pmcid, 
      'variants' => $variants     
  };
	
  return bless $self, $class;
}


=head2 title

  Arg [1]    : string $newval (optional)
               The new value to set the title attribute to
  Example    : $title = $obj->title()
  Description: Getter/Setter for the title attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub title{
  my $self = shift;
  return $self->{'title'} = shift if(@_);
  return $self->{'title'};
}


=head2 authors

  Arg [1]    : string $newval (optional)
               The new value to set the authors attribute to
  Example    : $name = $obj->authors()
  Description: Getter/Setter for the authors attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub authors{
  my $self = shift;
  return $self->{'authors'} = shift if(@_);
  return $self->{'authors'};
}


=head2 pmid

  Arg [1]    : string $newval (optional)
               The new value to set the pmid attribute to
  Example    : $name = $obj->pmid()
  Description: Getter/Setter for the PubMed ID attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub pmid{
  my $self = shift;
  return $self->{'pmid'} = shift if(@_);
  return $self->{'pmid'};
}


=head2 pmcid

  Arg [1]    : string $newval (optional)
               The new value to set the pmcid attribute to
  Example    : $name = $obj->pmcid()
  Description: Getter/Setter for the PubMed Central ID  attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub pmcid{
  my $self = shift;
  return $self->{'pmcid'} = shift if(@_);
  return $self->{'pmcid'};
}


1;

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
use Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor;

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
  Arg [-DOI :
    string - the doi for this publication
			
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
	my ($dbID, $adaptor, $title, $authors, $pmid, $pmcid, $year, $doi, $ucsc_id,$variants ) = 
			rearrange([qw(dbID ADAPTOR TITLE AUTHORS PMID PMCID YEAR DOI UCSC_ID VARIANTS)], @_);

  $self = {
      'dbID'     => $dbID,
      'adaptor'  => $adaptor,
      'title'    => $title,
      'authors'  => $authors,
      'pmid'     => $pmid,
      'pmcid'    => $pmcid, 
      'year'     => $year,
      'doi'      => $doi,
      'ucsc_id'  => $ucsc_id,
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
  Example    : $author_list = $obj->authors()
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
  Example    : $pmid = $obj->pmid()
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
  Example    : $pmcid = $obj->pmcid()
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

=head2 year

  Arg [1]    : string $newval (optional)
               The new value to set the pmcid attribute to
  Example    : $publication_year = $obj->year()
  Description: Getter/Setter for the the publication year attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub year{
  my $self = shift;
  return $self->{'year'} = shift if(@_);
  return $self->{'year'};
}

=head2 doi

  Arg [1]    : string $newval (optional)
               The new value to set the doi attribute to
  Example    : $doi = $obj->doi()
  Description: Getter/Setter for the the publication doi attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub doi{
  my $self = shift;
  return $self->{'doi'} = shift if(@_);
  return $self->{'doi'};
}


=head2 ucsc_id

  Arg [1]    : string $newval (optional)
               The new value to set the ucsc_id attribute to
  Example    : $ucsc_id = $obj->ucsc_id()
  Description: Getter/Setter for the the publication USCS external id attribute
               (enables link to UCSC website)
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub ucsc_id{
  my $self = shift;
  return $self->{'ucsc_id'} = shift if(@_);
  return $self->{'ucsc_id'};
}
=head2 variations

  Arg [1]    : array ref [optional]
               an array of variations cited in the publication
  Example    : $variant_objects = $obj->variations()
  Description: Getter/Setter for variatio_citations
  Returntype : arrayref
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub variations{
   my $self = shift;
  return $self->{'variants'} = shift if(@_);
 
   my $variation_adaptor = $self->adaptor->db->get_VariationAdaptor();      
   $self->{'variants'} = $variation_adaptor->fetch_all_by_Publication($self);

   return $self->{'variants'};
}

1;

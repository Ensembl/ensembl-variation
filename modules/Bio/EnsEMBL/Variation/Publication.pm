=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

# valid source names
# this must correspond the citation sources defined in attrib table
our %SOURCES = (
  'dbSNP'   => 1,
  'EPMC'    => 1,
  'UCSC'    => 1,
  'ClinVar' => 1,
  'dbGaP'   => 1,
  'GWAS'    => 1,
);

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
  Description: Getter/Setter for variation_citations
  Returntype : arrayref
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub variations{
   my $self = shift;
  return $self->{'variants'} = shift if(@_);
 
   my $variation_adaptor = $self->adaptor->db->get_VariationAdaptor();      
   $self->{'variants'} = $variation_adaptor->fetch_all_by_publication($self);

   return $self->{'variants'};
}

=head2 set_variation_id_to_source

  Arg [1]    : Integer $variation_id
               String $source
  Example    : $obj->set_variation_id_to_source(1234, 'dbSNP')
  Description: Setter for the publication sources
  Returntype : none
  Exceptions : throw if variation_id is not valid
               throw if source is not valid
  Caller     : general
  Status     : At Risk

=cut

sub set_variation_id_to_source{
  my $self         = shift;
  my $variation_id = shift;
  my $source       = shift;

  if($variation_id !~ /^\d+$/ || !defined($variation_id)) {
    throw('variation_id is not valid');
  }

  if(!defined($source) || !defined($SOURCES{$source})) {
    throw('Source is not valid');
  }

  $self->{'variation_id_to_source'}->{$variation_id}->{$source} = 1;
}

=head2 get_all_sources_by_Variation

  Arg [1]    : Variation object
  Example    : my $sources = $obj->get_all_sources_by_Variation($variation)
  Description: Get all publication sources for a given variation object
  Returntype : listref
               example: ['EPMC', 'dbSNP']
  Exceptions : throw if variation object is not valid
  Caller     : general
  Status     : At Risk

=cut

sub get_all_sources_by_Variation{
  my $self      = shift;
  my $variation = shift;

  if( !defined($variation) || (defined($variation) && ref($variation ne 'Variation')) ) {
    throw('Variation argument is required');
  }

  my $variation_id = $variation->dbID();
  if($self->{'variation_id_to_source'}->{$variation_id}) {
    my @sources = keys %{$self->{'variation_id_to_source'}->{$variation_id}};
    return \@sources;
  }
  else {
    return [];
  }
}

1;

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

# Ensembl module for Bio::EnsEMBL::Variation::Source
#
#


=head1 NAME

Bio::EnsEMBL::Variation::Source - Ensembl representation of a source.

=head1 SYNOPSIS

    # Source
		$source = Bio::EnsEMBL::Variation::Source->new
       (-name        => 'dbSNP',
        -version     => 138,
        -description => 'Variants (including SNPs and indels) imported from dbSNP',
				-url         => 'http://www.ncbi.nlm.nih.gov/projects/SNP/'
			 );
    ...


=head1 DESCRIPTION

This is a class representing a source from the ensembl-variation database.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Source;

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
    string - name of the source
  Arg [-VERSION] :
    int - source version
  Arg [-DESCRIPTION] :
    string - source description
	Arg [-URL] :
    string - url of the source websit
	Arg [-TYPE] :
    string - type of the source (e.g. chip, lsdb)
	Arg [-SOMATIC_STATUS] :
    string - information whether the data source is germline, somatic or mixed (both)
	Arg [-DATA_TYPES] :
    array ref - list of the source data type in the variation database (e.g. variation, variation_synonym, structural_variation, phenotype_feature, study)
			
  Example    :
		
    $source = Bio::EnsEMBL::Variation::Source->new
       (-name        => 'dbSNP',
        -version     => 138,
        -description => 'Variants (including SNPs and indels) imported from dbSNP',
				-url         => 'http://www.ncbi.nlm.nih.gov/projects/SNP/'
			 );

  Description: Constructor. Instantiates a new Source object.
  Returntype : Bio::EnsEMBL::Variation::Source
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
	my ($dbID,$adaptor,$source_name,$source_version,$source_description,$source_url,$source_type,$source_somatic_status,$source_data_types) = 
			rearrange([qw(dbID ADAPTOR NAME VERSION DESCRIPTION URL TYPE SOMATIC_STATUS DATA_TYPES)], @_);

  $self = {
			'dbID'           => $dbID,
			'adaptor'        => $adaptor,
  		'name'           => $source_name,
  		'version'        => $source_version,
			'description'    => $source_description,
			'url'            => $source_url,
			'type'           => $source_type,
			'somatic_status' => $source_somatic_status,
			'data_types'     => $source_data_types
	};
	
	return bless $self, $class;
}

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
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


=head2 version

  Arg [1]    : string $newval (optional)
               The new value to set the version attribute to
  Example    : $version = $obj->version()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub version{
  my $self = shift;
  return $self->{'version'} = shift if(@_);
  return $self->{'version'};
}


=head2 formatted_version

  Arg [1]    : string $newval (optional)
               The new value to set the version attribute to
  Example    : $formatted_version = $obj->formatted_version()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub formatted_version{
  my $self = shift;
  $self->{'version'} = shift if(@_);
  
  my $version = $self->{'version'};
  my $year_prefix = 20;
  
  if ($version =~ /^($year_prefix\d{2})(\d{2})(\d{2})$/) { # Year+Month+Day
    $version = "$3/$2/$1"; # Day/Month/Year
  }
  elsif ($version =~ /^($year_prefix\d{2})(\d{2})$/) { # Year+Month
    $version = "$2/$1"; # Month/Year
  }
  elsif ($version =~ /^($year_prefix\d{2})(\d)$/) { # HGMD (Year+Version)
    $version = "$1.$2"; # Year.Version
  }
  
  return $version;
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


=head2 somatic_status

  Arg [1]    : string $newval (optional)
               The new value to set the somatic status attribute to
  Example    : $name = $obj->somatic_status()
  Description: Getter/Setter for the somatic status attribute.
               e.g. germline, somatic, mixed (both)
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub somatic_status{
  my $self = shift;
  return $self->{'somatic_status'} = shift if(@_);
  return $self->{'somatic_status'};
}


=head2 get_all_data_types

  Arg [1]    : none
  Example    : my @data_types = @{$s->get_all_data_types()};
  Description: Retrieves all data types for this source. Current
               possible data types are variation, variation_synonym, 
               structural_variation, phenotype_feature, study.
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_data_types {
    my $self = shift;
    return $self->{'data_types'};

}

1;

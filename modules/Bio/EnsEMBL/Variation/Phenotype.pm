=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

# Ensembl module for Bio::EnsEMBL::Variation::Phenotype
#
#


=head1 NAME

Bio::EnsEMBL::Variation::Phenotype - Ensembl representation of a phenotype.

=head1 SYNOPSIS

    my $phenotype = Bio::EnsEMBL::Variation::Phenotype->new(-NAME => 'Type I diabetes');

=head1 DESCRIPTION

This is a class representing a phenotype from the ensembl-variation database.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Phenotype;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);

our @ISA = ('Bio::EnsEMBL::Storable');


=head2 new
    Arg [-DESCRIPTION] :
    phenotype description
			
  Example    :
		
    $phenotype = Bio::EnsEMBL::Variation::Phenotype->new(-DESCRIPTION => 'Hemostatic factors and hematological phenotypes');

  Description: Constructor. Instantiates a new Phenotype object.
  Returntype : Bio::EnsEMBL::Variation::Phenotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
    my $caller = shift;
    my $class  = ref($caller) || $caller;
    my $self = $class->SUPER::new(@_);
    my ($dbID, $description, $name) = rearrange([qw(dbID DESCRIPTION NAME)], @_);
    $self = {
        'dbID'        => $dbID,
        'description' => $description,
        'name'        => $name,
    };
    return bless $self, $class;
}

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}

=head2 dbID

  Example    : $name = $obj->dbID()
  Description: Getter/Setter for the dbIDattribute
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub dbID {
    my $self = shift;
    return $self->{'dbID'} = shift if(@_);
    return $self->{'dbID'};
}

=head2 name

  Example    : $name = $obj->name()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
    my $self = shift;
    return $self->{'name'} = shift if(@_);
    return $self->{'name'};
}


=head2 description

  Example    : $name = $obj->description()
  Description: Getter/Setter for the description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
    my $self = shift;
    return $self->{'description'} = shift if(@_);
    return $self->{'description'};
}


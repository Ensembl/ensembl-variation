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

#
# EnsEMBL module for Bio::EnsEMBL::Variation::DBSQL::MetaContainer
#
# Cared for by Daniel Rios
#
# POD documentation - main docs before the code

=head1 NAME

  Bio::EnsEMBL::Variation::DBSQL::MetaContainer - 
  Encapsulates all access to variation database meta information

=head1 SYNOPSIS

  my $meta_container = $db_adaptor->get_MetaContainer();

  my $default_population = $meta_container->get_default_LDPopulation();

=head1 DESCRIPTION

  An object that encapsulates specific access to variation db meta data

=head1 METHODS

=cut

package Bio::EnsEMBL::Variation::DBSQL::MetaContainer;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseMetaContainer;
use Bio::EnsEMBL::Utils::Exception;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseMetaContainer);


sub get_schema_version {
  my $self = shift;

  my $arrRef = $self->list_value_by_key( 'schema_version' );

  if( @$arrRef ) {
    my ($ver) = ($arrRef->[0] =~ /^\s*(\d+)\s*$/);
    if(!defined($ver)){ # old style format
      return 0;
    }
    return $ver;
  } else {
    warning("Please insert meta_key 'schema_version' " .
         "in meta table at variation db.\n");
  }
  return 0;
}

sub ploidy {
  my $self = shift;
  
  my $values = $self->list_value_by_key('ploidy');
  
  # default to 2
  return scalar @$values ? $values->[0] : 2;
}

1;

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

# Ensembl module for Bio::EnsEMBL::Variation::GenotypeCode
#
#


=head1 NAME

Bio::EnsEMBL::Variation::GenotypeCode - Class for a genotype code object

=head1 SYNOPSIS

    print (join "|", $genotypecode->genotype), "\n";

=head1 DESCRIPTION

This is class representing a GenotypeCode. It is intended for internal use only.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::GenotypeCode;

use Bio::EnsEMBL::Storable;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}

=head2 genotype

  Example    : $gt = $genotypecode->genotype()
  Description: Getter for the genotype represented by this GenotypeCode object.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub genotype {
  return $_[0]->{genotype};
}

=head2 phased

  Example    : $p = $genotypecode->phased()
  Description: Getter for the phased status of this GenotypeCode object.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub phased {
  return $_[0]->{phased};
}

1;

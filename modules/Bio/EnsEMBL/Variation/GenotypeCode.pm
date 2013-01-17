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

# Ensembl module for Bio::EnsEMBL::Variation::GenotypeCode
#
# Copyright (c) 2004 Ensembl
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

1;

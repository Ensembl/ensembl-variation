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

# Ensembl module for Bio::EnsEMBL::Variation::Phenotype
#
#


=head1 NAME

Bio::EnsEMBL::Variation::QTL - Ensembl representation of a quantitative trait
locus (QTL).

=head1 SYNOPSIS

    my $qtl = Bio::EnsEMBL::Variation::QTL->new(-stable_id => 'QTL00123');

=head1 DESCRIPTION

This is a class representing a QTL.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::QTL;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);

our @ISA = ('Bio::EnsEMBL::Storable');


=head2 new
  Arg [-stable_id] : stable identifier or name for QTL
  Example : my $qtl = Bio::EnsEMBL::Variation::QTL->new(-stable_id => 'QTL00123');

  Description: Constructor. Instantiates a new QTL object.
  Returntype : Bio::EnsEMBL::Variation::QTL
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
    my $caller = shift;
    my $class  = ref($caller) || $caller;
    my $self = $class->SUPER::new(@_);
    my ($stable_id) = rearrange([qw(STABLE_ID)], @_);
    $self = {
        'stable_id' => $stable_id,
    };
    return bless $self, $class;
}

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}

sub stable_id {
    my $self = shift;
    return $self->{stable_id} = shift if @_;
    return $self->{stable_id};
}

1;
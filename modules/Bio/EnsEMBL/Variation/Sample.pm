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

=head1 NAME
Bio::EnsEMBL::Variation::Sample
=head1 SYNOPSIS
=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Variation::Sample;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Storable;

our @ISA = ('Bio::EnsEMBL::Storable');

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my (
    $dbID,
    $adaptor,
    $individual_id,
    $desc,
    $study_id,
    $display_flag,
    $has_coverage) = rearrange([qw(dbID adaptor individual_id description study_id display has_coverage)], @_);

  $display_flag ||= 'UNDISPLAYABLE'; 

  return bless {
    'dbID'    => $dbID,
    'adaptor' => $adaptor,
    'individual_id' => $individual_id,
    'description' => $desc,
    'study_id' => $study_id,
    'display' => $display_flag,
    'has_coverage' => $has_coverage,
  }, $class;
}


sub name {
  my $self = shift;
  return $self->{'name'} = shift if (@_);
  return $self->{'name'};
}

sub description {
  my ($self, $description) = @_;
  $self->{description} = $description if defined $description;
  return $self->{'description'};
}

sub display {
  my $self = shift;
  if (@_) {
    my $display = shift;
    $display = uc($display);
    unless (grep $_ eq $display, ('UNDISPLAYABLE', 'REFERENCE', 'DISPLAYABLE', 'DEFAULT')) {
      throw('Display flag must be one of "REFERENCE", "DEFAULT", "DISPLAYABLE", "UNDISPLAYABLE"');
    }
    $self->{'display'} = $display;
  }
  return $self->{'display'};
}

sub has_coverage {
  my $self = shift;
  if (@_){
    my $has_coverage = shift;
    return $self->{'has_coverage'} = $has_coverage;
  }
  return $self->{'has_coverage'} || 0;
}

sub individual {
  my $self = shift;
}

sub get_all_Populations {
  my $self = shift;
  if (!defined($self->{populations})) {
    if (defined ($self->{'adaptor'})) {
      my $pop_adaptor = $self->{'adaptor'}->db()->get_PopulationAdaptor();
      $self->{populations} = $pop_adaptor->fetch_all_by_Sample($self);
    }
    $self->{populations} ||= [];
  }
  return $self->{populations};
}

1;

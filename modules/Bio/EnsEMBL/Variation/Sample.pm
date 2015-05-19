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
    $individual,
    $individual_id,
    $name,
    $desc,
    $study,
    $study_id,
    $display_flag,
    $has_coverage) = rearrange([qw(dbID adaptor individual individual_id name description study study_id display has_coverage)], @_);

  $display_flag ||= 'UNDISPLAYABLE'; 

  return bless {
    'dbID'    => $dbID,
    'adaptor' => $adaptor,
    'individual' => $individual,
    'individual_id' => $individual_id,
    'name' => $name,
    'description' => $desc,
    'study' => $study,
    'study_id' => $study_id,
    'display' => $display_flag,
    'has_coverage' => $has_coverage,
  }, $class;
}

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
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
  if (@_) {
    if (!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Individual')) {
      throw("Bio::EnsEMBL::Variation::Individual argument expected");
    }
    $self->{'individual'} = shift;
  }
  elsif (!defined($self->{'individual'}) && $self->adaptor() && defined($self->{'individual_id'})) {
    # lazy-load from database on demand
    my $ia = $self->adaptor->db()->get_IndividualAdaptor();
    $self->{'individual'} = $ia->fetch_by_dbID($self->{'individual_id'});
  }

  return $self->{'individual'};
}

sub study {
  my $self = shift;
  if (@_) {
    if (!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Study')) {
      throw("Bio::EnsEMBL::Variation::Study argument expected");
    }
    $self->{'study'} = shift;
  }
  elsif (!defined($self->{'study'}) && $self->adaptor() && defined($self->{'study_id'})) {
    # lazy-load from database on demand
    my $sa = $self->adaptor->db()->get_StudyAdaptor();
    $self->{'study'} = $sa->fetch_by_dbID($self->{'study_id'});
  }

  return $self->{'study'};
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

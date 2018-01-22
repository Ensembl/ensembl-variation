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

=head1 NAME
Bio::EnsEMBL::Variation::Sample
=head1 SYNOPSIS

  my $study_1000G = Bio::EnsEMBL::Variation::Study->new(
    -name => '1000G phase 3',
    -description => 'Whole-genome sequencing',
  );

  my $individual = Bio::EnsEMBL::Variation::Individual(
    -name => 'NA18967',
  );
  $individual = $individual_adaptor->store($individual);

  my $sample_1000G = Bio::EnsEMBL::Variation::Sample->new(
    -individual => $individual,
    -study => $study,
  );

=head1 DESCRIPTION
This is a class representing a single Sample. A Sample belongs to one Individual. A
Sample may be part of one population or several. 
=cut

package Bio::EnsEMBL::Variation::Sample;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Storable;

our @ISA = ('Bio::EnsEMBL::Storable');

=head2 new
  Arg [-dbID]                : int - unique internal identifier
  Arg [-ADAPTOR]             : Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor
  Arg [-INDIVIDUAL]          : Bio::EnsEMBL::Variation::Individual 
                               - Individual to which this sample belongs to
  Arg [-INDIVIDUAL_ID]       : int - set the internal id for the individual object.
                               Individual object can be retrieved on demand.
  Arg [-NAME]                : string - name of this sample
  Arg [-DESCRIPTION]         : string description - description of this sample
  Arg [-STUDY]               : Bio::EnsEMBL::Variation::Study
                               - study object storing additional information e.g.
                                 on the study in which genotypes for this sample
                                 have been computed
  Arg [-STUDY_ID]            : int - set the internal id for the study object,

  Example                    : $sample = Bio::EnsEMBL::Variation::Sample->new(
                                          -name => 'WI530.07',
                                          -description => 'african',
                                          -gender => 'Male',
                                          -father_individual => $father_ind,
                                          -mother_individual => $mother_ind,
                                          -type_individual => 'outbred',
                                          -type_description => 'a single organism which breeds freely'
                                         );
  Description                 : Constructor instantiates a Sample object.
  Returntype                  : Bio::EnsEMBL::Variation::Sample
  Exceptions                  : Throw if 
  Caller                      : general
  Status                      : Stable
=cut
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
    $synonyms,
    $display_flag,
    $has_coverage) = rearrange([qw(dbID adaptor individual individual_id name description study study_id synonyms display has_coverage)], @_);

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
    'synonyms' => $synonyms || {},
    'display' => $display_flag,
    'has_coverage' => $has_coverage,
  }, $class;
}

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}

=head2 name
  Arg [1]    : String $name (optional)
               The new value to set the name attribute to
  Example    : $name = $sample->name()
  Description: Getter/Setter for the name attribute
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable
=cut

sub name {
  my $self = shift;
  return $self->{'name'} = shift if (@_);
  return $self->{'name'};
}

# get_all_synonyms
sub get_all_synonyms {
  my $self = shift;
  my $source = shift;

  if ($source) {
    $source = [$source];
  } else {
    $source = $self->get_all_synonym_sources();
  }

  my @synonyms;
  map { push ( @synonyms, keys ( %{$self->{synonyms}{$_} } ) ) } @{$source};

  return \@synonyms;
}

# get_all_synonym_sources
sub get_all_synonym_sources {
  my $self = shift;
  my @sources = keys %{$self->{'synonyms'}};
  return \@sources;
}

# add_synonym

sub add_synonym {
  my $self = shift;
  my $source = shift;
  my $synonym = shift;

  throw("source argument is required") if (!$source);
  throw("synonym argument is required") if (!$synonym);

  $self->{'synonyms'}{$source}{$synonym}++;

  return;
}


=head2 description
  Arg [1]    : String $description (optional)
               The new value to set the description attribute to
  Example    : $description = $individual->description()
  Description: Getter/Setter for the description attribute
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable
=cut

sub description {
  my ($self, $description) = @_;
  $self->{description} = $description if defined $description;
  return $self->{'description'};
}

=head2 display
  Arg [1]    : String $display (optional)
  Example    : $display = $display->display()
  Description: Getter/Setter for the display attribute which can be one of UNDISPLAYABLE,
               REFERENCE, DISPLAYABLE or DEFAULT.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable
=cut

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

=head2 has_coverage
  Arg [1]     : int $has_coverage (optional)
  Example     : $has_coverage = $sample->has_coverage();
  Description : Getter/Setter for the flag indicating if this sample has
                read coverage data available
  Returntype  : int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
=cut

sub has_coverage {
  my $self = shift;
  if (@_){
    my $has_coverage = shift;
    return $self->{'has_coverage'} = $has_coverage;
  }
  return $self->{'has_coverage'} || 0;
}

=head2 individual
  Arg [1]    : Bio::EnsEMBL::Variation::Individual (optional)
               The new value to set the individual attribute to
  Example    : $individual = $sample->individual()
  Description: Getter/Setter for Individual object. If this
               has not been set manually and this Sample has an attached
               adaptor, an attempt will be made to lazy-load it from the
               database using the individual_id.
  Returntype : Bio::EnsEMBL::Variation::Individual
  Exceptions : none
  Caller     : general
  Status     : Stable
=cut

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

=head2 study
  Arg [1]    : Bio::EnsEMBL::Variation::Study (optional)
               The new value to set the study attribute to
  Example    : $study = $sample->study()
  Description: Getter/Setter for Study object. If this
               has not been set manually and this Sample has an attached
               adaptor, an attempt will be made to lazy-load it from the
               database using the study_id if set.
  Returntype : Bio::EnsEMBL::Variation::Study
  Exceptions : none
  Caller     : general
  Status     : Stable
=cut

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

=head2 get_all_Populations
   Args        : none
   Example     : $populations = $sample->get_all_Populations();
   Description : Get all populations for this sample. Returns
                 empty list if there are none.
   ReturnType  : reference to list of Bio::EnsEMBL::Population objetcs
   Exceptions  : none
   Caller      : general
   Status      : Stable
=cut

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

=head2 add_Population
   Args        : Bio::EnsEMBL::Variation::Population
   Example     : $sample->add_Population($population);
   Description : Define populations to which the sample belongs
   ReturnType  : None
   Exceptions  : Throw on wrong argument
   Caller      : General
   Status      : Stable
=cut

sub add_Population {
  my $self = shift;
  my $population = shift;
  if (!$population || !$population->isa('Bio::EnsEMBL::Variation::Population')) {
    throw("Bio::EnsEMBL::Variation::Population argument expected");
  }
  push @{$self->{populations}}, $population;
  return;
}


1;

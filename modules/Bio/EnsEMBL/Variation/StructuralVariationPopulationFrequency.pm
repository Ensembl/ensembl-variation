=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

# Ensembl module for Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency
#
#


=head1 NAME

Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency - A module to represent the population frequency 
of a structural variant

=head1 SYNOPSIS

    # Structural variation population frequency representing a CNV
    $svpf = Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency->new
       (-population         => $pop, # Bio::EnsEMBL::Population object
        -structural_variant => $sv,  # Bio::EnsEMBL::StructuralVariation object
        -name               => "1000GENOMES:phase_3:GBR",
        -description        => "British in England and Scotland",
        -size               => 91,
        -samples_class      => { 'copy_number_gain' => { '1000GENOMES:phase_3:HG00253' => 'homozygous', 
                                                         '1000GENOMES:phase_3:HG00106' => 'heterozygous'}
                               }
                          
       );

    ...

    print $svpf->name(), "-", $svpf->size(), '(', $svpf->frequency(), ')', "\n";
    
=head1 DESCRIPTION

This is a class representing the frequency of a given structural variant in a population


=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(check_ref assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);

=head2 new

  Arg [-ADAPTOR]                  : Bio::EnsEMBL::Variation::DBSQL::StructuralVariationPopulationFrequencyAdaptor
  Arg [-NAME]                     : string - name of the population
  Arg [-DESCRIPTION]              : string - description of the population
  Arg [-SIZE]                     : int - the size of the population
  Arg [-SAMPLES_CLASS]            : hashref of sample IDs by SO terms. This also contains the zygosity of each sample
  Arg [-_POPULATION_ID]           : int - the internal id of the population object associated with this identifier.
                                    This may be provided instead of a population object so that the population may be 
                                    lazy-loaded from the database on demand if we need more information about the population.
  Arg [-_STRUCTURAL_VARIATION_ID] : int - the internal id of the structural variation object associated with this
                                    identifier. This may be provided instead of a structural variation object so that
                                    the structural variation may be lazy-loaded from the database on demand.
 
  Example               : $svpf = Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency->new
                          ( -population         => $pop, # Bio::EnsEMBL::Population object
                            -structural_variant => $sv,  # Bio::EnsEMBL::StructuralVariation object
                            -name               => "1000GENOMES:phase_3:GBR",
                            -description        => "British in England and Scotland",
                            -size               => 91,
                            -samples_class      => { 'copy_number_gain' => { '1000GENOMES:phase_3:HG00253' => 'homozygous', 
                                                                             '1000GENOMES:phase_3:HG00106' => 'heterozygous'}
                                                   }
                          );
       );
  Description           : Constructor. Instantiates a new StructuralVariationPopulationFrequency object
  Returntype            : Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency
  Exceptions            : none
  Caller                : general
  Status                : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($adaptor, $structural_variation_id, $population_id, $name, $desc, $size, $samples_class) = rearrange([
    'ADAPTOR', '_STRUCTURAL_VARIATION_ID', '_POPULATION_ID', 'NAME', 'DESCRIPTION', 'SIZE', 'SAMPLES_CLASS'], @_);

  return bless {
    'adaptor'                  => $adaptor,
    '_structural_variation_id' => $structural_variation_id,
    '_population_id'           => $population_id,
    'name'                     => $name,
    'description'              => $desc,
    'size'                     => $size,
    'samples_class'            => $samples_class,
  }, $class;
}


=head2 population

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Population $pop
  Example    : $pop = $svpf->population();
  Description: Getter/Setter for the population associated with this structural variant.
               If not set, and this StructuralVariationPopulationFrequency has an associated adaptor
               an attempt will be made to lazy-load the population from the database.
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub population {
  my ($self, $population) = @_;
  
  # set
  if (defined($population)) {
    
    assert_ref($population, 'Bio::EnsEMBL::Variation::Population');
 
    $self->{population} = $population;
  }
  
  # get
  if(!defined($self->{population}) && defined($self->{_population_id})) {
    my $pa = $self->adaptor->db->get_PopulationAdaptor();
    
    $self->{population} = $pa->fetch_by_dbID($self->{_population_id});
  }
  
  return $self->{population};
}


=head2 name

  Arg [1]    : String $name (optional)
               The new value to set the name attribute to
  Example    : $name = $svpf->name()
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

=head2 description

  Arg [1]    : String $description (optional) 
               The new value to set the description attribute to
  Example    : $description = $svpf->description()
  Description: Getter/Setter for the description attribute
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub description {
  my $self = shift;
  return $self->{'description'} = shift if (@_);
  return $self->{'description'};
}

=head2 size

  Arg [1]    : int $size (optional) 
               The new value to set the population size attribute to
  Example    : $size = $svpf->size()
  Description: Getter/Setter for the population size attribute
  Returntype : Int. Returns undef if information on size is not given.
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub size {
  my $self = shift;
  return $self->{'size'} = shift if (@_);
  return $self->{'size'};
}


=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::Variation::DBSQL::StructuralVariationPopulationFrequencyAdaptor $adaptor (optional)
               Set the adaptor for this StructuralVariationPopulationFrequency
  Example    : my $adaptor = $svpf->adaptor()
  Description: Getter/Setter for the adaptor.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::StructuralVariationPopulationFrequencyAdaptor
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub adaptor {
  my $self = shift;
  $self->{adaptor} = shift if @_;
  return $self->{adaptor};
}


=head2 frequency

  Arg [1]    : float $frequency (optional)
               The new value to set the SV frequency attribute to
  Example    : $frequency = $svpf->frequency()
  Description: Getter/Setter for the SV frequency attribute
  Returntype : float
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub frequency {
  my $self = shift;
  if (@_) {
    $self->{frequency} = shift @_;
  }
  elsif (scalar(keys(%{$self->{samples_class}})) && $self->{size}) {
    my $samples_count = 0;
    # Loop over SO terms
    foreach my $SO_term (keys(%{$self->{samples_class}})) {
      # Loop over Sample IDs
      foreach my $sample_id (keys(%{$self->{samples_class}{$SO_term}})) {
        $samples_count += ($self->{samples_class}{$SO_term}{$sample_id} eq 'homozygous') ? 2 : 1;
      }
    }
    $self->{frequency} = $samples_count / ($self->{size}*2); 
  }
  
  return $self->{frequency};
}


=head2 frequencies_by_class_SO_term

  Arg [1]    : String $class_SO_term (optional)
               The new value to get the frequency by with the SO term attribute
  Example    : $SO_term_frequencies = $svpf->frequencies_by_class_SO_term($class_SO_term)
  Description: Getter for the frequencies by SO term
               If a SO term is given as parameter, it returns an hashref with the corresponding frequency for the population, e.g.
               $SO_term_frequencies = $svpf->frequencies_by_class_SO_term('copy_number_gain');
               $SO_term_frequencies will contain {'copy_number_gain' => 0.02}
               If no SO terms are given as parameter, it returns an hashref with all the associated SO terms - frequencies for the population, e.g.
               $SO_term_frequencies = $svpf->frequencies_by_class_SO_term();
               $SO_term_frequencies will contain {'copy_number_gain' => 0.02, 'copy_number_loss' => 0.05}
  Returntype : hashref
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub frequencies_by_class_SO_term {
  my $self = shift;
  my $class_SO_term = shift;
  
  my %class_SO_term_freq;
  
  if ($class_SO_term) {
    if ($self->{samples_class}->{$class_SO_term}) {
      my $samples_count = 0;
      # Loop over Sample IDs
      foreach my $sample_id (keys(%{$self->{samples_class}->{$class_SO_term}})) {
        $samples_count += ($self->{samples_class}{$class_SO_term}{$sample_id} eq 'homozygous') ? 2 : 1;
      }
      $class_SO_term_freq{$class_SO_term} = $samples_count / ($self->{size}*2); 
    }
    else {
      $class_SO_term_freq{$class_SO_term} = "No data found";
    }
  }
  elsif (scalar(keys(%{$self->{samples_class}})) && $self->{size}) {
    # Loop over SO terms
    foreach my $SO_term (keys(%{$self->{samples_class}})) {
      my $samples_SO_count = 0;
      # Loop over Sample IDs
      foreach my $sample_id (keys(%{$self->{samples_class}->{$SO_term}})) {
        $samples_SO_count += ($self->{samples_class}{$SO_term}{$sample_id} eq 'homozygous') ? 2 : 1;
      }
      # Frequency by SO term
      $class_SO_term_freq{$SO_term} = $samples_SO_count / ($self->{size}*2);
    }
  }
  
  return \%class_SO_term_freq;
}


=head2 get_all_Samples

  Arg [1]    : none
  Example    : @samples = @{$p->get_all_Samples()};
  Description: Retrieves all the Sample objects belonging to this Population and having information for the StructuralVariation.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Sample objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Samples {
  my $self = shift;
  my $sa = $self->adaptor->db->get_SampleAdaptor;
  return (defined $sa ? $sa->fetch_all_by_dbID_list($self->{samples_class}) : []);
}


=head2 structural_variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::StructuralVariation or 
               Bio::EnsEMBL::Variation::SupportingStructuralVariation $structural_variation
  Example    : $sv = $svpf->structural_variation();
  Description: Getter/Setter for the structural variant associated with this feature.
               If not set, and this StructuralVariationFeature has an associated adaptor
               an attempt will be made to lazy-load the structural variation from the
               database.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation or 
               Bio::EnsEMBL::Variation::SupportingStructuralVariation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub structural_variation {
  my $self = shift;

  if(@_) {
  
    unless (check_ref($_[0], 'Bio::EnsEMBL::Variation::StructuralVariation') || check_ref($_[0], 'Bio::EnsEMBL::Variation::SupportingStructuralVariation')) {
      throw("Bio::EnsEMBL::Variation::StructuralVariation or Bio::EnsEMBL::Variation::SupportingStructuralVariation argument expected");
    }
    $self->{'structural_variation'} = shift;
  }
  elsif(!defined($self->{'structural_variation'}) && $self->{'adaptor'} && defined($self->{'_structural_variation_id'})) {
    # lazy-load from database on demand
    my $sva = $self->{'adaptor'}->db()->get_StructuralVariationAdaptor();
    $self->{'structural_variation'} = $sva->fetch_by_dbID($self->{'_structural_variation_id'});
    if (!defined($self->{'structural_variation'})) {
      $sva = $self->{'adaptor'}->db()->get_SupportingStructuralVariationAdaptor();
      $self->{'structural_variation'} = $sva->fetch_by_dbID($self->{'_structural_variation_id'});
    }
  }

  return $self->{'structural_variation'};
}

1;

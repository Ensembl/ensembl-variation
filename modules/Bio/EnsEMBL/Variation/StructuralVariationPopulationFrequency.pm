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
        -region_name        => 'X',
        -samples_class      => { 'copy_number_gain' => { '1000GENOMES:phase_3:HG00253' => 2, 
                                                         '1000GENOMES:phase_3:HG00106' => 1 }
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
  Arg [-REGION_NAME]              : string - region/chromosome where the structural variant is located
  Arg [-SAMPLES_CLASS]            : hashref of sample IDs by SO terms. This also contains the numeric zygosity of each sample
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
                            -region_name        => 'X',
                            -samples_class      => { 'copy_number_gain' => { '1000GENOMES:phase_3:HG00253' => 2, 
                                                                             '1000GENOMES:phase_3:HG00106' => 1 }
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

  my ($adaptor, $structural_variation_id, $population_id, $name, $desc, $size, $samples_class, $region_name) = rearrange([
    'ADAPTOR', '_STRUCTURAL_VARIATION_ID', '_POPULATION_ID', 'NAME', 'DESCRIPTION', 'SIZE', 'SAMPLES_CLASS', 'REGION_NAME'], @_);

  return bless {
    'adaptor'                  => $adaptor,
    '_structural_variation_id' => $structural_variation_id,
    '_population_id'           => $population_id,
    'name'                     => $name,
    'description'              => $desc,
    'size'                     => $size,
    'samples_class'            => $samples_class,
    'region_name'              => $region_name
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


=head2 allele_count

  Arg [1]    : int $count (optional)
               The new value to set the SV allele count attribute to
  Example    : $count = $svpf->allele_count()
  Description: Getter/Setter for the SV allele count attribute
  Returntype : Int
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub allele_count {
  my $self = shift;
  if (@_) {
    $self->{allele_count} = shift @_;
  }
  elsif (scalar(keys(%{$self->{samples_class}}))) {
    my $allele_count_by_class_SO_term = $self->allele_count_by_class_SO_term();
    my $allele_count = 0;
    # Loop over allele count by SO terms
    foreach my $SO_term_allele_count (values(%{$allele_count_by_class_SO_term})) {
      $allele_count += $SO_term_allele_count;
    }
    $self->{allele_count} = $allele_count; 
  }
  
  return $self->{allele_count};
}


=head2 allele_count_by_class_SO_term

  Arg [1]    : String $class_SO_term (optional)
               The new value to get the allele count by with the SO term attribute
  Example    : $SO_term_frequencies = $svpf->allele_count_by_class_SO_term($class_SO_term)
  Description: Getter for the allele_count by SO term
               If a SO term is given as parameter, it returns an hashref with the corresponding allele count for the population, e.g.
               $SO_term_allele_count = $svpf->allele_count_by_class_SO_term('copy_number_gain');
               $SO_term_allele_count will contain {'copy_number_gain' => 3}
               If no SO terms are given as parameter, it returns an hashref with all the associated SO terms - allele count for the population, e.g.
               $SO_term_allele_count = $svpf->allele_count_by_class_SO_term();
               $SO_term_allele_count will contain {'copy_number_gain' => 3, 'copy_number_loss' => 7}
  Returntype : hashref
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub allele_count_by_class_SO_term {
  my $self = shift;
  my $class_SO_term = shift;
  
  my %class_SO_term_allele_count;
  
  # Defined SO term
  if ($class_SO_term) {
    if ($self->{samples_class}->{$class_SO_term}) {
      
      my $allele_count = 0;
      my @sample_ids = keys(%{$self->{samples_class}{$class_SO_term}});
      
      # Get sample objects if the SV falls into a special chromosome
      my $samples_list = ($self->{region_name} =~ /^(Y|X)$/) ? $self->get_Samples_by_dbID_list(\@sample_ids) : {};
      
      foreach my $sample_id (@sample_ids) {
        
        my $sample = ($samples_list) ? $samples_list->{$sample_id} : undef;
          
        if ($self->{region_name} eq 'X') {
          $allele_count += ($sample->individual->gender eq 'Male') ? 1 : $self->{samples_class}{$class_SO_term}{$sample_id};
        }
        elsif ($self->{region_name} eq 'Y') {
          $allele_count += 1 if ($sample->individual->gender eq 'Male');
        }
        else {
          $allele_count += $self->{samples_class}{$class_SO_term}{$sample_id};
        }
      }
    }
    else {
      $class_SO_term_allele_count{$class_SO_term} = "No data found";
    }
  }
  # Not defined SO term: returns the allele counts for each SO term available
  elsif (scalar(keys(%{$self->{samples_class}}))) {
    # Loop over SO terms
    foreach my $SO_term (keys(%{$self->{samples_class}})) {
    
      my $allele_SO_count = 0;
      my @sample_ids = keys(%{$self->{samples_class}->{$SO_term}});
      
      # Get sample objects if the SV falls into a special chromosome
      my $samples_list = ($self->{region_name} =~ /^(Y|X)$/) ? $self->get_Samples_by_dbID_list(\@sample_ids) : {};
      
      # Loop over Sample IDs
      foreach my $sample_id (@sample_ids) {
          
        my $sample = ($samples_list) ? $samples_list->{$sample_id} : undef;
        
        if ($self->{region_name} eq 'X') {
          $allele_SO_count += ($sample->individual->gender eq 'Male') ? 1 : $self->{samples_class}{$SO_term}{$sample_id};
        }
        elsif ($self->{region_name} eq 'Y') {
          $allele_SO_count += 1 if ($sample->individual->gender eq 'Male');
        }
        else {
          $allele_SO_count += $self->{samples_class}{$SO_term}{$sample_id};
        }
      }

      # Allele count by SO term
      $class_SO_term_allele_count{$SO_term} = $allele_SO_count;
    }
  }
  
  return \%class_SO_term_allele_count;
}


=head2 get_individuals_allele_count

  Arg [1]    : none
  Example    : $allele_count = $svpf->get_individuals_allele_count();
  Description: Retrieves all the total number of allele available in a population.
               Normally this corresponds to the number of individuals x ploidy (e.g. in human: 91 individuals * 2 = 182 alleles).
               However for human there are special cases with the chromosomes X and Y where the ploidy is not the same for male and female.
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_individuals_allele_count {
  my $self = shift;
  
  return $self->{individuals_allele_count} if ($self->{individuals_allele_count});
  
  # Specific to human
  if ($self->{region_name} =~ /^(Y|X)$/) {
    my $count = 0;
    my $samples = $self->population->get_all_Samples();
    foreach my $sample (@$samples) {
      if ($self->{region_name} eq 'X') {
        $count += ($sample->individual->gender eq 'Male') ? 1 : 2;
      }
      elsif ($self->{region_name} eq 'Y') {
        $count += 1 if ($sample->individual->gender eq 'Male');
      }
    }
     $self->{individuals_allele_count} = $count;
  }
  else {
     $self->{individuals_allele_count} = $self->{size}*2; # Might need to replace the "2" by the ploidie of the species 
  }
  return $self->{individuals_allele_count};
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
    my $allele_count = $self->allele_count();
    my $total_allele_count = $self->get_individuals_allele_count();

    return 0 if ($total_allele_count == 0 || !$total_allele_count);

    $self->{frequency} = $allele_count / $total_allele_count;
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
  
  my $total_allele_count = $self->get_individuals_allele_count();
  
  my $class_SO_term_count;
  if ($class_SO_term) {
    $class_SO_term_count = $self->allele_count_by_class_SO_term($class_SO_term);
  }
  else {
    $class_SO_term_count = $self->allele_count_by_class_SO_term();
  }
  foreach my $SO_term (keys(%$class_SO_term_count)) {
    if ($class_SO_term_count->{$SO_term} ne "No data found" && $total_allele_count) {
    
      my $count = $class_SO_term_count->{$SO_term};

      $class_SO_term_freq{$SO_term} = $count / $total_allele_count;
    }
  }
  
  return \%class_SO_term_freq;
}


=head2 get_Samples_by_dbID_list

  Arg [1]    : none
  Example    : @samples = @{$svpf->get_Samples_by_dbID_list()};
  Description: Retrieves all the Sample objects belonging to this Population and having information for the StructuralVariation.
  Returntype : hashref with sample ID as key and its corresponding Bio::EnsEMBL::Variation::Sample object as value
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_Samples_by_dbID_list {
  my $self = shift;
  my $sample_ids = shift;
  my $sa = $self->adaptor->db->get_SampleAdaptor;
  my %samples_list;
  if ($sa) {
    my $samples = $sa->fetch_all_by_dbID_list($sample_ids);
    foreach my $sample (@$samples) {
      my $sample_id = $sample->dbID;
      $samples_list{$sample_id} = $sample;
    }
  }
  return \%samples_list;
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

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

# Ensembl module for Bio::EnsEMBL::Variation::Population
#
#


=head1 NAME

Bio::EnsEMBL::Variation::Population - A population represents a phenotypic 
group, ethnic group, set of individuals used in an assay, etc.

=head1 SYNOPSIS

    # Population
    $pop = Bio::EnsEMBL::Variation::Population->new(
        -name        => 'WEST AFRICA',
        -description => 'Sub-Saharan Nations bordering Atlantic north' .
                        ' of Congo River, and Central/Southern Atlantic' .
                        ' Island Nations.');

    # print out all sub populations of a population
    # same could work for super populations

    print_sub_pops($pop);

    sub print_sub_pops {
        my $pop = shift;
        my $level = shift || 0;
        my $sub_pops = $pop->get_all_sub_Populations();

        foreach my $sp (@$sub_pops) {
            print ' ' x $level++,
                  'name: ', $sp->name(),
                  'desc: ', $sp->description(),
                  'size: ', $sp->size(),"\n";
            print_sub_pops($sp, $level);
        }
    }

=head1 DESCRIPTION

This is a class representing a population.  A population may consist of any
grouping of individuals, including phenotypic groups (e.g. people with
diabetes), ethnic groups (e.g. caucasians), individuals used in an assay
(e.g. subjects in experiment X), etc.

Populations may be arranged into an arbitrary hierarchy of sub and super
populations.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Population;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);

our @ISA = ('Bio::EnsEMBL::Storable');

=head2 new

  Arg [-dbID]           : int - unique internal identifier of the population
  Arg [-ADAPTOR]        : Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor
  Arg [-NAME]           : string - name of the population
  Arg [-DESCRIPTION]    : string - description of the population
  Arg [-SIZE]           : int - the size of the population
  Arg [-COLLECTION]     : int - flag indication if population is a collection
                                of individuals (1) or defined by geography (0)
  Arg [-SUB_POPULATIONS]: listref of Bio::EnsEMBL::Population objects 
  Example               : $pop = Bio::EnsEMBL::Variation::Population->new(
                                    -name            => 'WEST AFRICA',
                                    -description     => 'Sub-Saharan Nations bordering Atlantic north' .
                                                        ' of Congo River, and Central/Southern Atlantic' .
                                                        ' Island Nations.',
                                    -collection      => 0,
                                    -sub_populations => \@sub_pops);
  Description            : Constructor. Instantiates a new Population object
  Returntype             : Bio::EnsEMBL::Variation::Population
  Exceptions             : none
  Caller                 : general
  Status                 : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbID, $adaptor, $name, $desc, $size, $freqs, $sub_pops, $display, $collection, $display_name, $display_priority) = rearrange([
    'DBID','ADAPTOR','NAME', 'DESCRIPTION', 'SIZE', 'FREQS', 'SUB_POPULATIONS', 'DISPLAY', 'COLLECTION','DISPLAY_GROUP_NAME','DISPLAY_GROUP_PRIORITY'], @_);

  return bless {
    'dbID'                   => $dbID,
    'adaptor'                => $adaptor,
    'name'                   => $name,
    'description'            => $desc,
    'size'                   => $size,
    'freqs'                  => $freqs,
    'sub_populations'        => $sub_pops,
    'display'                => $display,
    'display_group_name'     => $display_name,
    'display_group_priority' => $display_priority,
    'collection'             => $collection,
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
  Example    : $name = $population->name()
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
  Example    : $description = $population->description()
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
               The new value to set the size attribute to
  Example    : $size = $population->size()
  Description: Getter/Setter for the size attribute
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

=head2 collection

  Arg [1]    : int $collection (optional) 
  Example    : $size = $population->size()
  Description: Getter/Setter for the collection attribute
  Returntype : Int. Returns 1 if population is a collection of individuals and 0 if it is a population defined by geography
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub collection {
  my $self = shift;
  return $self->{'collection'} = shift if (@_);
  return $self->{'collection'};
}


sub display {
  my $self = shift;
  return $self->{'display'} = shift if (@_);
  return $self->{'display'};
}

=head2 get_all_sub_Populations

  Arg [1]    : none
  Example    : foreach my $sub_pop (@{$pop->get_all_sub_Populations}) {
                 my $sub_sub_pops = $sub_pop->get_all_sub_Populations();
               }
  Description: Retrieves all populations which are conceptually a sub set
               of this population.
  Returntype : listref of Bio::EnsEMBL::Variation::Population objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_sub_Populations {
  my $self = shift;

  if (!defined($self->{'sub_populations'}) && $self->{'adaptor'}) {
    # lazy-load from database
    $self->{'sub_populations'} = $self->{'adaptor'}->fetch_all_by_super_Population($self);
  }
  return $self->{'sub_populations'} || [];
}


=head2 get_all_super_Populations

  Arg [1]    : none
  Example    : foreach my $sup_pop (@{$pop->get_all_super_Populations}) {
                 my $sup_sup_pops = $sup_pop->get_all_super_Populations();
               }
  Description: Retrieves all populations which this population is a part of
               from the database.
               Super populations may not be directly added in order to avoid
               circular references and memory leaks.  You must add
               sub_Populations instead and store this in the database.
  Returntype : listref of Bio::EnsEMBL::Variation::Population objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_super_Populations {
  my $self = shift;
  return [] if (!$self->{'adaptor'});
  # load from database - do not cache to avoid circular references (mem leak)!
  return $self->{'adaptor'}->fetch_all_by_sub_Population($self);
}


=head2 add_sub_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : $pop->add_sub_Population($sub_pop);
               $sub_pop->add_super_Population($pop);
  Description: Adds a sub population to this population.
  Returntype : none
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub add_sub_Population {
  my $self = shift;
  my $pop = shift;

  if (!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected.');
  }

  if ($pop == $self) {
    throw("Cannot add self as sub population.");
  }

  $self->{'sub_populations'} ||= [];
  push @{$self->{'sub_populations'}}, $pop;

  return $pop;
}

=head2 get_all_synonyms

  Arg [1]    : (optional) string $source - the source of the synonyms to
               return.
  Example    : @dbsnp_syns = @{$p->get_all_synonyms('dbSNP')};
               @all_syns = @{$p->get_all_synonyms()};
  Description: Retrieves synonyms for this Population. If a source argument
               is provided all synonyms from that source are returned,
               otherwise all synonyms are returned.
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_synonyms {
  my $self = shift;
  my $source = shift;

  return [] if(!$self->adaptor()); #if there is no adaptor, return empty string

  return $self->adaptor()->fetch_synonyms($self->dbID(),$source);

}

=head2 get_all_Individuals

  Arg [1]    : none
  Example    : @individuals = @{$p->get_all_individuals()};
  Description: Retrieves all Individuals belonging to this Population.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Individual objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Individuals {
  my $self = shift;

  my $ia = $self->adaptor->db->get_IndividualAdaptor;
  
  return (defined $ia ? $ia->fetch_all_by_Population($self) : []);
}

=head2 get_all_Samples

  Arg [1]    : none
  Example    : @samples = @{$p->get_all_Samples()};
  Description: Retrieves all Samples belonging to this Population.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Sample objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Samples {
  my $self = shift;
  my $sa = $self->adaptor->db->get_SampleAdaptor;
  return (defined $sa ? $sa->fetch_all_by_Population($self) : []);
}

sub _freqs_from_gts {
  my $self = shift;
  $self->{freqs} = shift @_ if @_;
  return $self->{freqs};
}

=head2 display_group_priority

  Arg [1]    : none
  Example    : $priority = $p->display_group_priority()
  Description: Retrieves the priority of the group of populations this 
               population belongs to, enabling the ordering of frequency 
               data on the Population genetics page
  Returntype : 0 or the rank of the group
  Exceptions : none
  Caller     : webcode
  Status     : experimental

=cut

sub display_group_priority{
  my $self = shift;
  $self->{display_group_priority} = shift @_ if @_;
  return $self->{display_group_priority};
}

=head2 display_group_name

  Arg [1]    : none
  Example    : $name = $p->display_group_name()
  Description: Retrieves the name of the group of populations this 
               population belongs to, to serve as the table header
               on the Population genetics page
  Returntype : null or a string representing the name of the group
  Exceptions : none
  Caller     : webcode
  Status     : experimental

=cut
sub display_group_name{
  my $self = shift;
  $self->{display_group_name} = shift @_ if @_;
  return $self->{display_group_name};
}


1;

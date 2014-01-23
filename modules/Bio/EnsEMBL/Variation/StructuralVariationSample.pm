=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

# Ensembl module for Bio::EnsEMBL::Variation::StructuralVariationSample
#
#


=head1 NAME

Bio::EnsEMBL::Variation::StructuralVariationSample - Samples for a structural variant (sample and phenotype samples).

=head1 SYNOPSIS
  
  # A study object
  $study = $study_adaptor->fetch_by_name('estd199');
  $individual = $individual_adaptor->fetch_by_name('NA18517');

  $svs = Bio::EnsEMBL::Variation::StructuralVariationSample->new
        (-individual => $individual,
         -study => $study);
  ...
  
  $svs->structural_variation->variation_name(),":", $svs->individual->name();      

=head1 DESCRIPTION

This is a class representing the sample of a structural variant
from the ensembl-variation database. The actual structural variant information is
represented by an associated Bio::EnsEMBL::Variation::StructuralVariation object. 

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::StructuralVariationSample;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Variation::BaseStructuralVariation;
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);

=head2 new

  Arg [-dbID] :
    int - unique internal identifier for variation_sample
  
  Arg [-ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::StructuralVariationSampleAdaptor
    Adaptor which provides database connectivity for this StructuralVariationSample object
  
  Arg [-INDIVIDUAL] :
    object ref - the individual object associated with the structural variant.
  
  Arg [-STRAIN] :
    object ref - the individual object (used as a strain) associated with the structural variant.
  
  Arg [-STUDY] :
    object ref - the study object describing where the annotated structural variant comes from.  
  
  Arg [_STRUCTURAL_VARIATION_ID] :
    int _ the internal id of the structural variant object associated with this
    identifier. TUsing this identifier the structural variant may be lazy-loaded from 
    the database on demand.
  
  Example    :
    $study = $study_adaptor->fetch_by_name('nstd37');

    $sva = Bio::EnsEMBL::Variation::StructuralVariationSample->new
          (-individual => $individual,
           -strain     => $strain,
           -study      => $study);

  Description: Constructor. Instantiates a new StructuralVariationSample object.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariationSample
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($dbID,$adaptor,$structural_variation_id,$individual,$strain,$study) =
    rearrange([qw(dbID ADAPTOR _STRUCTURAL_VARIATION_ID INDIVIDUAL STRAIN STUDY)],@_); 

  $self->{'dbID'}                     = $dbID;
  $self->{'adaptor'}                  = $adaptor;
  $self->{'_structural_variation_id'} = $structural_variation_id;
  $self->{'individual'}               = $individual;
  $self->{'strain'}                   = $strain;
  $self->{'study'}                    = $study;
  return $self;
}



sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 structural_variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::StructuralVariation or 
               Bio::EnsEMBL::Variation::SupportingStructuralVariation $structural_variation
  Example    : $sv = $svs->structural_variation();
  Description: Getter/Setter for the structural variant associated with this feature.
               If not set, and this StructuralVariationFeature has an associated adaptor
               an attempt will be made to lazy-load the structural variation from the
               database.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation or Bio::EnsEMBL::Variation::SupportingStructuralVariation
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub structural_variation {
  my $self = shift;

  if(@_) {
    if(!ref($_[0]) || (!$_[0]->isa('Bio::EnsEMBL::Variation::StructuralVariation') &&
                       !$_[0]->isa('Bio::EnsEMBL::Variation::SupportingStructuralVariation')
    )) {
      throw("Bio::EnsEMBL::Variation::StructuralVariation or Bio::EnsEMBL::Variation::SupportingStructuralVariation argument expected");
    }
    $self->{'_structural_variation_id'} = shift;
  }
  elsif(!defined($self->{'structural_variation'}) && $self->{'adaptor'} &&
         defined($self->{'_structural_variation_id'})) {
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


=head2 study

  Arg [1]    : Bio::EnsEMBL::Variation::Study (optional)
  Example    : $study = $svs->study()
  Description: Getter/Setter for the study object
  Returntype : Bio::EnsEMBL::Variation::Study
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub study {
  my $self = shift;
  return $self->{'study'} = shift if(@_);
  return $self->{'study'};
}


=head2 individual

  Arg [1]    : Bio::EnsEMBL::Variation::Individual (optional)
  Example    : $individual = $svs->individual()
  Description: Getter/Setter for the individual object
  Returntype : Bio::EnsEMBL::Variation::Individual
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub individual {
  my $self = shift;
  return $self->{'individual'} = shift if(@_);
  return $self->{'individual'};
}


=head2 strain

  Arg [1]    : Bio::EnsEMBL::Variation::Individual (optional)
  Example    : $strain = $svs->strain()
  Description: Getter/Setter for the individual object as a strain
  Returntype : Bio::EnsEMBL::Variation::Individual
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub strain {
  my $self = shift;
  return $self->{'strain'} = shift if(@_);
  return $self->{'strain'};
}

1;

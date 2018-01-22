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

# Ensembl module for Bio::EnsEMBL::Variation::StructuralVariation
#
#


=head1 NAME

Bio::EnsEMBL::Variation::StructuralVariation - Ensembl representation of a structural variation.

=head1 SYNOPSIS
    # A Study object
    $study = $study_adaptor->fetch_by_name('estd59');
    
    # A Source object
    $source = $source_adaptor->fetch_by_name('DGVa');
    
    # Structural variation representing a CNV
    $sv = Bio::EnsEMBL::Variation::StructuralVariation->new
       (-variation_name => 'esv234231',
        -class_so_term => 'structural_variant',
        -source => $source,
        -study => $study,
        -is_evidence => 0,
        -is_somatic => 0);

    ...

    print $sv->variation_name(), ":", $sv->var_class();

=head1 DESCRIPTION

This is a class representing a structural variation from the
ensembl-variation database. A structural variant may have a copy number variation, a tandem duplication, 
an inversion of the sequence or others structural variations. 

The position of a StructuralVariation object on the Genome is represented
by the <Bio::EnsEMBL::Variation::StructuralVariationFeature> class.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::StructuralVariation;

use Bio::EnsEMBL::Variation::BaseStructuralVariation;

our @ISA = ('Bio::EnsEMBL::Variation::BaseStructuralVariation');


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = Bio::EnsEMBL::Variation::BaseStructuralVariation->new(@_);
  return(bless($self, $class));
}

=head2 get_all_SupportingStructuralVariants

  Example     : $sv->get_all_SupportingStructuralVariants();
  Description : Retrieves all SupportingStructuralVariation associated with this structural variation.
                Return empty list if there are none.
  Returntype  : reference to list of Bio::EnsEMBL::Variation::SupportingStructuralVariation objects
  Exceptions  : None
  Caller      : general
  Status      : Stable

=cut

sub get_all_SupportingStructuralVariants {
  my $self = shift;
  
  if (defined ($self->{'adaptor'})){
    my $ssv_adaptor = $self->{'adaptor'}->db()->get_SupportingStructuralVariationAdaptor();
    return $ssv_adaptor->fetch_all_by_StructuralVariation($self);
  }
  warn("No variation database attached");
  return [];
}


=head2 summary_as_hash

  Example       : $sv_summary = $sv->summary_as_hash();
  Description   : Retrieves a textual summary of this StructuralVariation object.
  Returns       : hashref of descriptive strings

=cut

sub summary_as_hash {
  my $self = shift;
  my %summary;
  $summary{'display_id'} = $self->display_id;
  $summary{'study_name'} = $self->study->name;
  $summary{'study_description'} = $self->study->description;
  $summary{'class'} = $self->var_class;
  return \%summary;

}

1;

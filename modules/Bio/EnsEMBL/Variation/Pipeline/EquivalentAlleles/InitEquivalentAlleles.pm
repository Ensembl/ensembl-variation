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




=head1 NAME 

Bio::EnsEMBL::Variation::Pipeline::EquivalentAlleles::InitEquivalentAlleles

=head1 DESCRIPTION

Initiation module for Equivalent Allele pipeline

=cut

package Bio::EnsEMBL::Variation::Pipeline::EquivalentAlleles::InitEquivalentAlleles;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);



sub fetch_input {
   
  my $self = shift;

  $self->delete_previous();
  
  my $core_dba = $self->get_species_adaptor('core');

  ## Fetch all top level slices and allocate into overlapping regions for analysis

  my $sa  = $core_dba->get_SliceAdaptor();
  my $slices = $sa->fetch_all('toplevel');

  my @locations;

  my $offset       = $self->required_param('region_size') - $self->required_param('overlap');

  $self->warning(scalar @{$slices} .'  slices found');

  foreach my $slice (@{$slices}){

    my $start  = 1;
    my $length = $slice->length();

    while($start < $length){

      my $end = $start + $self->required_param('region_size');

      ## store slice required
      my $location = $slice->seq_region_name() . ":" .$start. "-" . $end ;

      push @locations, { location => $location }; 
      $start  = $start + $offset; 
    }
  }
  $self->param('locations', \@locations);
}

## delete results of previous analysis
sub delete_previous{

  my $self = shift;

  my $var_dba = $self->get_species_adaptor('variation');
  $var_dba->dbc->do(qq[ delete from variation_attrib 
                        where attrib_id in(
                          select attrib_id from attrib where value ='co-located allele'
                          )
                       ]);

}

sub write_output {
    
  my $self = shift;

  ## Check each region 
  my $locations =  $self->param('locations');
  $self->warning(scalar @{$locations} .'  jobs to do');
  $self->dataflow_output_id($locations, 2);  


  ## run basic checks when everything is updated
  $self->dataflow_output_id($self->param('finish_equivalent_alleles'), 3);

  return;
}

1;



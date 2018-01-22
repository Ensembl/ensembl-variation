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

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::UnmappedVariant;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);


=head1 NAME

  Bio::EnsEMBL::Variation::Pipeline::VariantQC::UnmappedVariant

=head1 DESCRIPTION

This module extracts a list of ids for variations without a mapped location on the currents genomic sequence and enters them in the failed variation table.
It is run as a seperate independant process alongside the main variant QC 

=cut

sub run {
    
    my $self = shift;
    
    my $first = $self->required_param('start_id'); 
    my $last  = $first + $self->required_param('batch_size') -1; 
    if($first ==1){$last--;} 
    
    my $var_dba      = $self->get_species_adaptor('variation');

    my $fail_ins_sth = $var_dba->dbc->prepare(qq[insert into failed_variation_working
                                                 (variation_id, failed_description_id)
                                                 values (?,?)
                                                 ]);       
            
    
  
    my $mapfail_extr_sth = $var_dba->dbc->prepare(qq[select variation.variation_id 
                                                from variation 
                                                where variation.variation_id between $first and $last
                                                and  variation.variation_id not in( select variation_feature.variation_id from variation_feature) 
                                                ]);
    

    ## export current data 
    $mapfail_extr_sth->execute()||die;

    my $fails = $mapfail_extr_sth->fetchall_arrayref();

    return unless(defined $fails >[0]->[0]);

	
    ## write to fail table
    foreach my $row( @{$fails}){
	
	$fail_ins_sth->execute($row->[0], 5)|| die "ERROR inserting variation fails info\n";
    }
}



1;

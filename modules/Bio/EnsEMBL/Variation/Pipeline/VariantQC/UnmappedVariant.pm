=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

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

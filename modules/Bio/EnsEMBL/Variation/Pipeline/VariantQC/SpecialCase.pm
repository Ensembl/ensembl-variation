
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



=head1 NAME 

Bio::EnsEMBL::Variation::Pipeline::VariantQC::SpecialCase

=head1 DESCRIPTION

ammends failure status in special cases:
    - variants with pubmed ids should not be failed
    - variants in PAR regions should not be failed on multiple map locations

It is quicker to handle these in bulk at the end

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::SpecialCase;


use strict;
use warnings;


use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);


our $DEBUG   = 1;
 
=head2 run

 Check what needs doing & do it

=cut
sub run {

  my $self = shift;
   
  
  if ( $self->param('run_PAR_check') ==1 ){ $self->check_PAR_variants();}

  if ( $self->param('run_Pubmed_check') ==1 ){ $self->check_Pubmed_variants();}

}

=head2 check_PAR_variants

  find variants in PAR regions and remove 'failed on multiple map locations' status

=cut
sub check_PAR_variants{

    my $self = shift;

    my $var_dba = $self->get_species_adaptor('variation');

    ## check if internal production db is available
    my $int_dba ;
    eval{ $int_dba = $self->get_species_adaptor('intvar');};

    unless (defined $int_dba){
	$self->warning('No internal database connection found to look up PAR variants '); 
	return;
    }

    my $fail_delete_sth = $var_dba->dbc->prepare(qq[delete from failed_variation where failed_description_id = 19 and variation_id in
                                                   select vf_x.variation_id
                                                   from  variation_feature vf_x, variation_feature vf_y, seq_region sr_x, seq_region sr_y
                                                   where vf_x.map_weight = 2
                                                   and vf_x.seq_region_id = sr_x.seq_region_id  
                                                   and sr_x.name ='X' and vf_x.seq_region_start between ? AND ?
                                                   and vf_y.variation_id = vf_x.variation_id
                                                   and vf_y.seq_region_id = sr_y.seq_region_id  
                                                   and sr_y.name ='Y' and vf_y.seq_region_start ? and ? )]);


   my $par_ext_sth = $int_dba->dbc->prepare(qq[ select X_region_start, X_region_end, Y_region_start, Y_region_end 
                                                from PAR_location
                                                where species = ? and is_current = ?
                                              ]);

    $par_ext_sth->execute($self->required_param('species'), 1);
    my $regions = $par_ext_sth->fetchall_arayref();

    foreach my $l (@{$self->param('PAR')}){

	$fail_delete_sth->execute($l->[0], $l->[1], $l->[2], $l->[3] )||die "Failed to update PAR variants\n"; ;	
    }
}
	
=head2  check_Pubmed_variants

 Remove all fail statuses from cited variants for which dbSNP holds Pubmed ids

=cut
sub check_Pubmed_variants{

    my $self = shift;

    my $var_dba = $self->get_species_adaptor('variation');

    ## check if there are any cited variants
    my $table_check_sth = $var_dba->dbc->prepare(qq[ show tables like '%pubmed_variation%']);
    $table_check_sth->execute()||die;
    my $present =  $table_check_sth->fetchall_arrayref();

    if( defined $present->[0]->[0]){
        $self->warning('Setting cited variants to non-failed '); 
 
	$var_dba->dbc->do(qq[ delete from failed_variation_working 
                              where variation_id in (select variation_id from pubmed_variation)
                              ]);
    }
    else{
      $self->warning('No cited variants found to unfailed '); 
    }
}


1;

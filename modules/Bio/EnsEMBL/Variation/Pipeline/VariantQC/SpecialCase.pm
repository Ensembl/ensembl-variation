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


 
=head2 run

 Check what needs doing & do it

=cut
sub run {

  my $self = shift;

  if ( $self->param('run_PAR_check') ==1 ){ $self->check_PAR_variants();}

  if ( $self->param('run_Pubmed_check') ==1 ){ $self->check_Pubmed_variants();}

  if ( $self->required_param('species') =~/sus_scrofa/ ){ $self->add_synonym();}

  if ( $self->required_param('species') =~/mus_musculus|gallus_gallus|rattus_norvegicus|canis_familiaris|homo_sapiens/ ){ $self->set_display();}


}

=head2 check_PAR_variants

  find variants in PAR regions and remove 'failed on multiple map locations' status

=cut
sub check_PAR_variants{

    my $self = shift;

    my $var_dba = $self->get_species_adaptor('variation');

    ## check if internal production db is available
    my $int_dba ;
    eval{ $int_dba = $self->get_adaptor('multi', 'intvar');};

    unless (defined $int_dba){
        $self->warning('No internal database connection found to look up PAR variants '); 
        return;
    }

    my $fail_delete_sth = $var_dba->dbc->prepare(qq[delete from failed_variation where failed_description_id = 19 and variation_id in
                                                   (select vf_x.variation_id from  variation_feature vf_x, variation_feature vf_y, seq_region sr_x, seq_region sr_y
                                                   where vf_x.map_weight = 2
                                                   and vf_x.seq_region_id = sr_x.seq_region_id  
                                                   and sr_x.name ='X' and vf_x.seq_region_start between ? and ?
                                                   and vf_y.variation_id = vf_x.variation_id
                                                   and vf_y.seq_region_id = sr_y.seq_region_id  
                                                   and sr_y.name ='Y' and vf_y.seq_region_start between ? and ? )]);


   my $par_ext_sth = $int_dba->dbc->prepare(qq[ select X_region_start, X_region_end, Y_region_start, Y_region_end 
                                                from PAR_location
                                                where species = ? and is_current = ?
                                              ]);

    $par_ext_sth->execute($self->required_param('species'), 1);
    my $regions = $par_ext_sth->fetchall_arrayref();

    foreach my $l (@{$regions}){

	$fail_delete_sth->execute($l->[0], $l->[1], $l->[2], $l->[3] )||die "Failed to update PAR variants\n"; ;	
    }
}
	
=head2  check_Pubmed_variants

 Remove all fail statuses from cited variants for which dbSNP holds Pubmed ids
 Citations ar held on studies 

=cut
sub check_Pubmed_variants{

    my $self = shift;

    my $var_dba = $self->get_species_adaptor('variation');

    ## check if there are any cited variants
    $self->warning('Setting cited variants to non-failed ');
 
    $var_dba->dbc->do(qq[ delete from failed_variation_working 
                          where variation_id in (select variation_id from variation_citation)
                        ]);

}

=head2 set_display

  set display status on expected strains/individuals to enable mart filtering & read coverage viewing
  report missing or duplicated samples

=cut
sub set_display{

    my $self = shift;

    my $var_dba = $self->get_species_adaptor('variation');
    my $dir = $self->required_param('pipeline_dir');
    open my $report, ">", "$dir/Display_report.txt" || die "Failed to open Display_report.txt : $!\n";


    ## check if internal production db is available
    my $int_dba ;
    eval{ $int_dba = $self->get_adaptor('multi', 'intvar');};

    unless (defined $int_dba){
	$self->warning('No internal database connection found to look up expected strains '); 
	return;
    }


    my $display_update_sth = $var_dba->dbc->prepare(qq[update individual set display = ? where name = ? ]);

    #### ADAPT TO NEW SCHEMA
    ## check individuals are neither missing or duplicated
    my $individual_check_sth =  $var_dba->dbc->prepare(qq[ select count(*) from individual
                                                           where individual.name = ?
                                                          ]);


    my $individual_ext_sth = $int_dba->dbc->prepare(qq[ select name, display
                                                        from individual_display_info
                                                        where species = ? 
                                                      ]);

   


    $individual_ext_sth->execute($self->required_param('species'));
    my $individual = $individual_ext_sth->fetchall_arrayref();

    foreach my $l (@{$individual}){

	$individual_check_sth->execute( $l->[0] )||die "Failed to check displayable individuals \n";
	my $count = $individual_check_sth->fetchall_arrayref();
	if($count->[0]->[0] == 1){
            ## set display status on individual
	    print $report "Single individual seen - setting display for: $l->[0]\n";
        }
        elsif ($count->[0]->[0] ==0){
	    print $report "Error : individual $l->[0] missing from new import\n";
	    next;
	}
	else{
	    print $report "Error : individual $l->[0] duplicated (x $count->[0]->[0]) in new import - setting display for all entries\n";
	}
	$display_update_sth->execute($l->[1], $l->[0] )||die "Failed to update display status \n";
	
    }
    close $report;
}

## Pig consortium variation names are supported - copy to new releases
## mapping on ss ids incase of rs-demerging
 
sub add_synonym{


    my $self = shift;

    my $var_dba     = $self->get_species_adaptor('variation');
    my $var_adaptor = $var_dba->get_VariationAdaptor; 
    
    ## using ignore here to allow for re-runs of pipeline; could create new working table
    my $synonym_ins_sth   = $var_dba->prepare(qq[ insert ignore into variation_synonym (variation_id, name, source_id) values (?,?,?) ]);  

    my %source_data = ( "name"           => "Pig SNP Consortium",
			"version"        =>  '\N',
			"description"    => 'PorcineSNP60 BeadChip', 
			"somatic_status" => "germline"
	);

    my $source_id = get_source($var_dba ,\%source_data);

    my $all_synonym = $self->get_synonym();
    
     
    foreach my $l(@{$all_synonym}){

	my $var = $var_adaptor->fetch_by_name($l->[0]);
	if(defined $var){

	    $synonym_ins_sth->execute( $var->dbID, $l->[1], $source_id );
	}
	else{
	    $self->warning( 'No database id for expected pig variant $l->[0] ');
	}
    }
}


sub get_source{

    my $var_dba     = shift;
    my $source_data = shift;

    my $source_ins_sth    = $var_dba->prepare(qq[ insert into source (name, version, description) values (?,?,?) ]);
    my $source_ext_sth    = $var_dba->prepare(qq[ select source_id from source where name = ? ]);

    $source_ext_sth->execute($source_data->{name})||die;
    my $id = $source_ext_sth->fetchall_arrayref();

    unless(defined $id ->[0]->[0]){

	$source_ins_sth->execute($source_data->{name}, $source_data->{version}, $source_data->{description} )||die;
	$id = $source_ext_sth->fetchall_arrayref();
    }
    return  $id ->[0]->[0];
}




sub get_synonym{
    
    my $self = shift;

    my $int_dba  = $self->get_adaptor('multi', 'intvar');

    my $synonym_ext_sth = $int_dba->prepare(qq[ select ss_name, synonym_name from pig_synonym ]);

    $synonym_ext_sth->execute();
    my $all_syn = $synonym_ext_sth->fetchall_arrayref();

    return $all_syn;
}

1;

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

Bio::EnsEMBL::Variation::Pipeline::VariantQC::SpecialCase

=head1 DESCRIPTION

ammends failure status in special cases and imports static data from production db
    - variants in PAR regions should not be failed on multiple map locations
    - some samples need their display status setting
    - chip set information is imported


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

  if ( $self->required_param('species') =~/homo_sapiens/ ){ $self->check_PAR_variants();}


  if ( $self->param('run_Pubmed_check') ==1 ){ $self->check_Pubmed_variants();}

  if ( $self->required_param('species') =~/sus_scrofa/ ){ $self->add_synonym();}

  if ( $self->required_param('species') =~/mus_musculus|gallus_gallus|rattus_norvegicus|canis_familiaris/ ){ $self->set_display();}

  $self->add_chip_info();

  $self->check_individual_names()
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

 Set display statuses on cited variants for which dbSNP holds Pubmed ids

=cut
sub check_Pubmed_variants{

    my $self = shift;

    my $var_dba = $self->get_species_adaptor('variation');

    ## check if there are any cited variants
    $self->warning('Setting cited variants to displayable ');
 
    $var_dba->dbc->do(qq[ update variation_working 
                          set display = 1
                          where variation_id in (select variation_id from variation_citation)
                        ]);

}

=head2 set_display

  set display status on expected strains/samples to enable mart filtering & read coverage viewing
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


    my $display_update_sth = $var_dba->dbc->prepare(qq[update sample set display = ? where name = ? ]);

    ## check samples are neither missing or duplicated
    my $sample_check_sth =  $var_dba->dbc->prepare(qq[ select count(*) from sample
                                                           where sample.name = ?
                                                          ]);


    my $individual_ext_sth = $int_dba->dbc->prepare(qq[ select name, display
                                                        from individual_display_info
                                                        where species = ? 
                                                      ]);

   


    $individual_ext_sth->execute($self->required_param('species'));
    my $samples = $individual_ext_sth->fetchall_arrayref();

    foreach my $l (@{$samples}){

	$sample_check_sth->execute( $l->[0] )||die "Failed to check displayable samples \n";
	my $count = $sample_check_sth->fetchall_arrayref();
	if($count->[0]->[0] == 1){
            ## set display status on sample
	    print $report "Single sample seen - setting display for: $l->[0]\n";
        }
        elsif ($count->[0]->[0] ==0){
	    print $report "Error : sample $l->[0] missing from new import\n";
	    next;
	}
	else{
	    print $report "Error : sample $l->[0] duplicated (x $count->[0]->[0]) in new import - setting display for all entries\n";
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
    my $synonym_ins_sth   = $var_dba->dbc->prepare(qq[ insert ignore into variation_synonym (variation_id, name, source_id) values (?,?,?) ]);

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

    my $source_ins_sth    = $var_dba->dbc->prepare(qq[ insert into source (name, version, description) values (?,?,?) ]);
    my $source_ext_sth    = $var_dba->dbc->prepare(qq[ select source_id from source where name = ? ]);

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

    my $synonym_ext_sth = $int_dba->dbc->prepare(qq[ select ss_name, synonym_name from pig_synonym ]);

    $synonym_ext_sth->execute();
    my $all_syn = $synonym_ext_sth->fetchall_arrayref();

    return $all_syn;
}


sub add_chip_info{

    my $self = shift;

    my $int_dba  = $self->get_adaptor('multi', 'intvar');

    my $chip_ext_sth = $int_dba->dbc->prepare(qq[ select name, shortname, description
                                                  from genotyping_chip
                                                  where species = ?
                                                 ]);

    $chip_ext_sth->execute( $self->required_param('species') );
    my $dat = $chip_ext_sth->fetchall_arrayref();

    return unless defined $dat->[0]->[0];


    ## enter chip info if any is available

    my $var_dba = $self->get_species_adaptor('variation');

    my $attrib_ext_sth = $var_dba->dbc->prepare(qq[select attrib_id from attrib, attrib_type
                                              where attrib.value = ? 
                                              and attrib.attrib_type_id = attrib_type.attrib_type_id  
                                              and attrib_type.code = ?]);

    my $meta_ins_sth   = $var_dba->dbc->prepare(qq[insert into meta (species_id, meta_key, meta_value ) values (?,?,?) ]);

    my $set_ins_sth = $var_dba->dbc->prepare(qq[ insert into variation_set
                                            ( name, description, short_name_attrib_id ) 
                                             values (?,?,?)]);

    my $setstr_ins_sth = $var_dba->dbc->prepare(qq[ insert into variation_set_structure
                                               (variation_set_super,variation_set_sub) values (?,?)]);


    ## set up default menus first
    $meta_ins_sth->execute( 1, "web_config", "source#Sequence variants (dbSNP and all other sources)#All sequence variants#variation_feature_variation#var");
    $meta_ins_sth->execute( 1, "web_config", "source#dbSNP variants#variation_feature_variation_dbSNP#var");
    $meta_ins_sth->execute( 1, "web_config", "menu_sub#Sequence variants##var#variants");
    $meta_ins_sth->execute( 1, "web_config", "set#All failed variants#All failed variants#variation_set_fail_all#failed");
    $meta_ins_sth->execute( 1, "web_config", "menu#Failed variants##failed#");


    ## enter parent set 
    $attrib_ext_sth->execute("all_chips", "short_name");
    my $super_att = $attrib_ext_sth->fetchall_arrayref();
   
    return unless defined $super_att->[0]->[0];

    $set_ins_sth->execute("Genotyping chip variants", 
			  "Variants which have assays on commercial chips held in ensembl",
			  $super_att->[0]->[0] );
    
    my $super_set =$var_dba->dbc->db_handle->last_insert_id(undef, undef, qw(variation_set variation_set_id))
	|| die "no insert id for chip super set\n";

    $meta_ins_sth->execute( 1, "web_config",
			   "set#All variants on genotyping chips#All genotyping chips#variation_set_all_chips#var_other");

    $meta_ins_sth->execute( 1, "web_config",
			   "menu#Arrays and other##var_other#");
    




    ### enter individual chips as sets

    foreach my $chip (@{$dat}){

	$attrib_ext_sth->execute( $chip->[1], 9);
	my $att = $attrib_ext_sth->fetchall_arrayref();

	$set_ins_sth->execute( $chip->[0], $chip->[2], $att->[0]->[0] );
    
	my $sub_set =$var_dba->dbc->db_handle->last_insert_id(undef, undef, qw(variation_set variation_set_id))
	    || die "no insert id for chip sub set\n";

	$setstr_ins_sth->execute( $super_set , $sub_set);

	my $meta_string = "set#" . $chip->[0]. "#" . $chip->[0]. "#variation_set_" .  $chip->[1] . "#var_other";
	$meta_ins_sth->execute( 1, "web_config", $meta_string );

	populate_set($var_dba, $int_dba, $chip->[1], $sub_set);
    }


}


sub populate_set{

    my $var_dba   = shift; 
    my $int_dba   = shift; 
    my $chip_name = shift;
    my $set_id    = shift;


    ### check the chip content is loaded 
    my $check_avail = $int_dba->dbc->prepare(qq[show tables like '$chip_name']);
    $check_avail->execute();
    my $avail =  $check_avail->fetchall_arrayref();

    return unless defined $avail->[0]->[0];


    open my $set_log, ">>variation_set_report.txt"|| die "Failed to open variation_set_report.txt: $!";

    my $v_adaptor = $var_dba->get_VariationAdaptor;

    my $vset_ins_sth = $var_dba->dbc->prepare(qq[ insert into variation_set_variation
                                             (variation_id, variation_set_id) values( ?,?)]); 


    ## extract info from production db
    my $list_ext_sth =  $int_dba->dbc->prepare(qq[select name from  $chip_name]);
    $list_ext_sth->execute();
    while(my $var = $list_ext_sth->fetchrow_arrayref()){

	my $var_ob = $v_adaptor->fetch_by_name($var->[0]);

	unless (defined $var_ob){
	    ### report and skip if variant not found
	    print $set_log "$var->[0] not found for chip $chip_name\n";
	    next;
	}
	## add variant to chip set
	$vset_ins_sth->execute( $var_ob->dbID(), $set_id );
    }

}
## individuals do not have names in dbSNP, but samples do
## individual name is populated by the first sample name seen
## report and check any individuals with multiple names
sub check_individual_names{

    my $self = shift;

    my $var_dba = $self->get_species_adaptor('variation');
    my $name_ext_sth = $var_dba->dbc->prepare(qq[ select individual_id, name from sample]);
   
    $name_ext_sth->execute()||die;
    my $dat = $name_ext_sth->fetchall_arrayref();

    my %names;
    foreach my $l (@{$dat}){
        push @{$names{$l->[0]}}, $l->[1];
    }

    open my $ind_log, ">>individual_name_report.txt"|| die "Failed to open individual_name_report.txt: $!";
    print $ind_log "Individuals with more than one sample name to choose from: \n\n";
    foreach my $ind( keys %names){
       next if @{$names{$ind}} ==1;
       print $ind_log "$ind\t". join(",", @{$names{$ind}}) . "\n";
    }    
}

1;

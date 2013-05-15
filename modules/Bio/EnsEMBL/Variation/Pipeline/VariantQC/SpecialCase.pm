
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

  if ( $self->required_param('species') =~/sus_scrofa/ ){ $self->add_synonym();}

  if ( $self->required_param('species') =~/mus_musculus|gallus_gallus|rattus_norvegicus|canis_familiaris/ ){ $self->set_strain_display();}

  if ( $self->required_param('species') =~/gallus_gallus|canis_familiaris/ ){ $self->fake_read_coverage();}

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

=head2 set_strain_display

  set display status on expected strains to enable mart filtering & read coverage viewing
  report missing or duplicated samples

=cut
sub set_strain_display{

    my $self = shift;

    my $var_dba = $self->get_species_adaptor('variation');
    my $dir = $self->required_param('pipeline_dir');
    open my $report, ">", "$dir/Strain_report.txt" || die "Failed to open Strain_report.txt : $!\n";


    ## check if internal production db is available
    my $int_dba ;
    eval{ $int_dba = $self->get_adaptor('multi', 'intvar');};

    unless (defined $int_dba){
	$self->warning('No internal database connection found to look up expected strains '); 
	return;
    }

#    my $display_update_sth = $var_dba->dbc->prepare(qq[update sample set display = ? where name = ? ]);
    my $display_update_sth = $var_dba->dbc->prepare(qq[update individual set display = ? where name = ? ]);

    #### ADAPT TO NEW SCHEMA
    ## check strains are neither missing or duplicated
#    my $strain_check_sth =  $var_dba->dbc->prepare(qq[ select count(*) from sample, individual
#                                                       where sample.name =?
#                                                       and sample.sample_id = individual.sample_id
#                                                      ]);

    my $strain_check_sth =  $var_dba->dbc->prepare(qq[ select count(*) from individual
                                                       where individual.name = ?
                                                      ]);


    my $strain_ext_sth = $int_dba->dbc->prepare(qq[ select name, display
                                                    from strain_info
                                                    where species = ? 
                                                  ]);

   


    $strain_ext_sth->execute($self->required_param('species'));
    my $strain = $strain_ext_sth->fetchall_arrayref();

    foreach my $l (@{$strain}){
	print $report "Checking expected strain $l->[0]\n";
	$strain_check_sth->execute( $l->[0] )||die "Failed to check displayable individuals \n";
	my $count = $strain_check_sth->fetchall_arrayref();
	if($count->[0]->[0] == 1){
            ## set display status on individual
	    print $report "Single individual seen - setting display for strain $l->[0]\n";
        }
        elsif ($count->[0]->[0] ==0){
	    print $report "Error : strain $l->[0] missing from new import\n";
	    next;
	}
	else{
	    print $report "Error : strain $l->[0] duplicated (x $count->[0]->[0]) in new import - setting display for all entries\n";
	}
	$display_update_sth->execute($l->[1], $l->[0] )||die "Failed to update display status \n";
	
    }
    close $report;
}

### Three species have near genome-coverage data, but precise coverage not known
### Fake to cover entire length of each seq_region, pending view change

sub fake_read_coverage{


    my $self = shift;

    my $var_dba  = $self->get_species_adaptor('variation');
    my $core_dba = $self->get_species_adaptor('core');

    my $len_extr_sth = $core_dba->dbc->prepare(qq[ select seq_region_id,length from seq_region ]);


#    my $sam_extr_sth = $var_dba->dbc->prepare(qq[ select individual.sample_id from individual,sample
#                                                  where individual.sample_id = sample.sample_id
#                                                  and sample.display !='UNDISPLAYABLE']);
     my $sam_extr_sth = $var_dba->dbc->prepare(qq[ select individual_id from individual
                                                  where display !='UNDISPLAYABLE']);
   
    
    my $cov_ins_sth = $var_dba->dbc->prepare(qq[ insert into  read_coverage
                                                 (seq_region_id,seq_region_start,seq_region_end, level, individual_id )
                                                 values(?,1,?,1,?)]);
    
    my %seq;
    ## extract sequence lengths from core database
    $len_extr_sth->execute()||die "Failed to extrace seq ids and lengths\n";
    my $seq = $len_extr_sth->fetchall_arrayref();
    foreach my $l (@{$seq}){
        $seq{$l->[0]} = $l->[1];
    }
    
    ## extract samples with display options from variation database
    $sam_extr_sth->execute()||die "Failed to extrace seq ids and lengths\n";
    my $sam = $sam_extr_sth->fetchall_arrayref();

    ## add fake read coverage info
    foreach my $l (@{$sam}){     
        foreach my $seq_id(keys %seq){            
            
            $cov_ins_sth->execute($seq_id ,$seq{$seq_id}, $l->[0])||die "Failed to load read_coverage \n";
        }
    }
}

## Pig consortium variation names are supported - copy to new releases
 
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

    my $synonym_ext_sth = $int_dba->prepare(qq[ select rs_name, synonym_name from pig_synonym ]);

    $synonym_ext_sth->execute();
    my $all_syn = $synonym_ext_sth->fetchall_arrayref();

    return $all_syn;
}

1;

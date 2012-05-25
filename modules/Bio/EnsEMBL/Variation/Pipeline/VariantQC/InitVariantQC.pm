=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.




=head1 NAME 

Bio::EnsEMBL::Variation::Pipeline::VariantQC::InitVariantQC

=head1 DESCRIPTION

Initiation module for variant QC eHive process

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::InitVariantQC;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub fetch_input {
   
    my $self = shift;

 
    my $core_dba = $self->get_species_adaptor('core');
    my $var_dba  = $self->get_species_adaptor('variation');
    
    my $start_at_variation_id = $self->required_param('start_at_variation_id');

    ### new variation feature tables to write to
    unless($self->required_param('create_working_tables') == 0){
	create_working_tables( $var_dba);
    }
    if($self->required_param('create_map_table') ==1){
	### create temp table to hold number of mapping to reference genome
	create_map_weight_table($core_dba,$var_dba);
    }

    
    ### get max variation_id to set up batches for processing
    my $data_ext_sth = $var_dba->dbc->prepare(qq[ SELECT MAX(variation_id)
                                                  FROM variation    ]);       
    
    $data_ext_sth->execute() || die "Failed to extract variation_id\n";
    my $max_id = $data_ext_sth->fetchall_arrayref();


    ### bin detailed qc 
    
    my @qc_start_id;

    my $start_from  = int( $start_at_variation_id / $self->param('qc_batch_size'));
    my $qc_var_jobs = int( $max_id->[0]->[0] / $self->param('qc_batch_size') ); 

    $self->warning( 'Defining jobs, $start_from - $qc_var_jobs batch size: ' . $self->param('qc_batch_size') );

    for my $n ( $start_from .. $qc_var_jobs){
	my $start = $n *  $self->param('qc_batch_size');
	push @qc_start_id, {start_id => $start};
    }
    $self->param('qc_start_ids', \@qc_start_id);



    ## bin unmapped var check in larger chunks
    
    $self->warning( 'Running $jobs jobs, batch size: ' . $self->param('unmapped_batch_size') );
    my   @unmapped_start_id;
    
    my $start_unmapped_from = int( $start_at_variation_id / $self->param('unmapped_batch_size') );
    my $unmapped_var_jobs   = int( $max_id->[0]->[0]      / $self->param('unmapped_batch_size') ); 

    for my $n ( $start_unmapped_from .. $unmapped_var_jobs){
	
	my $start = $n *  $self->param('unmapped_batch_size');
	push @unmapped_start_id, {start_id => $start};
    }
    $self->param('unmapped_start_ids', \@unmapped_start_id);



}

sub create_working_tables{

    my $var_dba = shift;

    ## table to hold variation feature info after fliping & ref allele aissignment
    $var_dba->dbc->do(qq{CREATE TABLE variation_feature_working like variation_feature});

    ## table to hold allele info after fliping 
    $var_dba->dbc->do(qq{CREATE TABLE allele_working like allele});

    ## table to hold failed variations     ## TEMP FOR DEBUG
    $var_dba->dbc->do(qq{CREATE TABLE failed_variation_working like failed_variation});

    ## table to hold failed alleles        ## TEMP FOR DEBUG
    $var_dba->dbc->do(qq{CREATE TABLE failed_allele_working like failed_allele});

    ## table to hold list of flipped variation_ids     ## TEMP FOR DEBUG
    $var_dba->dbc->do(qq{CREATE TABLE variation_to_reverse_working (variation_id int(10) unsigned not null) });



}

## create a look-up table for the number of times a variant maps to the reference
sub create_map_weight_table{
    
    my $core_dba  = shift;
    my $var_dba   = shift;

    #### is it better to use not in () or temp columns??
    my $ref_ext_sth = $core_dba->dbc->prepare(qq [select sra.seq_region_id 
                                            from seq_region_attrib sra, attrib_type at 
                                            where sra.attrib_type_id=at.attrib_type_id 
                                            and at.name="Non Reference"]);
    
     $ref_ext_sth->execute()|| die "Failed to extract ref/non ref status for seq_regions\n";
     my $non_ref = $ref_ext_sth->fetchall_arrayref();
    
    $var_dba->dbc->do(qq{alter table seq_region add column is_reference Tinyint(1) default 1});
    my $sr_status_sth = $var_dba->dbc->prepare(qq[update seq_region set is_reference = 0 where seq_region_id  =?]);

    foreach my $srid (@{$non_ref}){
	$sr_status_sth->execute($srid->[0]) || die "ERROR updating seq regions\n";
    }


    ### Is this slow?  - 13mins for human; could chunk it..
    #create a temporary table to store the map_weight, that will be deleted by the last process
    $var_dba->dbc->do(qq[ CREATE TABLE tmp_map_weight_working
                SELECT variation_id, count(*) as count
                FROM   variation_feature,seq_region
                WHERE  variation_feature.seq_region_id = seq_region.seq_region_id
                AND    seq_region.is_reference =1
                GROUP BY variation_id]
	       );

    $var_dba->dbc->do(qq{ALTER TABLE tmp_map_weight_working 
		      ADD UNIQUE INDEX variation_idx(variation_id)});

    #add additional variation_ids only appear in haplotype chromosomes #removed IGNORE - test!!
    $var_dba->dbc->do(qq{INSERT INTO tmp_map_weight_working
		  SELECT variation_id, count(*) as count
		  FROM   variation_feature,seq_region
		  WHERE  variation_feature.seq_region_id = seq_region.seq_region_id
                  AND    seq_region.is_reference =0
		  GROUP BY variation_id});
    
    ## clean up seq_region table
    $var_dba->dbc->do(qq{alter table seq_region drop column is_reference});

    ### test above & not in (tmp_map_weight) approach

}

sub write_output {
    
    my $self = shift;
   
    ## No map fails - larger bins used as check is very quick

    unless ($self->param('run_unmapped_var') == 0){
	my $unmapped_start_ids =  $self->param('unmapped_start_ids');
	
	$self->warning(scalar @{$unmapped_start_ids} .' unmapped_variant_qc jobs to do');  
  
	$self->dataflow_output_id($unmapped_start_ids, 3);     
    }

    ## Variant QC - bin start positions supplied

    unless ($self->param('run_variant_qc') == 0){
	my $qc_start_ids =  $self->param('qc_start_ids');
	
	$self->warning(scalar @{$qc_start_ids} .' variant_qc jobs to do');    	

	$self->dataflow_output_id($qc_start_ids, 2);  
    }
    

    ## run basic checks when everything is updated
	
    $self->dataflow_output_id($self->param('finish_variation_qc'), 4);
   
    
    return;
}

1;

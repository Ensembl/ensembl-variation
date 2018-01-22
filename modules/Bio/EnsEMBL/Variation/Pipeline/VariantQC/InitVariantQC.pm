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

Bio::EnsEMBL::Variation::Pipeline::VariantQC::InitVariantQC

=head1 DESCRIPTION

Initiation module for variant QC eHive process

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::InitVariantQC;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use ImportUtils qw(dumpSQL create_and_load );

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


    ## look up variation_set_id for failed variants once
    my $failed_set_id = add_failed_variation_set( $var_dba );

    $self->required_param( 'failed_set_id', $failed_set_id );
    $self->warning( 'Got failed set id ',  $failed_set_id);

    ### get max variation_id to set up batches for processing
    my $data_ext_sth = $var_dba->dbc->prepare(qq[ SELECT MAX(variation_id)
                                                  FROM variation    ]);       
    
    $data_ext_sth->execute() || die "Failed to extract variation_id\n";
    my $max_id = $data_ext_sth->fetchall_arrayref();


    ### bin detailed qc 
    
    my @qc_start_id;

    my $start_from  = int( $start_at_variation_id / $self->param('qc_batch_size'));
    my $qc_var_jobs = int( $max_id->[0]->[0] / $self->param('qc_batch_size') ); 

    $self->warning( 'Defining jobs, '. $start_from .'-' . $qc_var_jobs .' batch size: ' . $self->param('qc_batch_size') );

    for my $n ( $start_from .. $qc_var_jobs){
        my $start = $n *  $self->param('qc_batch_size');
        push @qc_start_id, {start_id => $start};
    }
    $self->param('qc_start_ids', \@qc_start_id);



    ## bin unmapped var check in larger chunks
    
    my @unmapped_start_id;
    
    my $start_unmapped_from = int( $start_at_variation_id / $self->param('unmapped_batch_size') );
    my $unmapped_var_jobs   = int( $max_id->[0]->[0]      / $self->param('unmapped_batch_size') ); 

    for my $n ( $start_unmapped_from .. $unmapped_var_jobs){
        
        my $start = $n *  $self->param('unmapped_batch_size');
        push @unmapped_start_id, {start_id => $start};
    }
    $self->param('unmapped_start_ids', \@unmapped_start_id);



}
=head2 create_working_tables

copies of key tables created to hold post processed data 

=cut

sub create_working_tables{

  my $var_dba = shift;

  ## table to hold variation info after fliping & evidence assignment
  $var_dba->dbc->do(qq{ DROP TABLE IF EXISTS variation_working});
  $var_dba->dbc->do(qq{ CREATE TABLE variation_working like variation });
  $var_dba->dbc->do(qq{ ALTER TABLE variation_working DROP COLUMN snp_id }); ## tmp column not in released schema
  $var_dba->dbc->do(qq{ ALTER TABLE variation_working DISABLE KEYS});

  ## temp table to hold variants with minor alleles not in the variation_feature allele string
  $var_dba->dbc->do(qq{ DROP TABLE IF EXISTS tmp_failed_minor_allele});
  $var_dba->dbc->do(qq{ CREATE TABLE tmp_failed_minor_allele  (
                        variation_id int(11) unsigned NOT NULL,                        
                        KEY variation_idx (variation_id)   )});

  ## table to hold variation feature info after fliping & ref allele assignment
  $var_dba->dbc->do(qq{ DROP TABLE IF EXISTS variation_feature_working});
  $var_dba->dbc->do(qq{ CREATE TABLE variation_feature_working like variation_feature });
  $var_dba->dbc->do(qq{ ALTER TABLE variation_feature_working DISABLE KEYS});
  ## new mysql version errors with empty not null columns
  ## switched to null allowable for import then back to non null here 
  $var_dba->dbc->do("alter table variation_feature_working modify column map_weight int not null");
 

  ## table to hold coded allele info after fliping 
  $var_dba->dbc->do(qq{ DROP TABLE IF EXISTS allele_working});
  $var_dba->dbc->do(qq{ CREATE TABLE allele_working (
                        allele_id int(11) NOT NULL AUTO_INCREMENT,
                        variation_id int(11) unsigned NOT NULL,
                        subsnp_id int(11) unsigned DEFAULT NULL,
                        allele_code_id int(11) unsigned NOT NULL,
                        population_id int(11) unsigned DEFAULT NULL,
                        frequency float unsigned DEFAULT NULL,
                        count int(11) unsigned DEFAULT NULL,
                        frequency_submitter_handle int(10) DEFAULT NULL,
                        PRIMARY KEY (allele_id),
                        KEY variation_idx (variation_id),
                        KEY subsnp_idx (subsnp_id),
                        KEY population_idx (population_id))
                        engine=MyISAM
                      });
  $var_dba->dbc->do(qq{ ALTER TABLE allele_working DISABLE KEYS});
  
  ## add intial values to allele code table
  $var_dba->dbc->do(qq{TRUNCATE allele_code});# empty if re-running to allow genotype_code population for basics
  $var_dba->dbc->do(qq{INSERT IGNORE INTO `allele_code` (`allele_code_id`,`allele`)
                        VALUES
                             (1, 'T'),
                             (2, 'A'),
                             (3, 'G'),
                             (4, 'C'),
                             (5, '-'),
                             (6, 'N')});


   ## table to hold failed variations     ## TEMP FOR DEBUG
  $var_dba->dbc->do(qq{DROP TABLE IF EXISTS failed_variation_working});
  $var_dba->dbc->do(qq{CREATE TABLE failed_variation_working like failed_variation});

  ## table to hold failed alleles        ## TEMP FOR DEBUG
  $var_dba->dbc->do(qq{DROP TABLE IF EXISTS failed_allele_working});
  $var_dba->dbc->do(qq{CREATE TABLE failed_allele_working like failed_allele});

  ## table to hold non-coded population_genotype info after flipping 
  $var_dba->dbc->do(qq{ DROP TABLE IF EXISTS MTMP_population_genotype_working });
  $var_dba->dbc->do(qq{ CREATE TABLE MTMP_population_genotype_working like population_genotype });

  ## new version of population_genotype table with coded genotypes
  $var_dba->dbc->do(qq{DROP TABLE IF EXISTS population_genotype_working});
  $var_dba->dbc->do(qq{CREATE TABLE population_genotype_working (
                        population_genotype_id int(10) unsigned NOT NULL AUTO_INCREMENT,
                        variation_id int(11) unsigned NOT NULL,
                        subsnp_id int(11) unsigned DEFAULT NULL,
                        genotype_code_id int(11) DEFAULT NULL,
                        frequency float DEFAULT NULL,
                        population_id int(10) unsigned DEFAULT NULL,
                        count int(10) unsigned DEFAULT NULL,
                        PRIMARY KEY (population_genotype_id),
                        KEY population_idx (population_id),
                        KEY variation_idx (variation_id),
                        KEY subsnp_idx (subsnp_id)
                        )engine=MyISAM
                      });


  ## temp table for assigning genotype ids
  $var_dba->dbc->do(qq{ DROP TABLE IF EXISTS genotype_code_tmp});
  $var_dba->dbc->do(qq{ CREATE TABLE genotype_code_tmp (
                        genotype_code_id int(11) unsigned NOT NULL AUTO_INCREMENT,
                        allele_1 varchar(30000) NOT NULL,
                        allele_2 varchar(30000) NOT NULL,
                        phased tinyint(2) unsigned DEFAULT NULL,
                        PRIMARY KEY ( genotype_code_id ),
                        UNIQUE KEY genotype_idx (allele_1 (490),allele_2(490),phased)
                       )});

  # add basic genotypes to genotype_code_tmp first - smaller numbers compress better - phased for 1KG
  $var_dba->dbc->do(qq{ INSERT IGNORE INTO genotype_code_tmp (allele_1, allele_2, phased) 
                        SELECT ac1.allele, ac2.allele, 1 
                        FROM allele_code ac1, allele_code ac2 
                        WHERE ac1.allele_code_id < 5 AND ac2.allele_code_id < 5
                        ORDER BY ac1.allele_code_id, ac1.allele_code_id + ac2.allele_code_id
                       });

  # add basic genotypes to genotype_code_tmp first - smaller numbers compress better - phasing unknown
  $var_dba->dbc->do(qq{ INSERT IGNORE INTO genotype_code_tmp (allele_1, allele_2, phased ) 
                        SELECT ac1.allele, ac2.allele, 0
                        FROM allele_code ac1, allele_code ac2 
                        ORDER BY ac1.allele_code_id, ac1.allele_code_id + ac2.allele_code_id
                       });

 $var_dba->dbc->do(qq{ CREATE TABLE IF NOT EXISTS maf(
                        snp_id           int(11),
                        allele           text,
                        freq             float,
                        count            int(11),
                        is_minor_allele  int(11) )
                      });

}

=head2 add_failed_variation_set

 enter variation_set_variation entry for failed variant set
 needed in VariantQC for variation_feature updating

=cut
sub add_failed_variation_set{

    my $var_dba = shift;

    my $fail_attrib_ext_sth  = $var_dba->dbc->prepare(qq[ select at.attrib_id
                                                          from attrib at, attrib_type att
                                                          where att.code = 'short_name' 
                                                          and att.attrib_type_id = at.attrib_type_id
                                                          and at.value = 'fail_all'
                                                        ]);
 
    my $variation_set_ext_sth  = $var_dba->dbc->prepare(qq[ select variation_set_id
                                                            from variation_set
                                                            where name = ?
                                                           ]);

    my $variation_set_ins_sth  = $var_dba->dbc->prepare(qq[ insert into variation_set
                                                            ( name, description, short_name_attrib_id)
                                                            values (?,?,?)
                                                          ]);

    ## check if already present
    $variation_set_ext_sth->execute('All failed variations')  || die "Failed to extract failed variant set id\n";
    my $failed_set_id = $variation_set_ext_sth->fetchall_arrayref();

    unless(defined $failed_set_id->[0]->[0]){
        ## no set entered - look up attrib for short name and enter set

        $fail_attrib_ext_sth->execute() || die "Failed to extract failed set attrib reasons\n";
        my $attrib = $fail_attrib_ext_sth->fetchall_arrayref();

        die "Exiting: Error - attribs not found. Load attribs then re-run\n" unless defined $attrib->[0]->[0] ;

        ## if attribs loaded, enter set
        $variation_set_ins_sth->execute( 'All failed variations',
                                         'Variations that have failed the Ensembl QC checks' ,
                                          $attrib->[0]->[0] )|| die "Failed to insert failed set\n"; 

        ## and pull out id to return
        $variation_set_ext_sth->execute('All failed variations')  || die "Failed to extract failed variant set id\n";
        $failed_set_id = $variation_set_ext_sth->fetchall_arrayref();      
    }

    return $failed_set_id->[0]->[0];
   
}

=head2 create_map_weight_table

create a look-up table for the number of times a variant maps to the reference

=cut
sub create_map_weight_table{
    
    my $core_dba  = shift;
    my $var_dba   = shift;

    #### is it better to use not in () or temp columns??
    my $ref_ext_sth = $core_dba->dbc->prepare(qq [ select sr.seq_region_id 
                                                   from seq_region sr, coord_system cs
                                                   where sr.coord_system_id = cs.coord_system_id  
                                                   and cs.rank = 1
                                                   and seq_region_id not in (select seq_region_id from seq_region_attrib where attrib_type_id = 16 )
                                                  ]);
    
     $ref_ext_sth->execute()|| die "Failed to extract ref/non ref status for seq_regions\n";
     my $is_ref = $ref_ext_sth->fetchall_arrayref();
    
    if ($var_dba->dbc->do(qq{show columns from seq_region like 'is_reference';}) != 1) {
      $var_dba->dbc->do(qq{alter table seq_region add column is_reference Tinyint(1) default 0});
    }
    my $sr_status_sth = $var_dba->dbc->prepare(qq[update seq_region set is_reference = 1 where seq_region_id  =?]);

    foreach my $srid (@{$is_ref}){
       $sr_status_sth->execute($srid->[0]) || die "ERROR updating seq regions\n";
    }


    #create a temporary table to store the map_weight, that will be deleted by the last process
    $var_dba->dbc->do(qq[ DROP TABLE IF EXISTS tmp_map_weight_working]);
    $var_dba->dbc->do(qq[ CREATE TABLE tmp_map_weight_working
                           SELECT variation_id, count(*) as count
                           FROM   variation_feature,seq_region
                            WHERE  variation_feature.seq_region_id = seq_region.seq_region_id
                           AND    seq_region.is_reference =1
                           GROUP BY variation_id]
              );

    $var_dba->dbc->do(qq{ALTER TABLE tmp_map_weight_working 
                         ADD UNIQUE INDEX variation_idx(variation_id)});


}


sub write_output {
    
    my $self = shift;

    ## Initial quick check for obvious failures in dbSNP import

    unless ($self->param('run_check_dbSNP_import') == 0){
        
      $self->dataflow_output_id($self->param('check_dbSNP_import'), 2);        
    }
   
    unless ($self->param('run_create_seqdb') == 0){
        
      $self->dataflow_output_id($self->param('create_seqdb'), 3);        
    }
   
    ## No map fails - larger bins used as check is very quick

    unless ($self->param('run_unmapped_var') == 0){

       my $unmapped_start_ids =  $self->param('unmapped_start_ids');
       $self->warning(scalar @{$unmapped_start_ids} .' unmapped_variant_qc jobs to do');  
       $self->dataflow_output_id($unmapped_start_ids, 5);     
    }

    ## Variant QC - bin start positions supplied

    unless ($self->param('run_variant_qc') == 0){

        my $qc_start_ids =  $self->param('qc_start_ids');
        $self->warning(scalar @{$qc_start_ids} .' variant_qc jobs to do');
        $self->dataflow_output_id($qc_start_ids, 4);  
    }

    ## Complement alleles in population_genotype for variants which are being flipped
 
    unless ($self->param('run_flip_population_genotype') == 0){

       $self->dataflow_output_id( $self->param('flip_population_genotype'), 6);
    }

      
    ## Migrate raw genotype data to new coded schema  
    
    unless ($self->param('run_update_population_genotype') == 0){

       $self->dataflow_output_id( $self->param('update_population_genotype'), 7);
    }
    ## bulk updates to statuses for special cases

    $self->dataflow_output_id($self->param('special_cases'), 8);

    ## run basic checks when everything is updated

    $self->dataflow_output_id($self->param('finish_variation_qc'), 9);
   
   
    
    return;
}

1;

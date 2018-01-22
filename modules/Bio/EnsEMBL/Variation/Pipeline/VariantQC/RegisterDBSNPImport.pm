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

=head1 DESCRIPTION

This module is used in the internal variation quality control process. 

=cut

## Runs at end of dbSNP import on simple hash:
##
## $data{dbSNP_name} = 'mouse_10090';
## $data{species}    = 'mus_musculus';
## $data{ensdb_name} = 'databasename3';
## $data{registry}   = Bio::EnsEMBL::Registry
## $data{ensdb_version} = 30;
## $data{assembly_version} = 38;
## $data{pwd}    = `pwd`;

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::RegisterDBSNPImport;

use strict;
use warnings;

use base qw(Exporter);
use Bio::EnsEMBL::Variation::Utils::QCUtils qw(count_rows);
use dbSNP::DBManager;

our @EXPORT_OK = qw(register);


## register new ensembl database from dbSNP import
## check key counts against previous
## store results and overall verdict

sub register{

    my $details = shift;  

    ## find current variation db record for this species & set to non-current
    $details->{prev_db}   = get_previous_ensvardb($details);

    ## add details of new variation db to production db
    $details->{new_db}    = register_new_ensvardb($details);
   
   
    $details->{result_ad} = $details->{registry}->get_adaptor('multi', 'intvar', 'result');

    open my $report, ">$details->{pwd}/QC_report.txt" ||die "Failed to open report file : $!\n";

    ## check row counts on tables
    my $totals_ok = check_table_counts($details, $report);

    ## check variant counts by seq_contig
    my $seqs_ok   = check_variant_by_sequence($details, $report);
  

    ## exit if there is nothing to compare to
    return unless defined $details->{prev_db} ;

    ### summary result from comparison
    my $is_passed = 1;
    $is_passed = 0 if $totals_ok ==0 || $seqs_ok ==0;

    my $result =  Bio::EnsEMBL::IntVar::Result->new_fast({ ensvardb     => $details->{new_db},
                                                           result_value => $is_passed,
                                                           result_type  => 'passed_comparison_to_previous',
                                                           parameter    => 'dbSNP',
                                                           adaptor      => $details->{result_ad}
                                                         });
    
    $details->{result_ad}->store($result);
}   

## Find previous ensembl variation database to extract counts for
##   - set as no longer current 
sub get_previous_ensvardb{

    my $details = shift;

    my $ensvardb_adaptor  = $details->{registry}->get_adaptor('multi', 'intvar', 'ensvardb');
    my $prev_db = $ensvardb_adaptor->fetch_current_by_species( $details->{species} );
    unless (defined $prev_db){
	warn "No previous database found to compare to for species $details->{species}\n"  ;
	return;
    }
    ## set last ensembl database to non-current    
    $ensvardb_adaptor->non_current($prev_db);

    ## hack - last database entered may have only had a transcript_variation update
    my $db_ext_sth = $ensvardb_adaptor->dbc->prepare(qq[ select max(ensvar_db.ensvardb_id) 
                                                         from ensvar_db, check_result 
                                                         where ensvar_db.species =?
                                                         and ensvar_db.ensvardb_id = check_result.ensvardb_id
                                                         and check_result.check_id = 1
                                                       ]);
    $db_ext_sth->execute( $details->{species});
    my $last_checked_id = $db_ext_sth->fetchall_arrayref();
    if (defined $last_checked_id->[0]->[0] && $last_checked_id->[0]->[0] ne $prev_db->dbID()){
        ## pick up last checked ensembl db object instead of current if current does not have variation counts
        my $checked_db = $ensvardb_adaptor->fetch_by_dbID( $last_checked_id->[0]->[0] );
        print "using ".$checked_db->name() . " as previous ensembl db to check counts against (latest load, not most current release)\n";
        return  $checked_db;
    }
    else{
        print "using ".$prev_db->name() . " as previous ensembl db to check counts against\n";
        return $prev_db;
    }
}

## store required information for new ensembl variation database

sub register_new_ensvardb{

    my $details = shift;

    ## find  dbSNP database record to link to
    my $dbSNP_db_dba =  $details->{registry}->get_adaptor('multi', 'intvar', 'dbSNPdb');
    my $dbSNP_db     =  $dbSNP_db_dba->fetch_current_by_name($details->{dbSNP_name} );
    unless (defined $dbSNP_db){
        warn "Not adding dbSNP database link - no database of this name $details->{dbSNP_name} recorded\n";
    }

    ## check if ensembl db name novel
    my $ensvardb_dba =  $details->{registry}->get_adaptor('multi', 'intvar', 'EnsVardb');
    my $preexisting  =  $ensvardb_dba->fetch_by_name( $details->{ensdb_name} ); 

    if (defined $preexisting){
        warn "Not updating - ensembl database of this name  $details->{ensdb_name} already recorded\n";
        exit;
    }
   
    warn "Recording details for ensembl database  $details->{ensdb_name} \n";

    ## enter new db 
    my $ensvar_db = Bio::EnsEMBL::IntVar::EnsVardb->new_fast({ name            => $details->{ensdb_name} ,
                                                              dbsnpdb          => $dbSNP_db,
                                                              species          => $details->{species},
                                                              version          => $details->{ensdb_version},
                                                              genome_reference => $details->{assembly_version},
                                                              adaptor          => $ensvardb_dba,
                                                              status_desc      => 'dbSNP_imported'
                                                             });


    $ensvardb_dba->store($ensvar_db);

    return $ensvar_db;

}

## check row counts on tables not expected to decrease in size 

sub check_table_counts{

    my $details = shift;
    my $report = shift;

    my $var_dba = $details->{registry}->get_DBAdaptor($details->{species}, 'variation');

    my %tables_to_check  = ( 'variation'           =>  'dbSNP_variants',
                             'allele'              =>  'dbSNP_alleles',
                             'variation_feature'   =>  'dbSNP_variation_features', 
                             'population_genotype' =>  'dbSNP_population_genotype',
        );

    my $ok = 1;
    foreach my $table (keys  %tables_to_check){

        my $new_count = count_rows($var_dba , $table );

        ## can only check against previously if data held
	if (defined $details->{prev_db}){
	    my @prev_results = @{$details->{result_ad}->fetch_by_db_result_type($details->{prev_db},$tables_to_check{$table} )} ;
	    
	    if( defined $prev_results[0]){
		my $previous_row_count = $prev_results[0]->result_value();
		if($new_count  >= $previous_row_count){
		    print $report "OK: $table\t previously $previous_row_count now $new_count \n";
		}
		else{
		    print $report "ERROR: $table\t previously $previous_row_count now $new_count \n";
		    $ok =0;
		}
            }

        }
        ## store new result
        my $result =  Bio::EnsEMBL::IntVar::Result->new_fast({ ensvardb     => $details->{new_db},
                                                               result_value => $new_count,
                                                               result_type  => $tables_to_check{$table},
                                                               adaptor      => $details->{result_ad}
                                                               });
    
        $details->{result_ad}->store($result);

    }
    return $ok;
}

## check variant counts by sequence to spot missing chromosomes

sub check_variant_by_sequence{

    my $details = shift;
    my $report  = shift;

    ## find counts by sequence for new database
    my $new_seq_count = new_seq_count($details);

    ## find counts by sequence for previous database
    my $old_seq_count = old_seq_count($details) if defined $details->{prev_db};


    my $ok = 1;
    ### check no old sequences with variants are missing
    foreach my $seq (keys %{$old_seq_count}){

       if (defined $new_seq_count->{$seq} ){
           if($new_seq_count->{ $seq } >= $old_seq_count->{$seq} ){
               print $report "OK\tsequence $seq :\tvariation count previously $old_seq_count->{$seq} and now $new_seq_count->{$seq}\n";
           }
           else{
               print $report "ERROR\tsequence $seq :\tvariation count  previously $old_seq_count->{$seq} and now $new_seq_count->{$seq}\n";
               $ok = 0;
           }
       }
       else{
           print $report "ERROR\tsequence $seq : no variation count previously $old_seq_count->{$seq} \n";
           $ok = 0;
       }
    }

    ## store new variant counts 
    foreach my $seq (keys %{$new_seq_count}){

       my $result =  Bio::EnsEMBL::IntVar::Result->new_fast({ ensvardb     => $details->{new_db},
                                                              result_value => $new_seq_count->{$seq},
                                                              result_type  => 'dbSNP_variants_by_seq',
                                                              parameter    => $seq,
                                                              adaptor      => $details->{result_ad}
                                                            });
    
      $details->{result_ad}->store($result);
    }

    ### check for new sequences with variants
    foreach my $seq (keys %{$new_seq_count}){

       unless (defined $old_seq_count->{$seq} ){
           print $report "WARNING\t sequence $seq : showing $new_seq_count->{$seq} variants now, non previously\n";

       }
    }
    return $ok;

}

sub old_seq_count{

    my $details = shift;

    my %old_seq_count;
  
    my @prev_results = @{$details->{result_ad}->fetch_by_db_result_type($details->{prev_db}, 'dbSNP_variants_by_seq')};

    foreach my $res (@prev_results){
       $old_seq_count{ $res->parameter() } = $res->result_value();
    }

    return \%old_seq_count;
}

sub new_seq_count{

    my $details = shift;

    my %new_seq_count;

    my $var_dba = $details->{registry}->get_DBAdaptor($details->{species}, 'variation');
    my $data_ext_sth = $var_dba->dbc->prepare(qq[ select seq_region.name, count(*) 
                                                  from seq_region, variation_feature
                                                  where variation_feature.seq_region_id = seq_region.seq_region_id 
                                                  group by seq_region.name]);

    $data_ext_sth->execute()||die "Failed to count variation per sequence\n";

    my $count = $data_ext_sth->fetchall_arrayref();

    foreach my $l(@{$count}){

       $new_seq_count{$l->[0]} = $l->[1];
    }

    return (\%new_seq_count);

}



1;

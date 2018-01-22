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


=head1 CheckTranscriptVariation

This module runs at the end of the Transcript Effect pipeline and produces 
a summary report to check the results look reasonable.


=cut

package Bio::EnsEMBL::Variation::Pipeline::CheckTranscriptVariation;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub run {
    my $self = shift;
    
    ## retrieve old results from production db for comparison if available
    my $previous_counts = $self->get_old_results();  

    ## calculate new counts from new transcript_variation table
    my $new_counts      = $self->get_new_results();

    ## write report in working directory 
    $self->report_results($new_counts, $previous_counts);

    ## updated production db for later use
    $self->update_internal_db($new_counts);

}

=head2 get_new_results

Run all the counting SQL on the new database

=cut
sub get_new_results{

    my $self     = shift;
    my $previous = shift;  ## previous status only available if production db connection details supplied
  
    my $var_dba = $self->get_species_adaptor('variation');   
    my $dbc     = $var_dba->dbc;

    my %new;  ## hold some of the new counts to store
 
    ## counts on consequences relevant to all species

    my $tv_count_st       = qq[ select count(*) from transcript_variation];    
    my $tot_missense_st   = qq[ select count(*) from transcript_variation where consequence_types like '%missense_variant%' ];
    my $tot_synonymous_st = qq[ select count(*) from transcript_variation where consequence_types like '%synonymous_variant%' ];
    my $missing_hgvs_st   = qq[ select count(*) from transcript_variation where consequence_types like '%missense_var%' and hgvs_protein is null];
    my $missing_pep_st    = qq[ select count(*) from transcript_variation where consequence_types like '%missense_var%' and pep_allele_string is null];
    my $missing_con_st    = qq[ select count(*) from transcript_variation where consequence_types=''];

    
    ## counts on protein impact only relevant to some species

    my $pred_count_st = qq[ select attrib.value,count(*) 
                            from attrib, protein_function_predictions 
                            where protein_function_predictions.analysis_attrib_id = attrib.attrib_id 
                            group by attrib.value ];    
    
    my $sift_count_st = qq[ select count(*) from transcript_variation where sift_score is not null ];
    

    ## counts on polyphen protein impact only relevant to human

    my $pp_count_st       = qq[ select count(*) from transcript_variation where polyphen_score is not null ];
    
    my $agree_count_st    = qq[ select count(*) from transcript_variation 
                                where polyphen_prediction like '%damaging' and sift_prediction = 'deleterious'];
    
    my $disagree_count_st = qq[ select count(*) from transcript_variation 
                                where polyphen_prediction = 'benign' and sift_prediction = 'deleterious'];
    


    ## General results

    $new{transcript_variation_count}  = count_results($dbc, $tv_count_st );
    unless($new{transcript_variation_count} > 0){  ## report & die if total failure
        $self->warning('ERROR: no transcript variation results found');
        die;
    }

    $new{missense_variant_count}   = count_results($dbc, $tot_missense_st );
    $new{synonymous_variant_count} = count_results($dbc, $tot_synonymous_st );

    $new{missing_hgvs_protein}     = count_results($dbc, $missing_hgvs_st);
    $new{missing_peptide_alleles}  = count_results($dbc, $missing_pep_st);
    $new{missing_consequences}     = count_results($dbc, $missing_con_st); 

    

    ## Protein impact checks

    ### check species has Sift data
    my $pred_matrix_check_sth = $dbc->prepare(qq[show tables like 'protein_function_predictions']);
    $pred_matrix_check_sth->execute();
    my $present = $pred_matrix_check_sth->fetchall_arrayref();
    return unless defined  $present->[0]->[0];

    my $total_pred = count_results($dbc, $pred_count_st );    

    if( defined $total_pred->{sift} ){

        $new{Sift_predictions}  = count_results($dbc, $sift_count_st);
    }
    
    ### Polyphen available for human only

    if(defined $total_pred->{polyphen_humvar} ||$total_pred->{polyphen_humdiv} ){
        
        $new{Polyphen_predictions}     = count_results($dbc, $pp_count_st );
        
        ## check whether sift deleterious calls also look damaging with polyphen 
        $new{Sift_Polyphen_agree}    = count_results($dbc, $agree_count_st );
        $new{Sift_Polyphen_disagree} = count_results($dbc, $disagree_count_st);
    }

    return \%new;

}


=head2 report_results

Print a report in the working directory showing current and previous data

=cut
sub report_results{

    my $self     = shift;
    my $new      = shift;
    my $previous = shift;


    my $dir =  $self->required_param('pipeline_dir');
    open my $report, ">$dir/QC_report.txt"||die "Failed to open report file for summary info :$!\n";


    my $percent_missense        = format_percent( $new->{missense_variant_count},   $new->{transcript_variation_count} ) ;
    my $percent_synonymous      = format_percent( $new->{synonymous_variant_count}, $new->{transcript_variation_count} ) ;

    my $prev_percent_missense   = format_percent( $previous->{missense_variant_count},   $previous->{transcript_variation_count} ) ;
    my $prev_percent_synonymous = format_percent( $previous->{synonymous_variant_count}, $previous->{transcript_variation_count} ) ;


    my $text_out = "\nSummary of results from TranscriptVariation\n\n";
    $text_out.= "$new->{transcript_variation_count} transcript_variation entries";
    $text_out.= " (previously $previous->{transcript_variation_count} )" if defined  $previous->{transcript_variation_count} ;
    
    $text_out.=  "\n\n$new->{missense_variant_count} ($percent_missense%) missense variants";
    $text_out.= " (previously $prev_percent_missense%)" if defined  $prev_percent_missense ;

    $text_out.=  "\n$new->{synonymous_variant_count} ($percent_synonymous%) synonymous variants";
    $text_out.= " (previously $prev_percent_synonymous%)" if defined $prev_percent_synonymous;

    $text_out.=  "\n\n$new->{missing_hgvs_protein} missense variants don't have HGVS annotation\n";
    $text_out.=  "$new->{missing_peptide_alleles} missense variants don't have peptide allele strings\n";
    $text_out.=  "$new->{missing_consequences} variants don't have consequences\n";

    print $report $text_out;


    ### basic report ends here - only a few species have Sift predictions
    return unless defined $new->{Sift_predictions};

    my $percent_sift      = format_percent( $new->{Sift_predictions}, $new->{missense_variant_count} ) ;
    my $prev_percent_sift = format_percent( $previous->{Sift_predictions}, $previous->{missense_variant_count} ) ;
    
    print $report "\n$percent_sift% of misssense entries have Sift results (previously $prev_percent_sift%)\n";
    
 
    ### only human has Polyphen data
    return unless defined $new->{Polyphen_predictions};

    my $percent_poly      = format_percent( $new->{Polyphen_predictions}, $new->{missense_variant_count}) ;
    my $prev_percent_poly = format_percent( $previous->{Polyphen_predictions}, $previous->{missense_variant_count}) ;

    print $report "$percent_poly% of misssense entries have polyphen results (previously $prev_percent_poly%)\n";
        
    ## check whether sift deleterious calls also look damaging with polyphen 
    my $sift_del          = $new->{Sift_Polyphen_agree} + $new->{Sift_Polyphen_disagree};        
    my $agree_per         = format_percent( $new->{Sift_Polyphen_agree}, $sift_del ) ;
    my $disagree_per      = format_percent( $new->{Sift_Polyphen_disagree}, $sift_del);

    my $prev_sift_del     = $previous->{Sift_Polyphen_agree} + $previous->{Sift_Polyphen_disagree};         
    my $prev_agree_per    = format_percent( $previous->{Sift_Polyphen_agree}, $prev_sift_del ) ;
    my $prev_disagree_per = format_percent( $previous->{Sift_Polyphen_disagree}, $prev_sift_del);
    
            
    print $report "\n$sift_del combinations called deleterious by Sift also called by Polyphen: \n";
    print $report "$agree_per% Sift deleterious calls also called by Polyphen rated damaging (previously $prev_agree_per% from $prev_sift_del)\n";
    print $report "$disagree_per% Sift deleterious calls also called by Polyphen rated benign (previously $prev_disagree_per%  from $prev_sift_del)\n\n";  


}

=head2 get_old_results

Check internal production database for previous protein impact information
for this species 

=cut
sub get_old_results{

    my  $self = shift;

    my $int_dba ;
    eval{ $int_dba = $self->get_adaptor('multi', 'intvar');};

    unless (defined $int_dba){
        $self->warning('No internal database connection found to write status '); 
        return;
    }
    
    my %previous_result;

    my $result_adaptor = $int_dba->get_ResultAdaptor();
    my $res = $result_adaptor->fetch_all_current_by_species($self->required_param('species') );
   
    foreach my $result (@{$res}){

        $previous_result{ $result->result_type() } = $result->result_value();
    }
    return \%previous_result;

}


=head2 update_internal_db

Update internal production database with new statuses   

=cut
sub update_internal_db{

    my $self             = shift;
    my $new_counts       = shift;
  
    my $int_dba ;
    eval{ $int_dba = $self->get_adaptor('multi', 'intvar');};

    unless (defined $int_dba){
        $self->warning('No internal database connection found to write status '); 
        return;
    }

    my $var_dba = $self->get_species_adaptor('variation');
    my $ensdb_name = $var_dba->dbc->dbname;

    my $ensvardb_dba  =  $int_dba->get_EnsVardbAdaptor();
    my $result_dba    =  $int_dba->get_ResultAdaptor();

    my $ensdb = $ensvardb_dba->fetch_by_name($ensdb_name);
    unless (defined $ensdb ){

        ## add new db  *** NOT setting last to non-current - fix this ***
        my @b = split/\_/,$ensdb_name;
        pop @b;
        my $ens_version = pop @b;

        $ensdb = Bio::EnsEMBL::IntVar::EnsVardb->new_fast({ name             => $ensdb_name,
                                                            species          => $self->required_param('species'),
                                                            version          => $ens_version,
                                                            status_desc      => 'Created'
                                                          });
                
        $ensvardb_dba->store( $ensdb );
    }
    $ensvardb_dba->update_status( $ensdb, 'variation_consequence_run' );

        
    foreach my $type ( keys %{$new_counts}){

       ## set any previous results to non current for this check type and species
        $result_dba->set_non_current_by_species_and_type($self->required_param('species') , $type);
      
        my $result =  Bio::EnsEMBL::IntVar::Result->new_fast({ ensvardb     => $ensdb,
                                                               result_value => $new_counts->{$type},
                                                               result_type  => $type,
                                                               adaptor      => $result_dba
                                                             });
        
        $result_dba->store($result);
    }
    
}

=head2 count_results

 This takes SQL statements to count the rows in a table or
 count rows grouped by an attribute in the table.
 It returns either the total number of rows in the table or 
 a hash of attribute => row count depending on input.
=cut
sub count_results{

    my $dbc = shift;
    my $st  = shift;

    my $sth = $dbc->prepare($st);
    $sth->execute();
    my $dat = $sth->fetchall_arrayref();

    if(defined $dat->[0]->[1]){
        my %count;
        foreach my $l (@{$dat}){
            $count{$l->[0]} = $l->[1];
        }
        return \%count;
    }
    else{
        return $dat->[0]->[0];
    }
}


sub format_percent{

    my $part  = shift;
    my $whole = shift;

    return unless (defined $whole && $whole > 0);

    return  substr((100 * $part / $whole ),0,5) ;
}


1;


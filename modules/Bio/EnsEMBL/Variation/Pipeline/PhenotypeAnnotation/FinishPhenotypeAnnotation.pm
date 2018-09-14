=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

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


=head1 FinishPhenotypeAnnotation

This module runs at the end of the Phenotype Annotation pipeline and produces 
a summary report to check the results look reasonable.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::FinishPhenotypeAnnotation;

use strict;
use warnings;
use Data::Dumper; #TODO: remove once done testing

use base qw(Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation);

sub run {
    my $self = shift;
    
    ## retrieve old results from production db for comparison if available
    my $previous_counts = $self->get_old_results();  

    ## calculate new counts from new phenotype_feature table
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

    my $source = $self->param('source');
    
    my %new;  ## hold some of the new counts to store
 
    ## counts on phenotypes relevant to all species
    my $phenotype_count_st                    = qq[ select count(*) from phenotype];
    my $phenotype_feature_count_st            = qq[ select count(*) from phenotype_feature ];
    my $phenotype_feature_attrib_count_st     = qq[ select count(*) from phenotype_feature_attrib ];
    my $phenotype_ontology_accession_count_st = qq[ select count(*) from phenotype_ontology_accession ];  

    ## counts based on type and source
    my $phenotype_feature_grouped_count_st = qq[select s.name, pf.type, count(*)
                                            from phenotype_feature pf, source s
                                            where pf.source_id = s.source_id
                                            group by s.name, pf.type ];

    ## General results
    $new{phenotype_count}  = count_results($dbc, $phenotype_count_st );
    unless($new{phenotype_count} > 0){  ## report & die if total failure
        $self->warning('ERROR: no phenotypes found');
        die;
    }
    $new{phenotype_feature_count}   = count_results($dbc, $phenotype_feature_count_st );
    unless($new{phenotype_feature_count} > 0){  ## report & die if total failure
        $self->warning('ERROR: no phenotype feature results found');
        die;
    }

    $new{phenotype_feature_attrib_count} = count_results($dbc, $phenotype_feature_attrib_count_st );
    $new{phenotype_ontology_accession_count} = count_results($dbc, $phenotype_ontology_accession_count_st );

    # get grouped counts
    my $sth = $dbc->prepare($phenotype_feature_grouped_count_st);
    $sth->execute();
    my $dat = $sth->fetchall_arrayref();
    foreach my $l (@{$dat}){
      $new{phenotype_feature_count_details}{$l->[0]."_".$l->[1]} = $l->[2];
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

    my $dir =  $self->required_param('workdir');
    open my $report, ">$dir/REPORT_import.txt"||die "Failed to open report file for summary info :$!\n";

    my $text_out = "\nSummary of results from CheckPhenotypeAnnotation\n\n";
    $text_out.= "$new->{phenotype_count} phenotype entries";
    $text_out.= " (previously $previous->{phenotype_count})" if defined  $previous->{phenotype_count} ;
    $text_out.= "\n";
    
    $text_out.= "$new->{phenotype_ontology_accession_count} phenotype_ontology_accession entries";
    $text_out.= " (previously $previous->{phenotype_ontology_accession_count})" if defined  $previous->{phenotype_ontology_accession_count} ;
    $text_out.= "\n";

    $text_out.= "$new->{phenotype_feature_count} phenotype_feature entries";
    $text_out.= " (previously $previous->{phenotype_feature_count})" if defined  $previous->{phenotype_feature_count} ;
    $text_out.= "\n";

    $text_out.= "$new->{phenotype_feature_attrib_count} phenotype_feature_attrib entries";
    $text_out.= " (previously $previous->{phenotype_feature_attrib_count})" if defined  $previous->{phenotype_feature_attrib_count} ;
    $text_out.= "\n";
    
    foreach my $source_type ( keys %{$new->{phenotype_feature_count_details}}){
      $text_out.= "$new->{phenotype_feature_count_details}{$source_type} $source_type phenotype_feature entries";
      $text_out.= " (previously $previous->{phenotype_feature_count_details}{$source_type} )" if defined  $previous->{phenotype_feature_count_details}{$source_type} ;
      $text_out.= "\n";
    }

    print $report $text_out;
}

=head2 get_old_results

Check internal production database for previous phenotype annotation import information
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

      if ($result->parameter()){
        $previous_result{ $result->result_type()."_details" }{$result->parameter()} = $result->result_value();
      } else {
        $previous_result{ $result->result_type() } = $result->result_value();
      }
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
    $ensvardb_dba->update_status( $ensdb, 'phenotype_annotation_run' );

    foreach my $type ( keys %{$new_counts}){
      
      if ($type eq 'phenotype_feature_count_details'){
        ## set any previous results to non current for this check type and species
        $result_dba->set_non_current_by_species_and_type($self->required_param('species') , 'phenotype_feature_count');
    
        my $details = $new_counts->{$type};
        foreach my $source_type ( keys %{$details}){
          my $result =  Bio::EnsEMBL::IntVar::Result->new_fast({ ensvardb     => $ensdb,
                                                                 result_value => $details->{$source_type},
                                                                 result_type  => 'phenotype_feature_count',
                                                                 parameter => $source_type,
                                                                 adaptor      => $result_dba
                                                               });

          $result_dba->store($result);
        }
      } else {
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


1;


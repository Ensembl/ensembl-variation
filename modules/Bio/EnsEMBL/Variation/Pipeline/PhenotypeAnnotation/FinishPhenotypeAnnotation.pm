=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

This module runs at the end of the phenotype annotation import pipeline and produces
a summary report to check the results against the previous run counts.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::FinishPhenotypeAnnotation;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation);

sub fetch_input {
  my $self = shift;

  $self->variation_db_adaptor($self->get_species_adaptor('variation'));

}

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

  $self->update_meta();
}


=head2 get_new_results

  Arg [1]    : hashref $previous (optional)
               The previous status counts of the analysis.
  Example    : $new_counts = $obj->get_new_results()
  Description: Run all the counting SQL on the new database.
  Returntype : hashref of new counts
  Exceptions : none

=cut

sub get_new_results{
  my ($self, $previous) = @_; ## previous status only available if production db connection details supplied

  my $dbc     = $self->variation_db_adaptor->dbc;
  my $source = $self->param('source');

  my %new;  ## hold some of the new counts to store

  ## counts based on type and source
  my $phenotype_feature_grouped_count_st = qq[select s.name, pf.type, count(*)
                                          from phenotype_feature pf, source s
                                          where pf.source_id = s.source_id
                                          group by s.name, pf.type ];

  ## counts on phenotypes relevant to all species
  my @tables = ('phenotype', 'phenotype_feature', 'phenotype_feature_attrib', 'phenotype_ontology_accession');
  foreach my $table (@tables) {
    my $count_st = qq[ select count(*) from $table];
    $new{"$table\_count"}  = count_results($dbc, $count_st);
    unless($new{"$table\_count"} > 0){  ## report & die if total failure
        $self->warning("WARNING: no entries found in the table $table\n");
        die unless ($table eq 'phenotype_feature_attrib' ||
                    $table eq 'phenotype_ontology_accession');
    }
  }

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

  Arg [1]    : hashref $new
               The new counts of the analysis.
  Arg [2]    : hashref $previous
               The previous counts of the analysis.
  Example    : $obj->report_results($new, $previous)
  Description: Print a report in the working directory showing current and previous data.
  Returntype : none
  Exceptions : none

=cut

sub report_results{
  my ($self, $new, $previous) = @_;


  my $dir =  $self->required_param('workdir');
  my $report;
  open ($report, ">$dir/REPORT_import.txt") || die ("Failed to open report file for summary info :$!\n");

  my $text_out = "\nSummary of results from CheckPhenotypeAnnotation\n\n";

  my @tables = ('phenotype', 'phenotype_feature', 'phenotype_feature_attrib', 'phenotype_ontology_accession');
  foreach my $table (@tables) {
    $text_out.= $new->{"$table\_count"}." $table entries";
    $text_out.= " (previously ".$previous->{"$table\_count"}.")" if defined  $previous->{"$table\_count"} ;
    $text_out.= "\n";
  }

  foreach my $source_type ( keys %{$new->{phenotype_feature_count_details}}){
    $text_out.= "$new->{phenotype_feature_count_details}{$source_type} $source_type phenotype_feature entries";
    $text_out.= " (previously $previous->{phenotype_feature_count_details}{$source_type})" if defined  $previous->{phenotype_feature_count_details}{$source_type} ;
    $text_out.= "\n";
  }

  print $report $text_out;
  close $report;
}


=head2 get_old_results

  Example    : $obj->get_old_results()
  Description: Check internal production database for previous phenotype annotation import information
               for this species.
  Returntype : hashref of previous counts
  Exceptions : none

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

  Example    : $obj->update_internal_db($new)
  Description: Update internal production database with new statuses.
  Returntype : none
  Exceptions : none

=cut

sub update_internal_db{
  my ($self, $new_counts) = @_;

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

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBConnection $dbc
               The new variation database connection
  Arg [2]    : string $st
               The SQL statement to be run.
  Example    : $obj->count_results($dbc, $st)
  Description: Takes SQL statements to count the rows in a table or count rows grouped
               by an attribute in the table. It returns either the total number of rows in the
               table or a hash of attribute => row count depending on input.
  Returntype : interger or hashref
  Exceptions : none

=cut

sub count_results{
  my ($dbc, $st) = @_;

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


=head2 update_meta

  Example    : $obj->update_meta()
  Description: Store the pipeline name, date and imported source in the species meta table.
               key=PhenotypeAnnotation_run_date_<source_name> value=run_date
  Returntype : none
  Exceptions : none

=cut

sub update_meta{
  my $self = shift;

  my $source_info = $self->param("source");
  my $var_dbh = $self->variation_db_adaptor->dbc->db_handle;

  my $update_meta_sth = $var_dbh->prepare(qq[ insert ignore into meta
                                              ( meta_key, meta_value) values (?,?)
                                            ]);

  $update_meta_sth->execute('PhenotypeAnnotation_run_date_'.$source_info->{source_name}, $self->run_date() );

}

1;


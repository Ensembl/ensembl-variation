=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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
use POSIX qw(strftime);
use File::Path qw(make_path);
use Data::Dumper;

use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(OMIA ANIMALSET MOUSE NONE SPECIES);

use base qw(Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation);

sub fetch_input {
  my $self = shift;

  my $workdir = $self->param('pipeline_dir')."/FinalChecks";
  unless (-d $workdir) {
    my $err;
    make_path($workdir, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }

  open(my $logFH, ">>", $workdir."/REPORT_import_".$self->required_param('species').".txt") || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  
  open(my $pipelogFH, ">", $self->param('pipeline_dir')."/".'log_import_debug_pipe_finish') || die ("Failed to open file: $!\n");
  $self->pipelogFH($pipelogFH);
}

sub run {
  my $self = shift;

  ## before finishing, meta_coord needs updating
  $self->update_meta_coord();

  ## retrieve old results from production db for comparison if available
  my $previous_counts = $self->get_old_results();

  ## calculate new counts from new phenotype_feature table
  my $new_counts      = $self->get_new_results();

  ## write report in working directory
  $self->report_results($new_counts, $previous_counts);

  ## updated production db for later use
  $self->update_internal_db($new_counts);
   
  $self->clean_phenotype_tables();
  
  $self->process_output_ids();
  
  close($self->logFH) if defined $self->logFH ;
  close($self->pipelogFH) if defined $self->pipelogFH ;
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

  my $time = strftime("%Y-%m-%d %H:%M:%S", localtime);
  $self->print_logFH("Running time: $time\n");

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

  $self->print_logFH($text_out);
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

  my $species = $self->required_param('species');
  my $var_dba = $self->get_species_adaptor('variation');
  my $ensdb_name = $var_dba->dbc->dbname;
  my $genome_assembly = $self->get_assembly;

  my $ensvardb_dba  =  $int_dba->get_EnsVardbAdaptor();
  my $result_dba    =  $int_dba->get_ResultAdaptor();

  my $ensdb = $ensvardb_dba->fetch_by_name($ensdb_name);
  unless (defined $ensdb ){

    ## add new db
    my @b = split/\_/,$ensdb_name;
    pop @b;
    my $ens_version = pop @b;

    $ensdb = Bio::EnsEMBL::IntVar::EnsVardb->new_fast({ 
      name        => $ensdb_name,
      species     => $self->required_param('species'),
      version     => $ens_version,
      status_desc => 'Created'
    });
    $ensdb->genome_reference($genome_assembly) if defined $genome_assembly;
    $ensvardb_dba->store( $ensdb );
  }
  $ensvardb_dba->update_status( $ensdb, 'phenotype_annotation_run' );

  my @counts = ('phenotype_count', 'phenotype_feature_count',
                'phenotype_feature_attrib_count', 'phenotype_ontology_accession_count');
  foreach my $type (@counts) {
    if (!defined $new_counts->{$type}){
      $self->print_logFH("$type count not found in current db!\n");
      next;
    }
    ## set any previous results to non current for this check type and species
    if ($species eq 'homo_sapiens') {
      $result_dba->set_non_current_by_species_assembly_and_type($species, $genome_assembly, $type);
    } else {
      $result_dba->set_non_current_by_species_and_type($species , $type);
    }
    my $result =  Bio::EnsEMBL::IntVar::Result->new_fast({ ensvardb     => $ensdb,
                                                           result_value => $new_counts->{$type},
                                                           result_type  => $type,
                                                           parameter    => 'All',
                                                           adaptor      => $result_dba
                                                         });
    $result_dba->store($result);
  }

  # save grouped counts
  # previous release phenotype_feature_count details already set to not current in above loop
  my $type = 'phenotype_feature_count_details';
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
}



=head2 process_output_ids

  Example    : $obj->process_output_ids($new)
  Description: Process output_ids to create job for next stage (datacheck).
  Returntype : none
  Exceptions : none

=cut

sub process_output_ids {
  my $self = shift;
  
  my $run_type = $self->required_param('run_type');
  my %import_species = SPECIES;
  
  $self->warning("run type $run_type");
  unless ($run_type eq NONE) {
    my ($type, @dc_params);
    $type = ($run_type eq ANIMALSET) ? 'OMIA' : $run_type;
    $type = ($run_type eq MOUSE) ? 'IMPC' : $run_type;
    
    map {
      my $adaptor = $self->get_adaptor($_, 'variation');
      my $dbname = $adaptor->dbc()->dbname();
      
      my $old_server_uri = $self->param('old_server_uri'); 
      unless ( $old_server_uri ){
        my $host_domain = $adaptor->dbc()->hostname();
        my ($host) = ( split(/\./, $host_domain) );
        my $port = $adaptor->dbc()->port();
        my $user = $adaptor->dbc()->user();
        my $pass = $adaptor->dbc()->pass();
        my $old_release_number = ( $self->param('ensembl_release') - 1 );
        
        $old_server_uri ||= sprintf("mysql://%s:%s@%s:%s/%s", $user, $pass, $host, $port, $old_release_number);
      }
      
      push @dc_params, { dbname => $dbname, old_server_uri => [ $old_server_uri ] }; 
    
    } @{$import_species{$type}};
    
    $self->param('output_ids', \@dc_params );
    
    $self->print_pipelogFH("Setting up for datacheck import: ". join(", ",@{$import_species{$type}}). "\n") if $self->param('debug_mode') ;
  }
}



sub write_output {
    my $self = shift;

    $self->dataflow_output_id($self->param('output_ids'), 1);
}

1;


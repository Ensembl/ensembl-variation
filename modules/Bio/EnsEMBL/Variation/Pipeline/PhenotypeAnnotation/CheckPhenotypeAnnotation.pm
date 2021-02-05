=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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


=head1 CheckPhenotypeAnnotation

This module runs checks and produces reports after the phenotype import stage.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation;

use strict;
use warnings;
use POSIX qw(strftime);
use File::Path qw(make_path);
use Data::Dumper;

use base qw(Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation);
use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(SPECIES MOUSE MGI OMIA HUMAN HUMAN_VAR HUMAN_GENE CGC EGA ANIMALSET ANIMALQTL RGD ZFIN);

my $source;
my $workdir;
my $species;
my $report;
my $count_ok = 1;

my %groups_end_source = (
  HUMAN      => CGC,
  HUMAN_GENE => CGC,
  HUMAN_VAR  => EGA,
  ANIMALSET  => ANIMALQTL,
  MOUSE      => MGI,
);

sub fetch_input {
  my $self = shift;

  $species = $self->required_param('species');
  $source = $self->param('source');
  $workdir = $self->param('workdir');
  my $run_type =  $self->required_param('run_type');

  my ($logName, $logPipeName);
  if (defined $source) {
    $workdir ||= $self->required_param('pipeline_dir')."/".$source->{source_name}."/".$self->required_param('species');
    $logName = "REPORT_QC_".$source->{source_name}.".txt";
    $logPipeName = "log_import_debug_pipe_".$source->{source_name}."_".$self->param('species');
  } else {
    $workdir ||= $self->required_param('pipeline_dir')."/FinalChecks";
    unless (-d $workdir) {
      my $err;
      make_path($workdir, {error => \$err});
      die "make_path failed: ".Dumper($err) if $err && @$err;
    }
    $logName = "REPORT_QC_". $species .".txt";
    $logPipeName = "log_import_debug_pipe_".$self->param('species');
  }

  open(my $logFH, ">>", $workdir."/".$logName) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">>", $workdir."/".$logPipeName) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->pipelogFH($pipelogFH);

  $self->param('output_ids', [{run_type => $run_type, species => $species}]);

}

sub run {
  my $self = shift;

  my $dbh = $self->variation_db_adaptor->dbc;
  my $time = strftime("%Y-%m-%d %H:%M:%S", localtime);
  $self->print_logFH("Running time: $time\n");

  $self->print_logFH("\nRunning checks on phenotype, phenotype_feature".
                     "phenotype_feature_attrib and phenotype_ontology_accession data\n");
  $self->check_phenotype_description($dbh);
  $self->print_logFH("\n");
  $self->check_fk($dbh);
  $self->check_source();
  $self->update_meta();

}



sub write_output {
  my $self = shift;

  my $run_type = $self->param('run_type');

  #if these is a decrease in the number of entries, then stop the flow
  if (!$count_ok){
    $self->print_logFH("ERROR: check counts failed! No futher jobs will be triggerd!\n".
                       "PLEASE check import and redo import if needed!");
    close($self->logFH) if defined $self->logFH ;
    close($self->pipelogFH) if defined $self->pipelogFH ;
    die "Failed check_counts: decrease in one of the tables: phenotype, phenotype_feature, phenotype_feature_attrib, phenotype_ontology_accession a get DBA for $species \n";
  }

  #map of the species imported for each analysis
  my %import_species = SPECIES;

  #if parent job was source specific import: check if it is part of
  # source specific run_type or a group run_type
  if (defined $source){
    # source only check run OR final source in set
    if ( ($run_type eq $source->{source_name}) ||
      ($groups_end_source{$run_type}
      && $groups_end_source{$run_type} eq $source->{source_name}) ) {
      $self->print_pipelogFH("Finished source only import ($source->{source_name} ".
                       $self->param('species').
                       ")\n");
    }
    elsif ($groups_end_source{$run_type}
         && $groups_end_source{$run_type} ne $source->{source_name}) {

      #source check and flow to next import, if not already at the end
      my $fwd = 1;

      if ($run_type eq ANIMALSET && $source->{source_name} eq OMIA) {
        my %animalQTL_species = map { $_ => 1 } @{$import_species{ANIMALQTL}};
        # only forward jobs to AnimalQTL for AnimalQTL species
        $fwd = 0 if (!defined($animalQTL_species{$species}));
      }

      if ($fwd) {
        #generally another import_source module
        $self->dataflow_output_id($self->param('output_ids'), 2);
      } else {
        #either check_phenotypes or disappear into the ether
        $self->dataflow_output_id($self->param('output_ids'), 3);
      }
    }

    #If not RGD, ZFIN imports, then this is the end of source related check - module
    if ($source->{source_name} ne RGD && $source->{source_name} ne ZFIN){
      close($self->logFH) if defined $self->logFH ;
      close($self->pipelogFH) if defined $self->pipelogFH ;
      return;
    }
  }

  # don't proceed with ontology mapping if only human variants were imported
  return if $run_type eq HUMAN_VAR;

  # check job was post import - final species check run
  # if species is an 'ontology term species' then move to dataflow 2
  my %ontology_species = map { $_ => 1 } @{$import_species{'ontology'}};
  if (defined ($ontology_species{$self->param('species')})) {
    if ($self->param('debug_mode')) {
      my $source = $source->{source_name} // '';
      $self->print_pipelogFH("Passing $source import (".
                       $self->param('species').
                       ") for adding ontology accessions (ontology_mapping)\n");
    }
    $self->dataflow_output_id($self->param('output_ids'), 2);
  } else {
    if ($self->param('debug_mode')) {
      my $source = $source->{source_name} // '';
      $self->print_pipelogFH("Passing $source import (".
                       $self->param('species').
                       ") for summary counts (finish_phenotype_annotation)\n");
    }
    $self->dataflow_output_id($self->param('output_ids'), 3);
  }
  close($self->logFH) if defined $self->logFH ;
  close($self->pipelogFH) if defined $self->pipelogFH ;
}


=head2 check_phenotype_description

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBConnection $dbh
               The new variation database connection
  Example    : check_phenotype_description($dbh)
  Description: Check for unsupported characters in phenotype names and print to report.
  Returntype : none
  Exceptions : none

=cut

sub check_phenotype_description{
  my ($self,$dbh) = @_;

  my $pheno_ext_sth = $dbh->prepare(qq[ select phenotype_id, description from phenotype ]);
  $pheno_ext_sth->execute()||die;

  my $ph =  $pheno_ext_sth->fetchall_arrayref();
  foreach my $l (@{$ph}){
    $self->print_logFH("WARNING: Phenotype id:$l->[0] has no description!\n") unless defined $l->[1];
    next unless defined $l->[1];
    $self->print_logFH("WARNING: Phenotype id:$l->[0] has empty description!\n") if $l->[1] eq '';

    my $full = $l->[1];
    #  my @matches = $l->[1] =~ /\(|\)|\/|\.|\; |\+|\'|\:|\@|\*|\%/gm;
    # ' can be ok: example: Kupffer's vesicle
    # / can be ok: example: G1/S transition of mitotic cell cycle
    # : can be ok: example: UDP-glucose:hexose-1-phosphate uridylyltransferase activity
    # . can be ok: example: Blond vs. brown hair color (from gwas)
    # + can be ok: example: decreased KLRG1+ CD8 alpha beta T cell number
    # () are currently tolerated, but not desired
    my @matches = $l->[1] =~ /\;|\?|\@|\*|\%/gm;
    $self->print_logFH("WARNING: Phenotype : $full (id:$l->[0]) looks suspect!\n") if(scalar(@matches) >0);

    # check for characters which will be interpreted a new lines
    @matches = $l->[1] =~ /.*\n.*/;
    $self->print_logFH("WARNING: Phenotype : $full (id:$l->[0]) contains a newline \n") if(scalar(@matches) >0);

    # check for phenotype descriptions suggesting no phenotype; these are now accepted
    $self->print_logFH("WARNING: Phenotype : $full (id:$l->[0]) is not useful \n") if !checkNonTerms( $l->[1] ) && $self->param('debug_mode');

    # check for unsupported individual character
    my $unsupportedChar = getUnsupportedChar($l->[1]);
    $self->print_logFH("WARNING: Phenotype : $full (id:$l->[0]) has suspect start or unsupported characters: $unsupportedChar \n") if defined($unsupportedChar);

  }

}


=head2 checkNonTerms

  Arg [1]    : string $description
               The phenotype description to be checked.
  Example    : checkNonTerms($description)
  Description: Check for known default (not meaningful) phenotype descriptions,
               returns 1 if no term was found, 0 if at least one term was matched.
  Returntype : boolean
  Exceptions : none

=cut

sub checkNonTerms {
  my $desc = shift;

  my $is_ok = 1;
  my @junk = ("None", "Not provided", "not specified", "Not in OMIM", "Variant of unknown significance", "not_provided", "?",".", "ClinVar: phenotype not specified");

  for my $check (@junk){
    if ( lc($desc) eq lc($check) ) {
      $is_ok  = 0;
      return $is_ok;
    }
  }

  return $is_ok;
}


=head2 getUnsupportedChar

  Arg [1]    : string $description
               The phenotype description to be checked.
  Example    : getUnsupportedChar($description)
  Description: Get unsupported characters in the phenotype description,
               returns nothing if nothing was found or first unsupported char that was matched.
  Returntype : undef or char
  Exceptions : none

=cut

sub getUnsupportedChar {
  my $desc = shift;

  my $is_ok = 1;
  my $i = 0;
  for my $c (split //, $desc) {
    # get ascii code
    my $ascii_val = ord($c);

    # check code in supported range
    if($ascii_val < 32 || $ascii_val  > 126 || $ascii_val == 60 || $ascii_val == 62 ){
      $is_ok = 0;
    }

    # also check first character makes sense
    if($i == 0 && ( $ascii_val < 48 ||
      ($ascii_val  > 57 && $ascii_val < 65) ||
      ($ascii_val  > 90 && $ascii_val < 97) ||
      $ascii_val  > 122)){
      $is_ok = 0;
    }
    if (!$is_ok) {
      return chr($ascii_val);
    }
    $i++;
  }
  return;
}


=head2 check_fk

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBConnection $dbh
               The new variation database connection
  Example    : check_fk($dbh)
  Description: Basic foreign keys checks on all phenotype data and print to report.
  Returntype : none
  Exceptions : none

=cut

sub check_fk{
  my ($self, $dbh) = @_;

  # check FKs: that all phenotypes have a phenotype_feature (pf),
  # all pf have a phenotype_feature_attrib (pfa), and vice versa
  my $featless_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype
                                                  where phenotype_id not in (select phenotype_id from phenotype_feature) ]);
  my $attribless_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype_feature
                                                 where phenotype_feature_id not in (select phenotype_feature_id from phenotype_feature_attrib) ]);
  my $featless_attrib_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype_feature_attrib
                                                  where phenotype_feature_id not in (select phenotype_feature_id from phenotype_feature) ]);
  my $phenoless_feat_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype_feature
                                                  where phenotype_id not in (select phenotype_id from phenotype) ]);
  my $phenoless_acc_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype_ontology_accession
                                                where phenotype_id not in (select phenotype_id from phenotype) ]);

  $featless_count_ext_sth->execute()||die;
  my $featless_count = $featless_count_ext_sth->fetchall_arrayref();
  $self->print_logFH("$featless_count->[0]->[0] phenotype entries with no phenotype_feature entry (expected: 0)\n");

  $attribless_count_ext_sth->execute()||die;
  my $attrless_count = $attribless_count_ext_sth->fetchall_arrayref();
  $self->print_logFH("$attrless_count->[0]->[0] phenotype_feature entries with no phenotype_feature_attrib entry (can be valid cases)\n");

  $featless_attrib_count_ext_sth->execute()||die;
  my $featless_attrib_count = $featless_attrib_count_ext_sth->fetchall_arrayref();
  $self->print_logFH("$featless_attrib_count->[0]->[0] phenotype_feature_attrib entries with missing phenotype_feature entry (expected: 0)\n");

  $phenoless_feat_count_ext_sth->execute()||die;
  my $phenoless_count = $phenoless_feat_count_ext_sth->fetchall_arrayref();
  $self->print_logFH("$phenoless_count->[0]->[0] phenotype_feature entries with missing phenotype entry (expected: 0)\n");

  $phenoless_acc_count_ext_sth->execute()||die;
  my $phenoless_acc_count = $phenoless_acc_count_ext_sth->fetchall_arrayref();
  $self->print_logFH("$phenoless_acc_count->[0]->[0] phenotype_ontology_accession rows with missing phenotype entry (expected: 0)\n");

}

sub check_source {
  my $self = shift;

  my $source = $self->param('source');
  my $assembly = $self->get_assembly;

  ## retrieve old results from production db for comparison if available
  my $previous_counts = $self->get_old_results();

  ## calculate new counts from new phenotype_feature table
  my $new_counts      = $self->get_new_results();

  my $text_out = "\nSummary of results from CheckPhenotypeAnnotation (if less counts than previous)\n\n";

  # check only source specific if source specific check
  if (defined $source && defined $previous_counts->{phenotype_feature_count_details} ){

    my $source_name = $source->{source_name};
    $source_name = 'Animal_QTLdb_QTL' if $source_name eq 'AnimalQTLdb';
    my @check_names= grep {/$source_name/ } keys $previous_counts->{phenotype_feature_count_details};
    if (scalar(@check_names) != 1) {
      $text_out .= "WARNING:".scalar(@check_names) . " check results found for '" . $source->{source_name} . "', where 1 expected (0 is expected for EGA): " . @check_names . " \n";
      $self->print_logFH($text_out);
      return;
    }

    my $check_name =$check_names[0];
    if (defined $previous_counts->{phenotype_feature_count_details}{$check_name} &&
        $previous_counts->{phenotype_feature_count_details}{$check_name}  > $new_counts->{phenotype_feature_count_details}{$check_name}){
      $text_out .= "WARNING: " . $new_counts->{phenotype_feature_count_details}{$check_name} . " $check_name entries";
      $text_out .= " (previously " . $previous_counts->{phenotype_feature_count_details}{$check_name} . ")" ;
      $text_out .= "\n";
      #for grch37 do not fail the job as prev counts are retrieved by species and not by assembly and grch38 are always reported
      $count_ok = 0 unless $assembly eq 'GRCh37';
    }
  } else {
    my @tables = ('phenotype', 'phenotype_feature', 'phenotype_feature_attrib', 'phenotype_ontology_accession');
    foreach my $table (@tables) {
      my $check_name = "$table\_count";

      if (defined $previous_counts->{$check_name} &&
          $previous_counts->{$check_name}  > $new_counts->{$check_name}){
        $text_out .= "WARNING: " . $new_counts->{"$table\_count"} . " $table entries";
        $text_out .= " (previously " . $previous_counts->{"$table\_count"} . ")" if defined  $previous_counts->{"$table\_count"} ;
        $text_out .= "\n";
        #for grch37 do not fail the job as prev counts are retrieved by species and not by assembly and grch38 are always reported
        $count_ok = 0 unless $assembly eq 'GRCh37';
      }
    }
  }
  $self->print_logFH($text_out);
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

  # only source specific checks will update the meta table
  my $source_info = $self->param("source");
  return if !defined $source_info;

  my $source_name;
  $source_name = $source_info->{source_name} if defined $source_info;
  $source_name ||= $self->param("source_name");
  my $var_dbh = $self->variation_db_adaptor->dbc->db_handle;

  my $update_meta_sth = $var_dbh->prepare(qq[ insert ignore into meta
                                              ( meta_key, meta_value) values (?,?)
                                            ]);

  $update_meta_sth->execute('PhenotypeAnnotation_run_date_'.$source_name, $self->run_date() );

}

1;


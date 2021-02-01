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


=head1 ImportEGA

This module imports EGA (European Genome-phenome Archive) study data.
Note: only studies are imported, no phenotype feature data.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportEGA;

use warnings;
use strict;

use File::Path qw(make_path);
use POSIX qw(strftime);
use DBI qw(:sql_types);

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;


=head2 ega_db_handler

  Arg [1]    : DBI::db $db_handler (optional)
               The new ega_db_handler
  Example    : $ega_dbh = $obj->ega_db_handler()
  Description: Get/set the ega_db_handler
  Returntype : DBI::db
  Exceptions : none

=cut

sub ega_db_handler {
  my ($self, $db_handler) = @_;
  $self->{ega_dbh} = $db_handler if defined $db_handler;
  return $self->{ega_dbh};
}

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $conf_file    = $self->required_param('ega_database_conf');

  $self->debug($self->param('debug_mode'));

  %source_info = (source_description => 'Variants imported from the European Genome-phenome Archive with phenotype association',
                  source_url => 'https://www.ebi.ac.uk/ega/',
                  object_type => 'Variation',
                  source_status => 'germline',
                  source_name => 'EGA',        #source name in the variation db
                  source_name_short => 'EGA',  #source identifier in the pipeline
                  data_types => 'study',
                  );
  $source_info{source_version} = strftime("%Y%m", localtime); # it is current month

  my $workdir = $pipeline_dir."/".$source_info{source_name}."/".$species;
  make_path($workdir) or die "Failed to create $workdir $!\n";
  $self->workdir($workdir);
  my $file_ega = "ega.studies.csv";

  open(my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  #parse database connection details:
  my %database_conf;
  open(CONF,'<',$conf_file) || die ("Could not open $conf_file for reading\n");
  while (<CONF>) {
      chomp;                  # no newline
      s/#.*//;                # no comments
      s/^\s+//;               # no leading white
      s/\s+$//;               # no trailing white
      next unless length;     # anything left?
      my ($var, $value) = split(/\s*=\s*/, $_, 2);
      $database_conf{$var} = $value;
  }

  my $dsn = "dbi:mysql:$database_conf{DATABASE}:$database_conf{HOST}:$database_conf{PORT}";
  my $dbh = DBI->connect($dsn,$database_conf{USER},$database_conf{PASS});
  $self->ega_db_handler($dbh);

  print $logFH "Found files (".$workdir."/".$file_ega."), will skip new fetch\n" if -e $workdir."/".$file_ega;
  $self->get_input_file($file_ega) unless -e $workdir."/".$file_ega;

  $self->param('ega_file', $file_ega);
}

sub run {
  my $self = shift;

  my $file_ega = $self->required_param('ega_file');

  #get source id
  my $source_id = $self->get_or_add_source(\%source_info);
  $self->print_logFH("$source_info{source_name} source_id is $source_id\n") if ($self->debug);

  # get study data + save into study and associate_study (all in one method)
  my $results = $self->parse_and_add_input_file($file_ega, $source_id);
  $self->print_logFH("Got ".(scalar @{$results->{'studies'}})." new studies \n") if ($self->debug);

  my %param_source = (source_name => $source_info{source_name},
                      type => $source_info{object_type});
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species'),
                               run_type => $self->required_param('run_type'),
                             });
}

sub write_output {
  my $self = shift;

  $self->print_pipelogFH("Passing $source_info{source_name_short} import (".$self->required_param('species').") for checks (check_phenotypes)\n") if ($self->debug);
  close($self->logFH) if defined $self->logFH ;
  close($self->errFH) if defined $self->errFH ;
  close($self->pipelogFH) if defined $self->pipelogFH ;

  $self->dataflow_output_id($self->param('output_ids'), 2);
}


=head2 get_input_file

  Arg [1]    : String $outfile
              The output file path
  Example    : $self->get_input_file($outfile)
  Description: Fetches EGA data and writes the data comma separated into $outfile.
  Returntype : none.
  Exceptions : none

=cut

sub get_input_file {
  my ($self, $outfile) = @_ ;

  my $errFH1;
  open($errFH1, ">", $self->workdir."/".'log_import_err_'.$outfile) || die ("Could not open file for writing: $!\n");

  my $studies = $self->get_all_study_stable_ids();
  open(FILE, ">".$self->workdir."/".$outfile) || die ("Could not open file for writing: $!\n");

  foreach my $stable_id (sort {$a cmp $b } @$studies) {
    my $study_id = $self->get_study_id($stable_id);
    my $pmid = $self->get_publication($study_id);
    if (!$pmid) {
      print $errFH1 "WARNING: Cannot find a publication for $stable_id\n";
      next;
    }

    print FILE "$stable_id,$pmid,$source_info{source_url}studies/$stable_id\n";
  }

  close (FILE);
  close ($errFH1);
}


=head2 get_all_study_stable_ids

  Example    : $self->get_all_study_stable_ids()
  Description: Fetch stable_id(s) from study table in EGA
  Returntype : arrayref of stable_id(s)
  Exceptions : none

=cut

sub get_all_study_stable_ids {
  my $self = shift ;

  my @studies;

  my $sth = $self->ega_db_handler->prepare(qq[SELECT stable_id FROM study]);
  $sth->execute();

  while (my $id = $sth->fetchrow() ) {
    push @studies, $id;
  }

  $sth->finish();

  return \@studies;
}


=head2 get_study_id

  Arg [1]    : String $stable_id
  Example    : $self->get_study_id($stable_id)
  Description: Fetch study_id from study table in EGA for $stable_id
  Returntype : string study_id
  Exceptions : none

=cut

sub get_study_id {
  my ($self, $stable_id) = @_;

  my $sth = $self->ega_db_handler->prepare(qq[SELECT study_id FROM study WHERE stable_id = ?]);
  $sth->execute($stable_id);

  my $study_id = $sth->fetchrow();
  $sth->finish();

  return $study_id;
}


=head2 get_publication

  Arg [1]    : String $study_id
  Example    : $self->get_publication($study_id)
  Description: Fetch publication_id from study table in EGA for $study_id
  Returntype : string publication_id
  Exceptions : none

=cut

sub get_publication {
    my ($self, $study_id) = @_;

    my $sth = $self->ega_db_handler->prepare(qq[SELECT publication_id FROM study_publication WHERE study_id = ?]);
    $sth->execute($study_id);

    my $publication_id = $sth->fetchrow();
    $sth->finish();

    my $pmid;
    if ($publication_id) {
      $sth = $self->ega_db_handler->prepare(qq[SELECT pmid from publication WHERE publication_id = ?]);
      $sth->execute($publication_id);

      $pmid = $sth->fetchrow();
      $sth->finish();
    }

    return $pmid;
}


=head2 parse_and_add_input_file

  Arg [1]    : String $infile
  Arg [2]    : String $source_id
  Example    : $self->parse_and_add_input_file($infile, $source_id)
  Description: Parse EGA studies file and insert new data into study and associate_study tables in variatio db.
  Returntype : string study_id
  Exceptions : none

=cut

sub parse_and_add_input_file {
  my ($self, $infile, $source_id) = @_;

  my $variation_dba = $self->variation_db_adaptor;
  my $errFH2;
  open($errFH2, ">", $self->workdir."/".'log_import_err_parse_'.$infile) ;

  my $study_check_stmt = qq{
    SELECT
      study_id
    FROM
      study
    WHERE
      name=? AND source_id=$source_id
    LIMIT 1
  };
  my $nhgri_check_stmt = qq{
    SELECT
      st.study_id,st.study_type
    FROM
      study st, source s
    WHERE
      external_reference=? AND s.name like '%nhgri%'
      AND s.source_id=st.source_id
    LIMIT 1
  };
  my $study_ins_stmt = qq{
    INSERT INTO
      study (
      name,
      source_id,
      external_reference,
      url,
      study_type
      )
    VALUES (
      ?,
      $source_id,
      ?,
      ?,
      ?
    )
  };
  # NHGRI and EGA associated studies
  my $asso_study_check_stmt = qq{
    SELECT
      study1_id
    FROM
      associate_study
    WHERE
      study1_id = ? AND study2_id = ?
    LIMIT 1
  };
  my $asso_study_ins_stmt = qq{
    INSERT INTO
      associate_study (study1_id,study2_id)
    VALUES (?,?)
  };

  my $nhgri_check_sth = $variation_dba->dbc->prepare($nhgri_check_stmt);
  my $study_check_sth = $variation_dba->dbc->prepare($study_check_stmt);
  my $study_ins_sth   = $variation_dba->dbc->prepare($study_ins_stmt);
  my $asso_study_check_sth = $variation_dba->dbc->prepare($asso_study_check_stmt);
  my $asso_study_ins_sth   = $variation_dba->dbc->prepare($asso_study_ins_stmt);

  # Open the input file for reading
  open(IN,'<',$self->workdir."/".$infile) || die ("Could not open $infile for reading\n");

  # Read through the file and parse out the desired fields
  my @new_studies;
  while (<IN>) {
    chomp $_;
    my @attributes = split(",",$_);
    next if ($attributes[1] eq '');
    my $name = $attributes[0];
    my $pubmed = $self->get_pubmed_prefix.$attributes[1];
    my $url = $attributes[2];

    # NHGRI study
    my $nhgri_study_id;
    my $study_type;
    $nhgri_check_sth->bind_param(1,$pubmed,SQL_VARCHAR);
    $nhgri_check_sth->execute();
    $nhgri_check_sth->bind_columns(\$nhgri_study_id,\$study_type);
    $nhgri_check_sth->fetch();

    if (!defined($nhgri_study_id)) {
      print $errFH2 "No NHGRI-EBI study found for the EGA $name | $pubmed !\n";
      next;
    }

    # EGA study
    my $study_id;
    $study_check_sth->bind_param(1,$name,SQL_VARCHAR);
    $study_check_sth->execute();
    $study_check_sth->bind_columns(\$study_id);
    $study_check_sth->fetch();
    if (!defined($study_id)) {
      $study_ins_sth->bind_param(1,$name,SQL_VARCHAR);
      $study_ins_sth->bind_param(2,$pubmed,SQL_VARCHAR);
      $study_ins_sth->bind_param(3,$url,SQL_VARCHAR);
      $study_ins_sth->bind_param(4,$study_type,SQL_VARCHAR);
      $study_ins_sth->execute();
      $study_id = $variation_dba->dbc->db_handle->{'mysql_insertid'};
      push (@new_studies, $name)
    }

    my $is_associated;
    $asso_study_check_sth->bind_param(1,$nhgri_study_id,SQL_INTEGER);
    $asso_study_check_sth->bind_param(2,$study_id,SQL_INTEGER);
    $asso_study_check_sth->execute();
    $asso_study_check_sth->bind_columns(\$is_associated);
    $asso_study_check_sth->fetch();

    if (!defined($is_associated)) {
      $asso_study_ins_sth->bind_param(1,$nhgri_study_id,SQL_INTEGER);
      $asso_study_ins_sth->bind_param(2,$study_id,SQL_INTEGER);
      $asso_study_ins_sth->execute();
    }
  }
  close(IN);
  close($errFH2);

  my %result = ('studies' => \@new_studies);
  return \%result;
}

1;

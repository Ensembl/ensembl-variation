package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportGENCC;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation/;

use Text::CSV_XS;
use POSIX qw(strftime);
use File::Path qw(make_path);
use HTTP::Tiny;
use Bio::EnsEMBL::Registry;

sub fetch_input {
  my ($self) = @_;

  my $species = $self->param_required('species');
  my $pdir    = $self->param_required('pipeline_dir');

  my $workdir = $pdir . "/GenCC/" . $species;
  unless (-d $workdir) {
    my $err;
    make_path($workdir, { error => \$err });
    $self->throw("make_path failed to create $workdir") if $err && @$err;
  }
  $self->workdir($workdir);

  open(my $logFH,     ">", $workdir."/log_import_out_GenCC_".$species)
    || $self->throw("Failed to open log_import_out_GenCC_$species: $!");
  open(my $errFH,     ">", $workdir."/log_import_err_GenCC_".$species)
    || $self->throw("Failed to open log_import_err_GenCC_$species: $!");
  open(my $pipelogFH, ">", $workdir."/log_import_debug_pipe_GenCC_".$species)
    || $self->throw("Failed to open log_import_debug_pipe_GenCC_$species: $!");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  my $default_dir  = $pdir . "/GenCC";
  my $default_file = $default_dir . "/gencc.csv";
  my $file         = $self->param('gencc_file') // $default_file;

  # If file is missing/empty, fetch the latest CSV from GenCC
  unless (-s $file) {
    make_path($default_dir);
    my $url  = 'https://search.thegencc.org/download/action/submissions-export-csv';
    my $http = HTTP::Tiny->new(timeout => 180);
    my $res  = $http->get($url);
    $self->throw("GenCC download failed: $res->{status} $res->{reason}") unless $res->{success};
    open(my $OUT, ">:encoding(UTF-8)", $default_file) or $self->throw("Cannot write $default_file: $!");
    print $OUT $res->{content};
    close $OUT;
    $file = $default_file;
  }

  -s $file or $self->throw("GenCC input CSV not found or empty after download: $file");

  my $version = $self->param('gencc_version') // strftime("%Y%m%d", localtime);

  $self->param('input_file',    $file);
  $self->param('gencc_version', $version);
}

sub run {
  my ($self) = @_;

  my $species = $self->param_required('species');
  my $file    = $self->param_required('input_file');
  my $version = $self->param_required('gencc_version');

  my %skip = map { $_ => 1 } @{ $self->param('filter_submitters') || [ 'G2P' ] };

  # DB adaptors
  my $core_dba = $self->get_species_adaptor('core');
  my $var_dba  = $self->get_species_adaptor('variation');

  my $dbh_core = $core_dba->dbc->db_handle;
  my $dbh_var  = $var_dba->dbc->db_handle;

  ############################
  # Source info
  ############################
  my %source_info = (
    source_description => 'Gene Curation Coalition (gene–disease validity assertions)',
    source_url         => 'https://search.thegencc.org/',
    object_type        => 'Gene',
    source_status      => 'germline',
    source_name        => 'GenCC',
    source_name_short  => 'GenCC',
    source_version     => $version,
    data_types         => 'phenotype_feature',
  );
  $self->param('source_info', \%source_info);

  # Ensure source exists or get id
  my $source_id = $self->get_or_add_source(\%source_info);
  $self->print_logFH("$source_info{source_name} source_id is $source_id\n")
    if ($self->param('debug_mode'));

  ############################
  # Ensure attrib_types exist
  ############################
  my @required_codes = qw(
    gencc_submitter
    gencc_classification
    gencc_inherit_mode
    gencc_uuid
  );

  my %atid;
  for my $code (@required_codes) {
    my ($id) = $dbh_var->selectrow_array(
      'SELECT attrib_type_id FROM attrib_type WHERE code=?', undef, $code
    );
    $self->throw("Missing attrib_type '$code' in variation DB.") unless $id;
    $atid{$code} = $id;
  }

  ############################
  # SQL commands
  ############################
  my $dbq_hgnc_to_ensg = $dbh_core->prepare(q{
    SELECT g.stable_id
      FROM gene g
      JOIN object_xref ox ON ox.ensembl_id=g.gene_id AND ox.ensembl_object_type='Gene'
      JOIN xref x         ON x.xref_id=ox.xref_id
      JOIN external_db ed ON ed.external_db_id=x.external_db_id
     WHERE ed.db_name='HGNC' AND x.dbprimary_acc=? LIMIT 1
  });

  my $dbq_phen_sel = $dbh_var->prepare("SELECT phenotype_id FROM phenotype WHERE description=? LIMIT 1");
  my $dbq_phen_ins = $dbh_var->prepare("INSERT INTO phenotype (description) VALUES (?)");

  my $dbq_pf_exists = $dbh_var->prepare(q{
    SELECT 1 FROM phenotype_feature
     WHERE phenotype_id=? AND source_id=? AND type='Gene' AND object_id=? LIMIT 1
  });

  my $dbq_pf_ins = $dbh_var->prepare(q{
    INSERT INTO phenotype_feature (phenotype_id, source_id, type, object_id)
    VALUES (?,?, 'Gene', ?)
  });

  my $dbq_pfa_ins = $dbh_var->prepare(q{
    INSERT INTO phenotype_feature_attrib (phenotype_feature_id, attrib_type_id, value)
    VALUES (?,?,?)
  });

  ############################
  # Parse CSV and import
  ############################
  my $csv = Text::CSV_XS->new({ binary => 1, auto_diag => 1 });
  open(my $IN, "<:encoding(UTF-8)", $file) or $self->throw("Cannot open $file: $!");

  my $header = $csv->getline($IN) or $self->throw("Empty CSV: $file");
  $csv->column_names(@$header);

  my ($n_rows, $n_skip, $n_no_gene, $n_dupe, $n_loaded) = (0,0,0,0,0);

  ROW: while (my $r = $csv->getline_hr($IN)) {
    $n_rows++;

    # Filter submitter
    my $submitter = $r->{'submitter_title'} // '';
    $submitter =~ s/^\s+|\s+$//g;
    if ($submitter && $skip{$submitter}) { $n_skip++; next ROW; }

    # Minimal required field
    my $disease = $r->{'disease_title'} // '';
    next ROW unless $disease;

    # Optional attributes
    my $uuid  = $r->{'uuid'} // '';
    my $class = $r->{'classification_title'} // '';
    my $moi   = $r->{'moi_title'} // '';
    for ($class, $moi) { $_ =~ s/^\s+|\s+$//g if defined $_ }

    # Resolve gene: HGNC first, then symbol
    my $ensg;
    if ((my $curie = $r->{'gene_curie'} // '') =~ /^HGNC:(\d+)$/) {
      $dbq_hgnc_to_ensg->execute($1);
      ($ensg) = $dbq_hgnc_to_ensg->fetchrow_array;
    }
    unless ($ensg) {
      my $sym = $r->{'gene_symbol'} // '';
      if ($sym) {
        my $ga   = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'Gene');
        my $hits = $ga ? ($ga->fetch_all_by_external_name($sym) || []) : [];
        $ensg = @$hits ? $hits->[0]->stable_id : undef;
        unless ($ensg) {
          my $g = $ga ? eval { $ga->fetch_by_display_label($sym) } : undef;
          $ensg ||= $g && $g->stable_id;
        }
      }
    }
    unless ($ensg) { $n_no_gene++; next ROW; }

    # Ensure phenotype row (by description)
    $dbq_phen_sel->execute($disease);
    my ($phenotype_id) = $dbq_phen_sel->fetchrow_array;
    unless ($phenotype_id) {
      $dbq_phen_ins->execute($disease);
      $phenotype_id = $dbh_var->{mysql_insertid};
    }

    # Avoid duplicate Phenotype Feature rows
    $dbq_pf_exists->execute($phenotype_id, $source_id, $ensg);
    my ($exists) = $dbq_pf_exists->fetchrow_array;
    if ($exists) { $n_dupe++; next ROW; }

    # Insert link
    $dbq_pf_ins->execute($phenotype_id, $source_id, $ensg);
    my $pf_id = $dbh_var->{mysql_insertid};

    # Attributes
    $dbq_pfa_ins->execute($pf_id, $atid{gencc_submitter},      $submitter) if $submitter;
    $dbq_pfa_ins->execute($pf_id, $atid{gencc_classification}, $class)     if $class;
    $dbq_pfa_ins->execute($pf_id, $atid{gencc_inherit_mode},   $moi)       if $moi;
    $dbq_pfa_ins->execute($pf_id, $atid{gencc_uuid},           $uuid)      if $uuid;

    $n_loaded++;
  }
  close $IN;

  $self->warning(
    "ImportGenCC: rows=$n_rows, skipped_submitter=$n_skip, unresolved_gene=$n_no_gene, duplicate_links=$n_dupe, loaded=$n_loaded"
  );
  $self->param('loaded_count', $n_loaded);

  my %param_source = (
    source_name => $source_info{source_name_short},
    type        => $source_info{object_type}
  );
  $self->param('output_ids', {
    source   => \%param_source,
    species  => $species,
    run_type => 'GENCC',
  });
}

sub write_output {
  my $self = shift;

  my $debug = $self->param('debug_mode');
  my $si    = $self->param('source_info') || {};

  if ($debug && $self->pipelogFH && $si->{source_name_short}) {
    $self->print_pipelogFH(
      "Passing $si->{source_name_short} import (".$self->required_param('species').") for checks (check_gencc)\n"
    );
  }

  close($self->logFH)     if defined $self->logFH;
  close($self->errFH)     if defined $self->errFH;
  close($self->pipelogFH) if defined $self->pipelogFH;

  my $output_ids = $self->param('output_ids') || {};
  $self->dataflow_output_id($output_ids, 2);
}

1;
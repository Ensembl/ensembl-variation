package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportGENCC;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation/;

use Text::CSV_XS;
use POSIX qw(strftime);
use File::Path qw(make_path);
use HTTP::Tiny;
use Bio::EnsEMBL::Registry;
use Encode qw(encode decode);
use Unicode::Normalize qw(NFD);

# revert UTF-8 to Latin1 incorrect decoding
sub _demojibake {
  my ($s) = @_;
  return $s unless defined $s;
  if ($s =~ /Ã|Â|Ãƒ|â|€/) {
    my $bytes = encode('latin1', $s);
    my $fixed = eval { decode('utf8', $bytes, 1) };
    $s = $fixed if defined $fixed;
  }
  return $s;
}

# Clean phenotype text
sub _clean_text {
  my ($s) = @_;
  return undef unless defined $s;

  # remove BOM & control chars as these trigger "unsupported characters"
  $s =~ s/^\x{FEFF}//;
  $s =~ s/[\x00-\x08\x0B\x0C\x0E-\x1F\x7F]//g;
  $s =~ s/[\x{0080}-\x{009F}]//g;
  $s =~ s/\x{FFFD}//g;

  # undo incorrect decoding (e.g. StÃ¼ve to Stüve)
  $s = _demojibake($s);

  # normalise spaces & punctuation; keep semantics
  $s =~ s/[\x{00A0}\x{202F}\x{2000}-\x{200A}\x{205F}\x{3000}]/ /g;
  $s =~ s/[\x{2010}-\x{2015}\x{2212}]/-/g;    # dashes to hyphen
  $s =~ s/[\x{2018}\x{2019}\x{201B}]/'/g;     # smart single quotes to '
  $s =~ s/[\x{201C}\x{201D}]/"/g;             # smart double quotes to "

  # remove accents but keep the base letters
  $s = NFD($s);
  $s =~ s/\pM//g;

  # map Latin letters that didn't get decoded to ASCII
  $s =~ tr/ßØøÆæŒœŁłÐðÞþ/ssOoAEaeOEoeLlDdThth/;

  # collapse interior spaces and trim ends (otherwise dups appear that aren't real dups)
  $s =~ s/\s{2,}/ /g;
  $s =~ s/^\s+//;
  $s =~ s/\s+$//;

  # capitalise first letter if it was lower-case
  $s =~ s/^([a-z])/\U$1/;

  # final check to ensure ASCII only (post accent stripping)
  $s =~ s/[^\x00-\x7F]//g;

  return $s;
}

sub fetch_input {
  my ($self) = @_;

  my $species = $self->param_required('species');
  my $pdir    = $self->param_required('pipeline_dir');

  # Make a working directory for this species
  my $workdir = $pdir . "/GenCC/" . $species;
  unless (-d $workdir) {
    my $err;
    make_path($workdir, { error => \$err });
    $self->throw("make_path failed to create $workdir") if $err && @$err;
  }
  $self->workdir($workdir);

  # Log files as UTF-8
  open(my $logFH, ">", $workdir."/log_import_out_GenCC_".$species)
    || $self->throw("Failed to open log_import_out_GenCC_$species: $!");
  open(my $errFH, ">", $workdir."/log_import_err_GenCC_".$species)
    || $self->throw("Failed to open log_import_err_GenCC_$species: $!");
  open(my $pipelogFH, ">", $workdir."/log_import_debug_pipe_GenCC_".$species)
    || $self->throw("Failed to open log_import_debug_pipe_GenCC_$species: $!");
  binmode($logFH,     ':encoding(UTF-8)');
  binmode($errFH,     ':encoding(UTF-8)');
  binmode($pipelogFH, ':encoding(UTF-8)');

  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  # Locate or download the GenCC CSV
  my $default_dir  = $pdir . "/GenCC";
  my $default_file = $default_dir . "/gencc.csv";
  my $file         = $self->param('gencc_file') // $default_file;

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

  # Allow excluding submitters, i.e. G2P is the default
  my %skip = map { $_ => 1 } @{ $self->param('filter_submitters') || [ 'G2P' ] };

  # Get DB handles
  my $core_dba = $self->get_species_adaptor('core');
  my $var_dba  = $self->get_species_adaptor('variation');
  my $dbh_core = $core_dba->dbc->db_handle;
  my $dbh_var  = $var_dba->dbc->db_handle;

  # Fetch Gene adaptor
  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'Gene');

  # Source metadata
  my %source_info = (
    source_description => 'Gene Curation Coalition (gene disease validity assertions)',
    source_url         => 'https://search.thegencc.org/',
    object_type        => 'Gene',
    source_status      => 'germline',
    source_name        => 'GenCC',
    source_name_short  => 'GenCC',
    source_version     => $version,
    data_types         => 'phenotype_feature',
  );
  $self->param('source_info', \%source_info);

  # Ensure source row exists
  my $source_id = $self->get_or_add_source(\%source_info);
  $self->print_logFH("$source_info{source_name} source_id is $source_id\n")
    if ($self->param('debug_mode'));

  # Ensure required attrib types exist
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

  # SQL
  my $dbq_hgnc_to_ensg = $dbh_core->prepare(q{
    SELECT g.stable_id
      FROM gene g
      JOIN object_xref ox ON ox.ensembl_id=g.gene_id AND ox.ensembl_object_type='Gene'
      JOIN xref x         ON x.xref_id=ox.xref_id
      JOIN external_db ed ON ed.external_db_id=x.external_db_id
     WHERE ed.db_name='HGNC' AND x.dbprimary_acc=? LIMIT 1
  });
  my $dbq_phen_sel = $dbh_var->prepare('SELECT phenotype_id FROM phenotype WHERE description=? LIMIT 1');
  my $dbq_phen_ins = $dbh_var->prepare('INSERT INTO phenotype (description) VALUES (?)');
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

  # Read CSV, clean fields, and load rows
  my $csv = Text::CSV_XS->new({ binary => 1, auto_diag => 1 });
  open(my $IN, "<:encoding(UTF-8)", $file) or $self->throw("Cannot open $file: $!");

  my $header = $csv->getline($IN) or $self->throw("Empty CSV: $file");
  $csv->column_names(@$header);

  my ($n_rows, $n_skip, $n_no_gene, $n_dupe, $n_loaded) = (0,0,0,0,0);

  ROW: while (my $r = $csv->getline_hr($IN)) {
    $n_rows++;

    my $submitter = _clean_text($r->{'submitter_title'}      // '');
    my $disease   = _clean_text($r->{'disease_title'}        // '');
    my $uuid      = _clean_text($r->{'uuid'}                 // '');
    my $class     = _clean_text($r->{'classification_title'} // '');
    my $moi       = _clean_text($r->{'moi_title'}            // '');

    next ROW if $submitter && $skip{$submitter};
    next ROW unless $disease;

    # Find the gene (HGNC id first, fallback to symbol)
    my $ensg;
    if ((my $curie = $r->{'gene_curie'} // '') =~ /^HGNC:(\d+)$/) {
      $dbq_hgnc_to_ensg->execute($1);
      ($ensg) = $dbq_hgnc_to_ensg->fetchrow_array;
    }
    unless ($ensg) {
      my $sym = $r->{'gene_symbol'} // '';
      if ($sym) {
        my $hits = $gene_adaptor ? ($gene_adaptor->fetch_all_by_external_name($sym) || []) : [];
        $ensg = @$hits ? $hits->[0]->stable_id : undef;
        unless ($ensg) {
          my $g = $gene_adaptor ? eval { $gene_adaptor->fetch_by_display_label($sym) } : undef;
          $ensg ||= $g && $g->stable_id;
        }
      }
    }
    unless ($ensg) { $n_no_gene++; next ROW; }

    # Create/find phenotype by cleaned description
    $dbq_phen_sel->execute($disease);
    my ($phenotype_id) = $dbq_phen_sel->fetchrow_array;
    unless ($phenotype_id) {
      $dbq_phen_ins->execute($disease);
      $phenotype_id = $dbh_var->{mysql_insertid};
    }

    # Avoid duplicate links
    $dbq_pf_exists->execute($phenotype_id, $source_id, $ensg);
    my ($exists) = $dbq_pf_exists->fetchrow_array;
    if ($exists) { $n_dupe++; next ROW; }

    # Link and add attributes
    $dbq_pf_ins->execute($phenotype_id, $source_id, $ensg);
    my $pf_id = $dbh_var->{mysql_insertid};
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

  # Go to downstream checks
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
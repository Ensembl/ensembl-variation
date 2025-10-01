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
sub _utf8tolatin {
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
  $s = _utf8tolatin($s);

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

  # locate or download the GenCC CSV
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

  # Allow excluding submitters
  my %skip = map { $_ => 1 } @{ $self->param('filter_submitters') || [ 'G2P', 'Orphanet' ] };

  # Get DB handles
  my $core_dba = $self->get_species_adaptor('core');
  my $var_dba  = $self->get_species_adaptor('variation');

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

  # Dump and clean pre-existing phenotype features
  $self->dump_phenotypes($source_info{source_name}, 1);

  # Get phenotype data + save it (all in one method)
  my $results = $self->parse_input_file($file, \%skip, $gene_adaptor);

  # Save phenotypes
  $self->save_phenotypes(\%source_info, $results);

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

sub parse_input_file {
  my ($self, $infile, $skip_href, $gene_adaptor) = @_;

  my $dbh_core = $self->get_species_adaptor('core')->dbc->db_handle;

  # HGNC id to ENSG (accept 32700 or HGNC:32700), exclude LRG_*
  my $sth_hgnc_to_ensg = $dbh_core->prepare(q{
    SELECT g.stable_id
      FROM gene g
      JOIN object_xref ox ON ox.ensembl_id = g.gene_id AND ox.ensembl_object_type = 'Gene'
      JOIN xref x         ON x.xref_id = ox.xref_id
      JOIN external_db ed ON ed.external_db_id = x.external_db_id
     WHERE ed.db_name = 'HGNC'
       AND (x.dbprimary_acc = ? OR x.dbprimary_acc = CONCAT('HGNC:', ?))
       AND g.stable_id NOT LIKE 'LRG\_%'
     LIMIT 1
  });

  my $csv = Text::CSV_XS->new({ binary => 1, auto_diag => 1 });
  open(my $IN, "<:encoding(UTF-8)", $infile) or $self->throw("Cannot open $infile: $!");

  my $header = $csv->getline($IN) or $self->throw("Empty CSV: $infile");
  $csv->column_names(@$header);

  my @phenotypes;
  my ($n_rows, $n_skip, $n_unmatched) = (0,0,0);

  ROW: while (my $r = $csv->getline_hr($IN)) {
    $n_rows++;

    my $submitter = _clean_text($r->{'submitter_title'}      // '');
    my $disease   = _clean_text($r->{'disease_title'}        // '');
    my $uuid      = _clean_text($r->{'uuid'}                 // '');
    my $class     = _clean_text($r->{'classification_title'} // '');
    my $moi       = _clean_text($r->{'moi_title'}            // '');
    my $curie     = $r->{'gene_curie'}  // '';
    my $symbol    = $r->{'gene_symbol'} // '';

    # Skip excluded submitters
    if ($submitter && $skip_href->{$submitter}) { $n_skip++; next ROW; }

    # Skip empty disease with log
    unless ($disease) {
      $self->print_logFH(join("\t","UNMATCHED","no_identifier",($uuid//''),($submitter//''),($curie//''),($symbol//''),'')."\n");
      $n_unmatched++; next ROW;
    }

    # Resolving the gene mapping - HGNC id is preferred, or HGNC symbol (exclude LRG)
    my $ensg;
    if ($curie =~ /^HGNC:(\d+)$/) {
      $sth_hgnc_to_ensg->execute($1, $1);
      ($ensg) = $sth_hgnc_to_ensg->fetchrow_array;
    }
    unless ($ensg && $ensg !~ /^LRG_/) {
      if ($symbol && $gene_adaptor) {
        my $genes = $gene_adaptor->fetch_all_by_external_name($symbol, 'HGNC') || [];
        @$genes = grep { $_->stable_id !~ /^LRG_/ } @$genes;  # drop LRGs
        if (@$genes > 1) {
          my @tmp = grep { (($_->external_name // '') eq $symbol) } @$genes; # exact symbol tie-break
          $genes = \@tmp if @tmp;
        }
        $ensg = @$genes ? $genes->[0]->stable_id : undef;
      }
    }

    unless ($ensg) {
      my $reason = $curie =~ /^HGNC:/ ? 'HGNC_not_found'
                 : $symbol            ? 'symbol_not_found'
                 :                      'no_identifier';
      $self->print_logFH(join("\t","UNMATCHED",$reason,($uuid//''),($submitter//''),$curie,$symbol,$disease)."\n");
      $n_unmatched++; next ROW;
    }

    my @acc = ();
    if (defined $r->{'disease_curie'} && $r->{'disease_curie'} ne '') {
      my $acc = $r->{'disease_curie'};
      $acc =~ s/_/:/; # e.g. MONDO_0012345 to MONDO:0012345
      push @acc, $acc;
    }

    # build phenotype hash for save_phenotypes()
    push @phenotypes, {
      id          => $ensg,
      type        => 'Gene',
      description => $disease,
      gencc_submitter      => $submitter,
      gencc_classification => $class,
      gencc_inherit_mode   => $moi,
      gencc_uuid           => $uuid,
      accessions           => \@acc,
    };
  }
  close $IN;

  $self->print_logFH("Parsed rows=$n_rows, skipped_submitter=$n_skip, unmatched_gene=$n_unmatched, to_save=".scalar(@phenotypes)."\n")
    if $self->param('debug_mode');

  return { phenotypes => \@phenotypes };
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
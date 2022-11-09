#!/usr/bin/perl
use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;
use Digest::MD5 qw(md5_hex);

my ($species, $port, $host, $user, $pass, $dbname,
    $peptide, $output_file, $model) = @ARGV;

# Extract model name
my $model_name;
if ( $model =~ /humdiv/i ) {
  $model_name = "humdiv";
} elsif ( $model =~ /humvar/i ) {
  $model_name = "humvar";
} else {
  die "ERROR: PolyPhen-2 data model not recognised based on filename $model\n";
}

# parse and store the results
my @fields;
my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
    -analysis           => 'polyphen',
    -sub_analysis       => $model_name,
    -peptide_length     => length($peptide),
    -translation_md5    => md5_hex($peptide),
);
my $any_results = 0;

open (RESULT, "<$output_file") or die "Failed to open output file: $!";
while (<RESULT>) {
  if (/^#/) {
    s/#//g;
    @fields = split /\s+/;
    next;
  }
  die "No header line in result file $output_file?" unless @fields; 

  my @values = split /\t/;

  # trim whitespace
  map { $_ =~ s/^\s+//; $_ =~ s/\s+$// } @values; 

  # parse the results into a hash
  my %results = map { $fields[$_] => $values[$_] } (0 .. @fields-1);
  my $alt_aa      = $results{o_aa2};
  my $prediction  = $results{prediction};
  my $prob        = $results{pph2_prob};
  my $position    = $results{o_pos};

  next unless $position && $alt_aa && $prob;
  $any_results = 1;
  
  $pred_matrix->add_prediction(
    $position,
    $alt_aa,
    $prediction, 
    $prob,
  );
}

# save the predictions to the database unless they are null matrices
if ( $any_results ){
  my $var_dba = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(
      '-species' => $species,
      '-port'    => $port,
      '-host'    => $host,
      '-user'    => $user,
      '-pass'    => $pass,
      '-dbname'  => $dbname
    );
  my $pfpma = $var_dba->get_ProteinFunctionPredictionMatrixAdaptor
      or die "Failed to get matrix adaptor";
  $pfpma->store($pred_matrix);
  $var_dba->dbc and $var_dba->dbc->disconnect_if_idle();
} else {
  warn "Skipping: no results to store\n";
}

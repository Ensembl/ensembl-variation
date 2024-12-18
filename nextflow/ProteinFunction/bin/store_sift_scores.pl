#!/usr/bin/perl
use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;
use Digest::MD5 qw(md5_hex);

my ($species, $port, $host, $user, $pass, $dbname,
    $offline, $sqlite,
    $peptide, $res_file) = @ARGV;

# parse the results file
my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
    -analysis           => 'sift',
    -peptide_length     => length($peptide),
    -translation_md5    => md5_hex($peptide),
);

# parse and store the results
my %evidence_stored;
my $results_available = 0;
open (RESULTS, "<$res_file") or die "Failed to open $res_file: $!";

while (<RESULTS>) {
  chomp;
  next if /WARNING/;
  next if /NOT SCORED/;

  my ($subst, $prediction, $score, $median_cons, $num_seqs, $blocks) = split;
  my ($ref_aa, $position, $alt_aa) = $subst =~ /([A-Z])(\d+)([A-Z])/;
  next unless $ref_aa && $alt_aa && defined $position;
  $results_available = 1;

  my $low_quality = $median_cons > 3.25 || $num_seqs < 10;

  $pred_matrix->add_prediction(
    $position,
    $alt_aa,
    $prediction, 
    $score,
    $low_quality
  );
  unless ($evidence_stored{$position} == 1) {
    ## add attribs by position
    $pred_matrix->add_evidence( 'sequence_number',  $position, $num_seqs );
    $pred_matrix->add_evidence( 'conservation_score', $position, $median_cons );
    $evidence_stored{$position} = 1;
  }
}

# save the predictions to the database
if ($results_available == 1 ){
  if (!$offline){
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
  }

  if ($sqlite){
    my $dbh = DBI->connect("dbi:SQLite:dbname=$db","","");
    my $sth = $dbh->prepare("INSERT INTO predictions VALUES(?, ?, ?)");
    $sth->execute($pred_matrix->translation_md5, 267, $pred_matrix->serialize)
  }
} else {
  warn "Skipping: no results to store\n";
}

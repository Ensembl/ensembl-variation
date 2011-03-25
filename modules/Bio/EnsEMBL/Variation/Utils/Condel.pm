package Bio::EnsEMBL::Variation::Utils::Condel;

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(get_condel_prediction);

use Bio::EnsEMBL::Variation::Utils::CondelConstants qw($CONDEL_CONFIG $CONDEL_SIFT_DATA $CONDEL_POLYPHEN_DATA);

sub get_condel_prediction {

    my ($original_sift_score, $original_polyphen_score) = @_;
    
    my %config      = %$CONDEL_CONFIG;
    my %sift        = %$CONDEL_SIFT_DATA;
    my %polyphen    = %$CONDEL_POLYPHEN_DATA;

    my %WAS;
    my %class;
    
    $WAS{'sift'} = $original_sift_score;
    $WAS{'pph2'} = $original_polyphen_score;
   
    my $base = 0;
    my $int_score = 0;
    my $sift_score = sprintf("%.3f", $WAS{'sift'});
    if ($sift_score <= $config{'cutoff.HumVar.sift'}){
      $int_score += sprintf("%.3f", (1 - $sift_score/$config{'max.HumVar.sift'})*(1-$sift{'tn'}{"$sift_score"}));
      $base += 1-$sift{'tn'}{"$sift_score"};
      $class{'sift'} = 'deleterious';
    }
    else {
      $int_score += sprintf("%.3f", (1 - $sift_score/$config{'max.HumVar.sift'})*(1-$sift{'tp'}{"$sift_score"}));
      $base += 1-$sift{'tp'}{"$sift_score"};
      $class{'sift'} = 'neutral';
    }
    my $polyphen_score = sprintf("%.3f", $WAS{'pph2'});
    if ($polyphen_score >= $config{'cutoff.HumVar.polyphen'}){
      $int_score += sprintf("%.3f", $polyphen_score/$config{'max.HumVar.polyphen'}*(1-$polyphen{'tn'}{"$polyphen_score"}));
      $base += 1-$polyphen{'tn'}{"$polyphen_score"};
      $class{'polyphen'} = 'deleterious';
    }
    else {
      $int_score += sprintf("%.3f", $polyphen_score/$config{'max.HumVar.polyphen'}*(1-$polyphen{'tp'}{"$polyphen_score"}));
      $base += 1-$polyphen{'tp'}{"$polyphen_score"};
      $class{'polyphen'} = 'neutral';
    }
    if ($base == 0){
      $int_score = -1;
      $class{'condel'} = 'not_computable_was';
    }
    else {
      $int_score = sprintf("%.3f", $int_score/$base);
    }

    if ($int_score >= 0.469){
      $class{'condel'} = 'deleterious';
    }
    elsif ($int_score > 0 && $int_score < 0.469) {
      $class{'condel'} = 'neutral';
    }
    
    my $score = $int_score;
    my $class = $class{'condel'};

    return wantarray ?  ($class, $score) : $class;
}

1;

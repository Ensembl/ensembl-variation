=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Variation::Utils::Condel

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::Utils::Condel qw(get_condel_prediction);
    
    my ($prediction, $score) = get_condel_prediction($sift_score, $polyphen_score);

    print "condel prediction: $prediction score: $score\n";

=head1 DESCRIPTION

This module provides a single subroutine get_condel_prediction which calculates the Condel
(Consensus Deleteriousness) weighted average score for a missense mutation that has both a 
SIFT and PolyPhen-2 score. Condel is developed by the Biomedical Genomics Group of the 
Universitat Pompeu Fabra, at the Barcelona Biomedical Research Park (bg.upf.edu/Condel).
The code in this module is based on a script provided by this group and slightly
reformatted to fit into the Ensembl API.

Various constants used by Condel can be found in Bio::EnsEMBL::Variation::Utils::CondelConstants.

=cut

package Bio::EnsEMBL::Variation::Utils::Condel;

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(get_condel_prediction);

use Bio::EnsEMBL::Variation::Utils::CondelConstants qw($CONDEL_CONFIG $CONDEL_SIFT_DATA $CONDEL_POLYPHEN_DATA);

our $USE_V2 = 1;

=head2 get_condel_prediction

  Arg[1]      : float SIFT score
  Arg[2]      : float PolyPhen score
  Example     : my ($prediction, $score) = get_condel_prediction($sift_score, $polyphen_score);
  Description : returns the Condel consensus prediction given SIFT and PolyPhen scores for a 
                missense mutation
  ReturnType  : if called in scalar context, the qualitative Condel prediction (as a string which
                will be either 'neutral' or 'deleterious'), if called in list context a 2 element
                list of the qualitative prediction and the Condel score as a float between 0 and 1

=cut

sub get_condel_prediction {

    my ($sift_score, $polyphen_score) = @_;
    
    my %config      = %$CONDEL_CONFIG;
    my %sift        = %$CONDEL_SIFT_DATA;
    my %polyphen    = %$CONDEL_POLYPHEN_DATA;

    my %class;
    
    my $base = 0;
    my $int_score = 0;

    $sift_score     = sprintf("%.3f", $sift_score);
    $polyphen_score = sprintf("%.3f", $polyphen_score);
    
    if ($sift_score <= $config{'cutoff.HumVar.sift'}){
      $int_score += sprintf("%.3f", (1 - $sift_score/$config{'max.HumVar.sift'})*(1-$sift{'tn'}{"$sift_score"}));
      $base += $USE_V2 ? 1 : 1-$sift{'tn'}{"$sift_score"};
      $class{'sift'} = 'deleterious';
    }
    else {
      $int_score += sprintf("%.3f", (1 - $sift_score/$config{'max.HumVar.sift'})*(1-$sift{'tp'}{"$sift_score"}));
      $base += $USE_V2 ? 1 : 1-$sift{'tp'}{"$sift_score"};
      $class{'sift'} = 'neutral';
    }
    
    if ($polyphen_score >= $config{'cutoff.HumVar.polyphen'}){
      $int_score += sprintf("%.3f", $polyphen_score/$config{'max.HumVar.polyphen'}*(1-$polyphen{'tn'}{"$polyphen_score"}));
      $base += $USE_V2 ? 1 : 1-$polyphen{'tn'}{"$polyphen_score"};
      $class{'polyphen'} = 'deleterious';
    }
    else {
      $int_score += sprintf("%.3f", $polyphen_score/$config{'max.HumVar.polyphen'}*(1-$polyphen{'tp'}{"$polyphen_score"}));
      $base += $USE_V2 ? 1 : 1-$polyphen{'tp'}{"$polyphen_score"};
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
    elsif ($int_score >= 0 && $int_score < 0.469) {
      $class{'condel'} = 'neutral';
    }

    # if the user wants an array, return the class and score, otherwise just return the class
        
    return wantarray ?  ($class{'condel'}, $int_score) : $class{'condel'};
}

1;

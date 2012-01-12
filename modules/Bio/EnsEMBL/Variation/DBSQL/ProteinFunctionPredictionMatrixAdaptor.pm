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

Bio::EnsEMBL::Variation::DBSQL::ProteinFunctionPredictionMatrixAdaptor

=head1 DESCRIPTION

This adaptor lets you fetch compressed binary formatted protein
function prediction matrices from the variation databases.

=cut
 
use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::ProteinFunctionPredictionMatrixAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;

use DBI qw(:sql_types);

use base qw(Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor);

sub new_fake {
    my $class = shift;
    my $species = shift;

    my $self = bless {}, $class;

    return $self;
}

=head2 fetch_by_analysis_transcript_stable_id

  Arg [1]    : string $analysis - the name of the prediction tool
  Arg [2]    : string $transcript_stable_id - the stable id of the transcript
  Description: Fetch the prediction matrix for the given tool and transcript
  Returntype : Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix instance or undef
  Status     : At risk

=cut

sub fetch_by_analysis_transcript_stable_id {
    my ($self, $analysis, $transcript_stable_id) = @_;

    $self->{_cols} = [$analysis.'_predictions'];

    $self->{_analysis} = $analysis;
    
    my $constraint = "transcript_stable_id = ?";
    
    $self->bind_param_generic_fetch($transcript_stable_id, SQL_VARCHAR);

    my ($matrix) = @{ $self->generic_fetch($constraint) };
    
    return $matrix;
}

=head2 fetch_by_analysis_translation_stable_id

  Arg [1]    : string $analysis - the name of the prediction tool
  Arg [2]    : string $translation_stable_id - the stable id of the translation
  Description: Fetch the prediction matrix for the given tool and translation
  Returntype : Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix instance or undef
  Status     : At risk

=cut

sub fetch_by_analysis_translation_stable_id {
    my ($self, $analysis, $translation_stable_id) = @_;

    $self->{_cols} = [$analysis.'_predictions'];

    $self->{_analysis} = $analysis;
    
    my $constraint = "translation_stable_id = ?";
    
    $self->bind_param_generic_fetch($translation_stable_id, SQL_VARCHAR);

    my ($matrix) = @{ $self->generic_fetch($constraint) };
    
    return $matrix;
}

=head2 fetch_sift_predictions_by_transcript_stable_id

  Arg [1]    : string $transcript_stable_id - the stable id of the transcript 
  Description: Fetch the sift prediction matrix for the given transcript
  Returntype : Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix instance or undef
  Status     : At risk

=cut

sub fetch_sift_predictions_by_transcript_stable_id {
    my ($self, $transcript_stable_id) = @_;
    return $self->fetch_by_analysis_transcript_stable_id('sift', $transcript_stable_id);
}

=head2 fetch_polyphen_predictions_by_transcript_stable_id

  Arg [1]    : string $transcript_stable_id - the stable id of the transcript 
  Description: Fetch the polyphen prediction matrix for the given transcript
  Returntype : Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix instance or undef
  Status     : At risk

=cut

sub fetch_polyphen_predictions_by_transcript_stable_id {
    my ($self, $transcript_stable_id) = @_;
    return $self->fetch_by_analysis_transcript_stable_id('polyphen', $transcript_stable_id);
}

=head2 fetch_sift_predictions_by_translation_stable_id

  Arg [1]    : string $translation_stable_id - the stable id of the translation 
  Description: Fetch the sift prediction matrix for the given translation
  Returntype : Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix instance or undef
  Status     : At risk

=cut

sub fetch_sift_predictions_by_translation_stable_id {
    my ($self, $translation_stable_id) = @_;
    return $self->fetch_by_analysis_translation_stable_id('sift', $translation_stable_id);
}

=head2 fetch_polyphen_predictions_by_translation_stable_id

  Arg [1]    : string $translation_stable_id - the stable id of the translation
  Description: Fetch the polyphen prediction matrix for the given translation
  Returntype : Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix instance or undef
  Status     : At risk

=cut

sub fetch_polyphen_predictions_by_translation_stable_id {
    my ($self, $translation_stable_id) = @_;
    return $self->fetch_by_analysis_translation_stable_id('polyphen', $translation_stable_id);
}

sub _columns {
    my $self = shift;
    return @{ $self->{_cols} };
}

sub _tables {
    return (['protein_function_predictions']);
}

sub _objs_from_sth {

    my ($self, $sth) = @_;

    my $matrix_string;

    $sth->bind_columns(\$matrix_string);

    my @matrices;

    while ($sth->fetch) {
        if ($matrix_string) {
            push @matrices, Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
                -analysis   => $self->{_analysis},
                -matrix     => $matrix_string
            );
        }
    }

    return \@matrices;
}

1;


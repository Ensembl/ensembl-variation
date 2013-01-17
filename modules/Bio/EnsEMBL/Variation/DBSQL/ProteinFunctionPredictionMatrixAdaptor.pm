=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

This adaptor lets you store and fetch compressed binary formatted protein
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

=head2 store

  Arg [1]    : ProteinFunctionPredictionMatrix $matrix - the matrix you want to store
  Description: Store the given matrix in the database
  Status     : Stable

=cut

sub store {
    my ($self, $matrix) = @_;

    # get our analysis attrib ID from the attrib table

    throw("You need to supply a translation MD5 to store a matrix in the database")
        unless $matrix->translation_md5;

    my $analysis = $matrix->analysis;

    $analysis .= '_'.$matrix->sub_analysis if defined $matrix->sub_analysis;

    my $analysis_attrib_id = $self->db->get_AttributeAdaptor->attrib_id_for_type_value(
        'prot_func_analysis', 
        $analysis
    );

    throw("No attrib_id for analysis $analysis?") unless defined $analysis_attrib_id;

    my $dbh = $self->dbc->db_handle;

    # first add the MD5 to the translation_md5 table if necessary

    my $md5_sth = $dbh->prepare(qq{INSERT IGNORE INTO translation_md5 (translation_md5) VALUES (?)});

    $md5_sth->execute($matrix->translation_md5);

    # then add the matrix

    my $matrix_sth = $dbh->prepare(qq{
        INSERT INTO protein_function_predictions (translation_md5_id, analysis_attrib_id, prediction_matrix) 
        VALUES ((SELECT translation_md5_id FROM translation_md5 WHERE translation_md5 = ?),?,?)
    });
 
    $matrix_sth->execute($matrix->translation_md5, $analysis_attrib_id, $matrix->serialize);
}

=head2 fetch_by_analysis_translation_md5

  Arg [1]    : string $analysis - the name of the prediction tool
  Arg [2]    : string $translation_md5 - the hex MD5 hash of the translation sequence
  Description: Fetch the prediction matrix for the given tool and peptide sequence MD5
  Returntype : Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix instance or undef
  Status     : Stable

=cut

sub fetch_by_analysis_translation_md5 {
    my ($self, $analysis, $translation_md5) = @_;

    my $constraint = "t.translation_md5 = ? AND a.value = ?";
    
    $self->bind_param_generic_fetch($translation_md5, SQL_VARCHAR);
    $self->bind_param_generic_fetch($analysis, SQL_VARCHAR);

    my ($matrix) = @{ $self->generic_fetch($constraint) };
    
    return $matrix;
}

=head2 fetch_polyphen_predictions_by_translation_md5

  Arg [1]    : string $translation_md5 - the hex MD5 hash of the translation sequence
  Arg [2]    : string $model - the desired classifier model, either 'humvar' or 'humdiv', 
               the default is 'humvar'
  Description: Fetch the polyphen prediction matrix for the given translation sequence MD5
               and classifier model
  Returntype : Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix instance or undef
  Status     : Stable

=cut

sub fetch_polyphen_predictions_by_translation_md5 {
    my ($self, $translation_md5, $model) = @_;

    $model ||= 'humvar';

    $model = lc($model);

    throw("Unrecognised model for PolyPhen: '$model'") 
        unless (($model eq 'humvar') || ($model eq 'humdiv'));
    
    return $self->fetch_by_analysis_translation_md5('polyphen_'.$model, $translation_md5);
}

=head2 fetch_sift_predictions_by_translation_md5

  Arg [1]    : string $translation_md5 - the hex MD5 hash of the translation sequence
  Description: Fetch the sift prediction matrix for the given translation sequence MD5
  Returntype : Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix instance or undef
  Status     : Stable

=cut

sub fetch_sift_predictions_by_translation_md5 {
    my ($self, $translation_md5) = @_;
    return $self->fetch_by_analysis_translation_md5('sift', $translation_md5);
}

sub _columns {
    return qw(t.translation_md5 a.value p.prediction_matrix);
}

sub _tables {
    return (
        ['protein_function_predictions', 'p'],
        ['translation_md5', 't'],
        ['attrib', 'a']
    );
}

sub _default_where_clause {
    return join ' AND ', (
        'p.translation_md5_id = t.translation_md5_id',
        'p.analysis_attrib_id = a.attrib_id' 
    );
}

sub _objs_from_sth {

    my ($self, $sth) = @_;

    my $md5;
    my $analysis;
    my $matrix;

    $sth->bind_columns(\$md5, \$analysis, \$matrix);

    my @matrices;

    while ($sth->fetch) {
        if ($matrix) {
            my ($super_analysis, $sub_analysis) = split /_/, $analysis;
            
            push @matrices, Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
                -translation_md5    => $md5,
                -analysis           => $super_analysis,
                -sub_analysis       => $sub_analysis,
                -matrix             => $matrix,
            );
        }
    }

    return \@matrices;
}

1;


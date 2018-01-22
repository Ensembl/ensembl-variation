
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

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


   ## store attribs

   my $trans_id_sth = $dbh->prepare(qq{ SELECT translation_md5_id 
                                        FROM translation_md5 
                                        WHERE translation_md5 = ?
                                      });
    $trans_id_sth->execute($matrix->translation_md5);
    my $trans_id = $trans_id_sth->fetchall_arrayref();


   my $attrib_sth = $dbh->prepare(qq{
        INSERT INTO protein_function_predictions_attrib 
         (translation_md5_id, analysis_attrib_id, attrib_type_id ,  position_values)
        VALUES (?,?,?,?)
    });

    my %attrib_id;

    my %attribs;
    ##  compress position-value data by evidence type atrib 
    foreach my $ev ( @{$matrix->{evidence}} ){   ## array of arrays [type  pos value]
	$attrib_id{$ev->[0]}  = $self->db->get_AttributeAdaptor->attrib_id_for_type_code($ev->[0])
	    unless defined $attrib_id{$ev->[0]} ;
	die "No attrib available for evidence type $ev->[0]\n" unless defined $attrib_id{$ev->[0]} ;

          ## store conservation score as int not float
	  $ev->[2] = $ev->[2]*100 if $ev->[0] eq 'conservation_score';

	$attribs{$ev->[0]} .=  pack("ww", $ev->[1],  $ev->[2]);
    }

    ## compressed by type - now store
    foreach my $evidence_attrib(keys %attribs){

        $attrib_sth->execute($trans_id->[0]->[0], 
                             $analysis_attrib_id, 
                             $attrib_id{ $evidence_attrib }, 
                             $attribs{$evidence_attrib}
                            );
    }
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

sub fetch_evidence_for_prediction{
   my ($self, $translation_md5, $analysis_type) = @_;

   my %evidence;

   my $dbh = $self->dbc->db_handle;

   ## look up evidence data
   my $evidence_sth = $dbh->prepare(qq{ SELECT attrib_type.code, 
                                               pfpa.position_values
                                        from protein_function_predictions_attrib pfpa, 
                                             attrib, 
                                             attrib_type, 
                                             translation_md5
                                        where translation_md5.translation_md5 = ?
                                        and pfpa.translation_md5_id = translation_md5.translation_md5_id
                                        and pfpa.analysis_attrib_id = attrib.attrib_id
                                        and attrib.value = ?
                                        and pfpa.attrib_type_id = attrib_type.attrib_type_id
                                        });

 
   $evidence_sth->execute($translation_md5, $analysis_type );
   my $ev_data = $evidence_sth->fetchall_arrayref();
  
   foreach my $evidence(@{$ev_data}){

       ## unpack position- value data
       my @attribs = unpack("(ww)*", $evidence->[1] );
       my %data;
       
       while(@attribs) {
          my $position = shift @attribs;
          my $value    = shift @attribs;
          ## stored conservation score as int not float
	  $value  = $value /100 if $evidence->[0] eq 'conservation_score';
	  $evidence{$evidence->[0]}{$position} = $value;
       }
   
   }
   return \%evidence;

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
            
           my $matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
                -translation_md5    => $md5,
                -analysis           => $super_analysis,
                -sub_analysis       => $sub_analysis,
                -matrix             => $matrix,
                -adaptor            => $self
            );


	    push @matrices, $matrix;
        }
    }

    return \@matrices;
}

1;


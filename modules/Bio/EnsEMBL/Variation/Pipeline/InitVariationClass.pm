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
package Bio::EnsEMBL::Variation::Pipeline::InitVariationClass;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use POSIX qw(ceil);

my $DEBUG = 0;

sub fetch_input {
    
    my $self = shift;

    my $num_chunks = $self->required_param('num_chunks');
    
    my $var_dba = $self->get_species_adaptor('variation');
        
    my $aa = $var_dba->get_AttributeAdaptor;
    
    my $dbc = $var_dba->dbc();
    
    # first set everything in variation (except HGMDs) to 'sequence_alteration' by default
    # because sometimes we miss them because there is no variation_feature
    # or any alleles (though this should become unnecessary as we move to the
    # new approach to failing for all species)

    my $default_attrib_id = $aa->attrib_id_for_type_value('SO_term', 'sequence_alteration');

    die "No attrib_id for 'sequence_alteration'" unless defined $default_attrib_id;

    $dbc->do(qq{
        UPDATE  variation v, source s
        SET     v.class_attrib_id = $default_attrib_id
        WHERE   v.source_id = s.source_id
        AND     s.name != 'HGMD-PUBLIC'
    });
    
    # now create some temp tables to store the class attribs

    my $temp_var_table = 'temp_variation_class';
    my $temp_var_feat_table = 'temp_variation_feature_class';
    
    $dbc->do(qq{DROP TABLE IF EXISTS $temp_var_table});
    $dbc->do(qq{DROP TABLE IF EXISTS $temp_var_feat_table});
    
    $dbc->do(qq{CREATE TABLE $temp_var_table LIKE variation});
    $dbc->do(qq{CREATE TABLE $temp_var_feat_table LIKE variation_feature});

    $dbc->do(qq{ALTER TABLE $temp_var_table DISABLE KEYS});
    $dbc->do(qq{ALTER TABLE $temp_var_feat_table DISABLE KEYS});

    # now get bins of variation_ids to update

    my $get_max_var_id_sth = $dbc->prepare(qq{
        SELECT max(variation_id) FROM variation
    });

    $get_max_var_id_sth->execute;
    my $max = $get_max_var_id_sth->fetchall_arrayref();

    my $chunk_size = ceil($max->[0]->[0] / $num_chunks);

    my @output_ids;
    my $start = 1;    
    while ($start < $max->[0]->[0]){
      my $stop =  $start + $chunk_size -1;
      push @output_ids, {
            variation_id_start  => $start,
            variation_id_stop   => $stop,
            temp_var_table      => $temp_var_table,
            temp_var_feat_table => $temp_var_feat_table,
        };
      $start = $stop + 1;
    }

    $self->param('chunk_output_ids', \@output_ids);

    $self->param(
        'finish_var_class', [{
            temp_var_table      => $temp_var_table,
            temp_var_feat_table => $temp_var_feat_table,
        }]
    );
}

sub write_output {
    my $self = shift;

    $self->dataflow_output_id($self->param('finish_var_class'), 1);
    $self->dataflow_output_id($self->param('chunk_output_ids'), 2);
}

1;


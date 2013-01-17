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

    # now get an ordered list of all the variation_ids
        
    my $get_var_ids_sth = $dbc->prepare(qq{
        SELECT variation_id FROM variation ORDER BY variation_id
    });

    $get_var_ids_sth->execute;

    my @var_ids;

    while (my ($var_id) = $get_var_ids_sth->fetchrow_array) {
        push @var_ids, $var_id;
    }

    # and split them up into as many chunks as requested

    my $num_vars = scalar @var_ids;

    my $chunk_size = ceil($num_vars / $num_chunks);

    my @output_ids;

    while (@var_ids) {

        my $start = $var_ids[0];
        my $stop  = $chunk_size <= $#var_ids ? $var_ids[$chunk_size - 1] : $var_ids[$#var_ids];

        splice(@var_ids, 0, $chunk_size);

        push @output_ids, {
            variation_id_start  => $start,
            variation_id_stop   => $stop,
            temp_var_table      => $temp_var_table,
            temp_var_feat_table => $temp_var_feat_table,
        };
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


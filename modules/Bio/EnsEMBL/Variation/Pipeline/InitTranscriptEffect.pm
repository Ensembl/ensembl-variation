=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

package Bio::EnsEMBL::Variation::Pipeline::InitTranscriptEffect;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub fetch_input {
   
    my $self = shift;

    my $core_dba = $self->get_species_adaptor('core');
    my $var_dba = $self->get_species_adaptor('variation');
    
    my %default_params = (
        ensembl_registry    => $self->param('ensembl_registry'),
        species             => $self->param('species'),
    );

    my $dbc = $var_dba->dbc();

    # truncate the table because we don't want duplicates

    $dbc->do("TRUNCATE TABLE transcript_variation");
   
    # disable the indexes on the table we're going to insert into as
    # this significantly speeds up the TranscriptEffect process

    $dbc->do("ALTER TABLE transcript_variation DISABLE KEYS");
    
    my $ga = $core_dba->get_GeneAdaptor
        or die "Failed to get slice adaptor";

    my @transcript_output_ids;
    
    my $gene_count = 0;

    for my $gene (@{ $ga->fetch_all }) {
        $gene_count++;
        
        for my $transcript (@{ $gene->get_all_Transcripts }) {

            push @transcript_output_ids, {
                transcript_stable_id  => $transcript->stable_id,
                %default_params
            };
        }            
        if ($DEBUG) {
            last if $gene_count >= 100;
        }
    }

    $self->param('transcript_output_ids', \@transcript_output_ids);

    $self->param(
        'rebuild_indexes', [{
            %default_params,
            tables              => ['transcript_variation'],
        }]
    );
    
    $self->param(
        'update_vf', [{%default_params}]
    );
    
    $self->param(
        'set_var_class', [{%default_params}]
    );

}

sub write_output {
    my $self = shift;

    $self->dataflow_output_id($self->param('rebuild_indexes'), 1);
    $self->dataflow_output_id($self->param('update_vf'), 2);
    $self->dataflow_output_id($self->param('set_var_class'), 3);
    $self->dataflow_output_id($self->param('transcript_output_ids'), 4);
}

1;

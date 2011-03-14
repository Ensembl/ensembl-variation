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

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::EnsEMBL::Hive::Process');

my $DEBUG = 0;

sub fetch_input {
   
    my $self = shift;

    my $reg_file = $self->param('ensembl_registry')
        or die "ensembl_registry is a required parameter";

    my $species = $self->param('species')
        or die "species is a required parameter";
    
    my $reg = 'Bio::EnsEMBL::Registry';
    
    $reg->load_all($reg_file, 0, 1);
  
    my $core_dba = $reg->get_DBAdaptor($species, 'core')
        or die "failed to get core DBA for $species";
   
    my $var_dba = $reg->get_DBAdaptor($species, 'variation')
        or die "failed to get variation DBA for $species";

    my $dbh = $var_dba->dbc->db_handle;

    # truncate the table because we don't want duplicates

    $dbh->do("TRUNCATE TABLE transcript_variation");
   
    # disable the indexes on the table we're going to insert into as
    # this significantly speeds up the TranscriptEffect process

    $dbh->do("ALTER TABLE transcript_variation DISABLE KEYS")
        or warn "Failed to disable keys on transcript_variation: ".$dbh->errstr;
    
    my $ga = $core_dba->get_GeneAdaptor
        or die "Failed to get slice adaptor";

    my @transcript_output_ids;
    
    my $gene_count = 0;

    for my $gene (@{ $ga->fetch_all }) {
        $gene_count++;
        
        for my $transcript (@{ $gene->get_all_Transcripts }) {

            push @transcript_output_ids, {
                transcript_stable_id    => $transcript->stable_id,
                ensembl_registry        => $reg_file,
                species                 => $species,
            };
        }            
        if ($DEBUG) {
            last if $gene_count >= 100;
        }
    }

    $self->param('transcript_output_ids', \@transcript_output_ids);

    $self->param(
        'rebuild_indexes', [{
            ensembl_registry    => $reg_file, 
            species             => $species, 
            tables              => ['transcript_variation'],
        }]
    );
    
    $self->param(
        'update_vf', [{
            ensembl_registry    => $reg_file, 
            species             => $species, 
        }]
    );
    
    $self->param(
        'set_var_class', [{
            ensembl_registry    => $reg_file, 
            species             => $species, 
        }]
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

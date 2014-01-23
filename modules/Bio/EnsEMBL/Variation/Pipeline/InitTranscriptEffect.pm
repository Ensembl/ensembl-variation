=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::Pipeline::InitTranscriptEffect;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub fetch_input {
   
    my $self = shift;

    my $include_lrg = $self->param('include_lrg');

    my $core_dba = $self->get_species_adaptor('core');
    my $var_dba = $self->get_species_adaptor('variation');
    
    my $dbc = $var_dba->dbc();


    # check for out of date seq_regions in variation database
    my $sequences_ok = $self->check_seq_region();
    die "Seq region ids are not compatible\n" unless $sequences_ok == 1;


    my $ga = $core_dba->get_GeneAdaptor or die "Failed to get gene adaptor";

    my @transcript_output_ids;

    
    my $gene_count = 0;

    # fetch all the regular genes

    my @genes = @{ $ga->fetch_all };

    if ($include_lrg) {
        # fetch the LRG genes as well
        
        push @genes, @{ $ga->fetch_all_by_biotype('LRG_gene') }
    }

    for my $gene (@genes) {
        $gene_count++;
        
        for my $transcript (@{ $gene->get_all_Transcripts }) {

            push @transcript_output_ids, {
                transcript_stable_id  => $transcript->stable_id,
            };
        }            
        if ($DEBUG) {
            last if $gene_count >= 100;
        }
    }



    if (@transcript_output_ids) {
        
        # check we actually found some transcripts

        # truncate the table because we don't want duplicates

        $dbc->do("TRUNCATE TABLE transcript_variation");

        # disable the indexes on the table we're going to insert into as
        # this significantly speeds up the TranscriptEffect process

        $dbc->do("ALTER TABLE transcript_variation DISABLE KEYS");

        # create table to use for web index generation
        # truncate tables incase TranscriptVariation is being updated for a pre-existing database

        $dbc->do("create table if not exists variation_hgvs(
                  variation_id int(10) unsigned not null,
                  hgvs_name varchar(255),
                  primary key(variation_id, hgvs_name)) ");

        $dbc->do("create table if not exists variation_genename (
                  variation_id int(10), 
                  gene_name varchar(255) )" );

        $dbc->do("TRUNCATE TABLE variation_hgvs");
        $dbc->do("TRUNCATE TABLE variation_genename ");


        $self->param('transcript_output_ids', \@transcript_output_ids);

        $self->param(
            'rebuild_indexes', [{
                tables => ['transcript_variation'],
            }]
        );

        # we need to kick off the update_vf analysis as well, 
        # but it doesn't have any parameters we need to set here

        $self->param(
            'update_vf', [{}]
        );
    }
}

sub write_output {
    my $self = shift;
    
    if (my $transcript_output_ids = $self->param('transcript_output_ids')) {
        $self->dataflow_output_id($self->param('rebuild_indexes'), 2);
        $self->dataflow_output_id($self->param('update_vf'), 3);
        $self->dataflow_output_id($transcript_output_ids, 4);
        $self->dataflow_output_id($self->param('check_transcript_variation'), 5);
    }

    return;
}

## check for out of date seq_regions in variation database
##    human patches can change between releases
##    such differences break TranscriptEffect
sub check_seq_region{

   my $self = shift;
   
   my $stmt = qq[ select seq_region_id, name from seq_region];
   
   my $core_dba = $self->get_species_adaptor('core');
   my $var_dba  = $self->get_species_adaptor('variation');
 
   my $core_seq_sth = $core_dba->dbc->prepare($stmt);
   $core_seq_sth->execute();
   my $core_ids = $core_seq_sth->fetchall_arrayref();

   my %expected_ids;
   foreach my $l(@{$core_ids}){
       $expected_ids{$l->[0]} = $l->[1];
   }

   my $var_seq_sth = $var_dba->dbc->prepare($stmt);
   $var_seq_sth->execute();
   my $var_ids = $var_seq_sth->fetchall_arrayref();

   my $OK = 1;
   foreach my $l(@{$var_ids}){
       unless (defined $expected_ids{$l->[0]} ){
           $self->warning( 'Seq_region_id in variation db is not in core: '.  $l->[0]. ' ' . $l->[1]);
           $OK = 0;
       }
   }    
   return $OK;
}

1;

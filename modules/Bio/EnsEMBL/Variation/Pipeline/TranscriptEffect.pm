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

package Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect;

use strict;

use Bio::EnsEMBL::Registry;

use Bio::EnsEMBL::Variation::TranscriptVariation;

use base ('Bio::EnsEMBL::Hive::ProcessWithParams');

my $DEBUG   = 0;

sub run {
    my $self = shift;

    my $transcript_id = $self->param('transcript_stable_id') 
        or die "transcript_stable_id is a required parameter";

    my $reg_file = $self->param('ensembl_registry')
        or die "ensembl_registry is a required parameter";
    
    my $species = $self->param('species')
        or die "species is a required parameter";

    my $disambiguate_sn_alleles = 
        $self->param('disambiguate_single_nucleotide_alleles'); 
 
    my $reg = 'Bio::EnsEMBL::Registry';
    
    $reg->load_all($reg_file, 0, 1);
    
    my $var_dba = $reg->get_DBAdaptor($species, 'variation')
        or die "failed to get variation DBA for $species";

    my $ta = $reg->get_adaptor($species, 'core', 'Transcript') 
        or die "failed to get transcript adaptor";

    my $sa = $reg->get_adaptor($species, 'core', 'Slice')
        or die "failed to get slice adaptor";

    my $tva = $reg->get_adaptor($species, 'variation', 'TranscriptVariation')
        or die "failed to get transcript variation adaptor";

    my $transcript = $ta->fetch_by_stable_id($transcript_id) 
        or die "failed to fetch transcript for stable id: $transcript_id";

    # we need to include failed variations

    $tva->db->include_failed_variations(1);

    my $slice = $sa->fetch_by_transcript_stable_id($transcript->stable_id, 5000)
        or die "failed to get slice around transcript: ".$transcript->stable_id;

    for my $vf ( @{ $slice->get_all_VariationFeatures }, 
                 @{ $slice->get_all_somatic_VariationFeatures } ) {

        my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
            -feature                 => $transcript,
            -variation_feature       => $vf,
            -adaptor                 => $tva,
            -disambiguate_single_nucleotide_alleles => $disambiguate_sn_alleles,
        );

        # if the variation has no effect on the transcript $tv will be undef

        if ($tv && ( scalar(@{ $tv->consequence_type }) > 0) ) {
            $tva->store($tv);
        }
    }
}

sub write_output {
    my $self = shift;
}

1;

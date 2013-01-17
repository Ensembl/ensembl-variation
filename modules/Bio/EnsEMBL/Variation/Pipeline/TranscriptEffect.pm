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

package Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::TranscriptVariation;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG   = 0;

sub run {
    my $self = shift;

    my $transcript_id = $self->required_param('transcript_stable_id'); 

    my $disambiguate_sn_alleles = 
        $self->param('disambiguate_single_nucleotide_alleles'); 
    
    my $variations_to_include;
    
    if (my $vars = $self->param('variations_to_include')) {
        # turn the list of variation names into a hash to speed up checking
        $variations_to_include = { map { $_ => 1 } @$vars };
    }

    my $core_dba = $self->get_species_adaptor('core');
    my $var_dba = $self->get_species_adaptor('variation');
    
    my $ta = $core_dba->get_TranscriptAdaptor;
    my $sa = $core_dba->get_SliceAdaptor;
    
    my $tva = $var_dba->get_TranscriptVariationAdaptor;

    my $transcript = $ta->fetch_by_stable_id($transcript_id) 
        or die "failed to fetch transcript for stable id: $transcript_id";

    # we need to include failed variations

    $tva->db->include_failed_variations(1);

    my $slice = $sa->fetch_by_transcript_stable_id(
        $transcript->stable_id, 
        MAX_DISTANCE_FROM_TRANSCRIPT
    ) or die "failed to get slice around transcript: ".$transcript->stable_id;

    for my $vf ( @{ $slice->get_all_VariationFeatures }, 
                 @{ $slice->get_all_somatic_VariationFeatures } ) {

        if (defined $variations_to_include) {
            next unless $variations_to_include->{$vf->variation_name};
        }

        my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
            -transcript         => $transcript,
            -variation_feature  => $vf,
            -adaptor            => $tva,
            -disambiguate_single_nucleotide_alleles => $disambiguate_sn_alleles,
        );

        # if the variation has no effect on the transcript $tv will be undef

        if ($tv && ( scalar(@{ $tv->consequence_type }) > 0) ) {
            $tva->store($tv);
        }
    }

    return;
}

sub write_output {
    my $self = shift;
}

1;

package Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect;

use strict;

use Bio::EnsEMBL::Registry;

use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(transcript_effect);

use base ('Bio::EnsEMBL::Hive::ProcessWithParams');

my $DEBUG   = 0;

sub run {
    my $self = shift;

    my $gene_id = $self->param('gene_stable_id') 
        or die "gene_stable_id is a required parameter";

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
  
    my $ga = $core_dba->get_GeneAdaptor 
        or die "Failed to get gene adaptor";

    my $sa = $core_dba->get_SliceAdaptor 
        or die "Failed to get slice adaptor";

    my $tva = $var_dba->get_TranscriptVariationNewAdaptor
        or die "Failed to get transcript variation adaptor";

    my $oca = $reg->get_adaptor('multi', 'variation', 'OverlapConsequence')
        or die "Failed to get consequence adaptor";
    
    my $cons = $oca->fetch_all or die "Failed to get consequences";

    my $gene = $ga->fetch_by_stable_id($gene_id) 
        or die "Failed to fetch gene for stable id: $gene_id";

    my $slice = $sa->fetch_by_gene_stable_id($gene_id, 5000)
        or die "Failed to get slice around gene: $gene_id";

    my @transcripts = @{ $gene->get_all_Transcripts };
    
    for my $tran (@transcripts) {
    
        # we have to transfer the transcript to the slice the
        # vfs live on or everything goes wrong...
        
        $tran = $tran->transfer($slice);

        for my $vf ( @{ $slice->get_all_VariationFeatures }, 
            @{ $slice->get_all_somatic_VariationFeatures} ) {

            my $worst_conseq_rank = 1000;
            my @worst_conseqs;

            if (my $tv = transcript_effect($vf, $tran, $cons)) {

                $tva->store($tv);
                
                # TODO: update variation feature table with worst consequence

#                for my $allele (@{ $tv->alt_alleles }) {
#                    for my $conseq (@{ $allele->consequences }) {
#                        if ($conseq->rank < $worst_conseq_rank) {
#                            $worst_conseqs{$conseq->SO_term}++;
#                            $worst_conseq_rank = $conseq->rank;
#                        }
#                    }
#                }
            }
        }  
    }

}

sub write_output {
    my $self = shift;
}

1;

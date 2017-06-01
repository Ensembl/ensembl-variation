
# change the transcript_variation & variation_feature tables to add new consequence protein_altering_variant
##################

## VARIATION_FEATURE
ALTER TABLE variation_feature CHANGE consequence_types consequence_types SET (
    'intergenic_variant',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_lost',
    'coding_sequence_variant',
    'missense_variant',
    'stop_gained',
    'synonymous_variant',
    'frameshift_variant',
    'non_coding_transcript_variant',
    'non_coding_transcript_exon_variant',
    'mature_miRNA_variant',
    'NMD_transcript_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'incomplete_terminal_codon_variant',
    'intron_variant',
    'splice_region_variant',
    'downstream_gene_variant',
    'upstream_gene_variant',
    'start_lost',
    'stop_retained_variant',
    'inframe_insertion',
    'inframe_deletion',
    'transcript_ablation',
    'transcript_fusion',
    'transcript_amplification',
    'transcript_translocation',
    'TFBS_ablation',
    'TFBS_fusion',
    'TFBS_amplification',
    'TFBS_translocation',
    'regulatory_region_ablation',
    'regulatory_region_fusion',
    'regulatory_region_amplification',
    'regulatory_region_translocation',
    'feature_elongation',
    'feature_truncation',
    'regulatory_region_variant',
    'TF_binding_site_variant',
    'protein_altering_variant',
    'start_retained_variant'
) DEFAULT 'intergenic_variant' NOT NULL;


## TRANSCRIPT VARIATION
ALTER TABLE transcript_variation CHANGE consequence_types consequence_types SET (
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_lost',
    'coding_sequence_variant',
    'missense_variant',
    'stop_gained',
    'synonymous_variant',
    'frameshift_variant',
    'non_coding_transcript_variant',
    'non_coding_transcript_exon_variant',
    'mature_miRNA_variant',
    'NMD_transcript_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'incomplete_terminal_codon_variant',
    'intron_variant',
    'splice_region_variant',
    'downstream_gene_variant',
    'upstream_gene_variant',
    'start_lost',
    'stop_retained_variant',
    'inframe_insertion',
    'inframe_deletion', 
    'transcript_ablation',
    'transcript_fusion',
    'transcript_amplification',
    'transcript_translocation',
    'TFBS_ablation',
    'TFBS_fusion',
    'TFBS_amplification',
    'TFBS_translocation',
    'regulatory_region_ablation',
    'regulatory_region_fusion',
    'regulatory_region_amplification',
    'regulatory_region_translocation',
    'feature_elongation',
    'feature_truncation',
    'protein_altering_variant',
    'start_retained_variant'
);


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_89_90_b.sql|add start_retained_variant to consequence_types in variation_feature and transcript_variation');


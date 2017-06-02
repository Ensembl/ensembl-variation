
# Remove regulatory and TFBS consequences from consequence_types in transcript_variation table
##################

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
    'feature_elongation',
    'feature_truncation',
    'protein_altering_variant',
    'start_retained_variant'
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_89_90_c.sql|remove regulatory and TFBS consequences from consequence_types in  transcript_variation');


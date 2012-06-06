# rename consequence_type to consequence_types in variation_feature
ALTER TABLE variation_feature CHANGE consequence_type consequence_types SET('intergenic_variant','splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','nc_transcript_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','initiator_codon_variant','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation')  NOT NULL  DEFAULT 'intergenic_variant';

# add index
CREATE INDEX consequence_type_idx ON variation_feature(consequence_types);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_67_68_h.sql|rename consequence_type to consequence_types in variation_feature and add index');
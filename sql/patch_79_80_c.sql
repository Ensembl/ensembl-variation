
# change the transcript_variation & variation_feature tables to add new consequence protein_altering_variant
##################

## VARIATION_FEATURE
# create tmp table
CREATE TABLE tmp_vf LIKE variation_feature;

# alter its consequence type field
ALTER TABLE tmp_vf CHANGE consequence_types consequence_types SET (
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
        'protein_altering_variant'
    ) DEFAULT 'intergenic_variant' NOT NULL;

# insert into tmp table, replacing old terms with new as we go
INSERT INTO tmp_vf SELECT variation_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, variation_id, allele_string, variation_name, map_weight, flags, source_id, validation_status,
    REPLACE(consequence_types, 'initiator_codon_variant', 'start_lost'),
variation_set_id, class_attrib_id, somatic, minor_allele, minor_allele_freq, minor_allele_count, alignment_quality, evidence_attribs, clinical_significance, display
FROM variation_feature;

DROP TABLE variation_feature;
RENAME TABLE tmp_vf TO variation_feature;


## TRANSCRIPT VARIATION
# create tmp table
CREATE TABLE tmp_tv LIKE transcript_variation;

# alter its consequence type field
ALTER TABLE tmp_tv CHANGE consequence_types consequence_types SET (
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
                                            'protein_altering_variant'
                                        );

# insert into tmp table, replacing old terms with new as we go
INSERT INTO tmp_tv SELECT transcript_variation_id, variation_feature_id, feature_stable_id, allele_string, somatic,
    REPLACE(consequence_types, 'initiator_codon_variant', 'start_lost'),
cds_start, cds_end, cdna_start, cdna_end, translation_start, translation_end, distance_to_transcript, codon_allele_string, pep_allele_string, hgvs_genomic, hgvs_transcript, hgvs_protein, polyphen_prediction, polyphen_score, sift_prediction, sift_score, display
FROM transcript_variation;

DROP TABLE transcript_variation;
RENAME TABLE tmp_tv TO transcript_variation;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_79_80_c.sql|change the column consequence_types in transcript_variation and variation_feature to add protein_altering_variant and change initiator_codon_variant to start_lost');


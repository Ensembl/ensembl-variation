
# Update consequence_types in motif_feature_variation and regulatory_feature_variation tables
# Only store SO_terms that can be assigned to a variant overlapping a regulatory region 


ALTER TABLE motif_feature_variation MODIFY consequence_types SET (
        'TF_binding_site_variant',
        'TFBS_ablation',
        'TFBS_fusion',
        'TFBS_amplification',
        'TFBS_translocation'
    ) DEFAULT NULL;

ALTER TABLE regulatory_feature_variation MODIFY consequence_types SET (
                                            'regulatory_region_variant',
                                            'regulatory_region_ablation',
                                            'regulatory_region_fusion',
                                            'regulatory_region_amplification',
                                            'regulatory_region_translocation'
                                        ) DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_79_80_d.sql|Reduce consequence_terms to the set of relevant SO_terms in motif_feature_variation and regulatory_feature_variation tables');

# Add index on the source_id column in variation_feature and phenotype_feature
ALTER TABLE variation_feature ADD key source_idx (source_id);
ALTER TABLE phenotype_feature ADD key source_idx (source_id);

#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_f.sql|Add index on the source_id column in variation_feature and phenotype_feature');

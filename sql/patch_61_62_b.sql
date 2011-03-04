# drop subsnp_id column on failed_variation
ALTER TABLE failed_variation DROP COLUMN subsnp_id;

# Drop the key on variation_id column and make it a composite key with failed_description_id instead
ALTER TABLE failed_variation DROP KEY variation_idx;
ALTER TABLE failed_variation ADD UNIQUE KEY variation_idx (variation_id,failed_description_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_b.sql|change structure of failed_variation table');

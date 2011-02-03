# drop subsnp_id column on failed_variation
ALTER TABLE failed_variation DROP COLUMN subsnp_id;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_b.sql|change structure of failed_variation table');

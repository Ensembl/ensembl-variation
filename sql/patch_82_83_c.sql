# Drop the column validation_status in variation and variation_feature
ALTER TABLE `variation` DROP COLUMN `validation_status`;
ALTER TABLE `variation_feature` DROP COLUMN `validation_status`;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_82_83_c.sql|Drop the column validation_status in variation and variation_feature');


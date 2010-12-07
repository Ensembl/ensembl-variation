# add new failed description entry
INSERT IGNORE INTO failed_description (failed_description_id,description) VALUES (6,'Variation has no genotypes');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_e.sql|add new failed_description entry');

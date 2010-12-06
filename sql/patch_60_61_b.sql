# change structure of failed_variation table to allow
# multiple entries per variation and subsnp_id entry
ALTER TABLE failed_variation DROP PRIMARY KEY;
ALTER TABLE failed_variation ADD failed_variation_id int(11) NULL auto_increment PRIMARY KEY FIRST;
ALTER TABLE failed_variation ADD subsnp_id int(10) UNSIGNED NULL DEFAULT NULL  AFTER variation_id;
ALTER TABLE failed_variation ADD INDEX variation_idx (variation_id);
ALTER TABLE failed_variation ADD INDEX subsnp_idx (subsnp_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_b.sql|change structure of failed_variation table');

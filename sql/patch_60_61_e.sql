# change id to auto_increment
ALTER TABLE `failed_description` CHANGE `failed_description_id` `failed_description_id` int(10) UNSIGNED NOT NULL  auto_increment;

# add new failed description entry
INSERT IGNORE INTO failed_description (description) VALUES ('Variation has no genotypes');
INSERT IGNORE INTO failed_description (description) VALUES ('Genotype frequencies do not add up to 1');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_e.sql|add new failed_description entry');

# add url column to source table
ALTER TABLE `source` ADD `url` varchar(255) NULL DEFAULT NULL  AFTER `description`;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_57_58_d.sql|new column url in source table');

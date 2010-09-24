# Add flipped column to variation table - indicates if a variation has been strand flipped during import
ALTER TABLE `variation` ADD `flipped`  tinyint(1) unsigned NULL DEFAULT NULL AFTER `ancestral_allele`;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_59_60_c.sql|add flipped column to variation');

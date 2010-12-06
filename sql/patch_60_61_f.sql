# add allele_string column to structural_variation
ALTER TABLE `structural_variation` ADD `allele_string` longtext NULL DEFAULT NULL  AFTER `bound_end`;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_f.sql|add allele_string column to structural variation');

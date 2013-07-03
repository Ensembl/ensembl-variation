## add phased column to genotype_code table
ALTER TABLE `genotype_code` ADD `phased` TINYINT(2)  UNSIGNED  NULL  DEFAULT NULL  AFTER `haplotype_id`;

## patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_72_73_c.sql|Add phased column to genotype_code table');

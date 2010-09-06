# Add count column to allele and population_genotype tables
ALTER TABLE `allele` ADD `count` int(10) UNSIGNED NULL DEFAULT NULL  AFTER `sample_id`;
ALTER TABLE `population_genotype` ADD `count` int(10) UNSIGNED NULL DEFAULT NULL  AFTER `sample_id`;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_59_60_b.sql|add count column to allele and population_genotype');

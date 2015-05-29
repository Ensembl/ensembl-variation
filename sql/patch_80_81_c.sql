# Rename individual to sample in relevant table names
# Update column names and indices from individual to sample

RENAME TABLE individual_population TO sample_population;
ALTER TABLE sample_population CHANGE individual_id sample_id int(10) unsigned NOT NULL;
DROP INDEX `individual_idx` ON sample_population;
CREATE INDEX sample_idx ON sample_population(sample_id);

RENAME TABLE individual_genotype_multiple_bp TO sample_genotype_multiple_bp;
ALTER TABLE sample_genotype_multiple_bp CHANGE individual_id sample_id int(10) unsigned DEFAULT NULL;
DROP INDEX `individual_idx` ON sample_genotype_multiple_bp;
CREATE INDEX sample_idx ON sample_genotype_multiple_bp(sample_id);

ALTER TABLE compressed_genotype_region CHANGE individual_id sample_id int(10) unsigned NOT NULL;
DROP INDEX `individual_idx` ON compressed_genotype_region;
CREATE INDEX sample_idx ON compressed_genotype_region(sample_id);

ALTER TABLE structural_variation_sample CHANGE individual_id sample_id int(10) unsigned DEFAULT NULL;
DROP INDEX `individual_idx` ON structural_variation_sample;
CREATE INDEX sample_idx ON structural_variation_sample (sample_id);

ALTER TABLE read_coverage CHANGE individual_id sample_id int(10) unsigned NOT NULL;
CREATE INDEX sample_idx ON read_coverage (sample_id);

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_80_81_c.sql|Update table, column and index names from individual to sample.');


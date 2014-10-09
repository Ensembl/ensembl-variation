# Add a column copy_number for CNV supporting structural variants

ALTER TABLE structural_variation ADD COLUMN copy_number TINYINT(2) DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_77_78_b.sql|Add a column copy_number for CNV supporting structural variants');

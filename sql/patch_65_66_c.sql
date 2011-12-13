# populate the column is_evidence to structural variation feature table
ALTER TABLE structural_variation_feature ADD COLUMN is_evidence tinyint(1) DEFAULT 0 NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_65_66_c.sql|populate the column is_evidence to structural variation feature table');

# Drop the column strain_id from structural_variation_sample
ALTER TABLE structural_variation_sample DROP COLUMN strain_id;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_80_81_e.sql|Drop the column strain_id from structural_variation_sample');


# Drop the column strain_id from structural_variation_sample
UPDATE meta set meta_key='sample.reference_strain' WHERE meta_key='individual.reference_strain';
UPDATE meta set meta_key='sample.default_strain' WHERE meta_key='individual.default_strain';
UPDATE meta set meta_key='sample.display_strain' WHERE meta_key='individual.display_strain';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_80_81_f.sql|Update meta. Rename sample to individual.');


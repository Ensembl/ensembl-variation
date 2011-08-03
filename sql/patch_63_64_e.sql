# remove the column somatic from the structural variation table

ALTER TABLE structural_variation DROP COLUMN somatic;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_63_64_e.sql|remove the column somatic from the structural variation table');

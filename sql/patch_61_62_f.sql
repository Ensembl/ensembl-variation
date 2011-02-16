# add a study_id column in the tables structural_variation and variation_annotation

ALTER TABLE structural_variation ADD COLUMN study_id int(10);
ALTER TABLE structural_variation ADD INDEX study_idx (study_id);

ALTER TABLE variation_annotation ADD COLUMN local_study_id int(10);
ALTER TABLE variation_annotation ADD INDEX local_study_idx (local_study_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_f.sql|add a study_id column in the tables structural_variation and variation_annotation');

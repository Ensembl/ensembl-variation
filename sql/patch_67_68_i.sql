# change the type of the column study_type in the study table
ALTER TABLE study CHANGE study_type study_type varchar(255) DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_67_68_i.sql|change the type of the column study_type in the study table');

# change the default value of the study name column

ALTER TABLE study CHANGE name name varchar(255) DEFAULT NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_62_63_e.sql|change the default value of the study name column');

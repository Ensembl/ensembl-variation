# Change the type of the column description in the table study.

ALTER TABLE study CHANGE description description TEXT DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_71_72_i.sql|Change the type of the column description in the table study.');

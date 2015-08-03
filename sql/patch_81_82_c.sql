# add new description to the failed_description table

INSERT INTO failed_description (failed_description_id,description) VALUES (20,'Variant at first base in sequence');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_81_82_c.sql|new entry in the failed_description table');


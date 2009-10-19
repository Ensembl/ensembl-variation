# update the schema_version entry in the meta table
##################
update meta set meta_value = 56 where meta_key = 'schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_55_56_a.sql|schema version');

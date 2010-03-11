# update the schema_version entry in the meta table
##################
update meta set meta_value = 58 where meta_key = 'schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_57_58_a.sql|schema version');

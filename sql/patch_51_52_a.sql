# update the schema_version entry in the meta table
##################
update meta set meta_value = 52 where meta_key = 'schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_51_52_a.sql|schema version');

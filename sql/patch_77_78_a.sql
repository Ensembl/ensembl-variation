

# update the schema_version entry in the meta table
##################
UPDATE meta SET meta_value = 78 WHERE meta_key = 'schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_77_78_a.sql|schema version');

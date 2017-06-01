# update the schema_version entry in the meta table
UPDATE meta SET meta_value = 90 WHERE meta_key = 'schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_89_90_a.sql|schema version');


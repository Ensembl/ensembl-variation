# update the schema_version entry in the meta table
##################
update meta set meta_value = 48 where meta_key = 'schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_47_48_a.sql|schema version');

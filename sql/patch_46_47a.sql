# update the schema_version entry in the meta table
##################
update meta set meta_value = 47 where meta_key = 'schema_version';

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_46_47a.sql|schema version');

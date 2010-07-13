# Add schema_type entry to meta table
INSERT INTO meta (species_id,meta_key, meta_value) VALUES (NULL,'schema_type', 'variation');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_58_59_c.sql|schema_type entry in meta table');

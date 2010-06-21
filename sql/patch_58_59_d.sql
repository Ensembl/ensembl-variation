# Add schema_type entry to meta table
ALTER TABLE source ADD COLUMN `somatic` tinyint(1) DEFAULT '0';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_58_59_d.sql|add somatic column to source table');

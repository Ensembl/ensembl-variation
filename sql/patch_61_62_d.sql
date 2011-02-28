# add a type column in the source table

ALTER TABLE source ADD COLUMN type ENUM ('chip');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_d.sql|add a type column in the source table');

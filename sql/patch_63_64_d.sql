# add 'lsdb' to the type enum in source table

ALTER TABLE source CHANGE COLUMN type type ENUM('chip','lsdb') DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_63_64_d.sql|add lsdb entry to enum in type column');

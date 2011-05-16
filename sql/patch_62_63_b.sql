# add a column for specifying a variation set short name as an attribute
##################
ALTER TABLE variation_set ADD COLUMN short_name_attrib_id INT(10) UNSIGNED DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_62_63_b.sql|add a column for specifying a variation set short name as an attribute');

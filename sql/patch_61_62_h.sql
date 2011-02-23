# use the attrib_id for the var_class

DROP TABLE IF EXISTS variation_class;

ALTER TABLE variation CHANGE class_so_id class_attrib_id int(10) unsigned default 0;
ALTER TABLE variation_feature CHANGE class_so_id class_attrib_id int(10) unsigned default 0;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_h.sql|use the attrib_id for the var_class');

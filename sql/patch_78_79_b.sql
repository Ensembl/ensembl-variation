
# change the column attrib_type_id by attrib_id in the variation_attrib table
##################
ALTER TABLE variation_attrib CHANGE attrib_type_id attrib_id INT(11) DEFAULT NULL;
ALTER TABLE variation_attrib DROP KEY type_value_idx;
ALTER TABLE variation_attrib ADD KEY attrib_value_idx (attrib_id,value);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_78_79_b.sql|change the column attrib_type_id by attrib_id in the variation_attrib table');

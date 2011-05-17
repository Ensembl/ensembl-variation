# add a unique constraint on the index for attrib_type_id and value columns in attrib table 
##################
ALTER TABLE attrib DROP INDEX type_val_idx;
ALTER TABLE attrib ADD UNIQUE INDEX type_val_idx (attrib_type_id, value(40));

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_62_63_c.sql|add a unique constraint on the index for attrib_type_id and value columns in attrib table');

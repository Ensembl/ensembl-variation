# Extend the index type_val_idx in the table attrib
ALTER TABLE attrib DROP KEY type_val_idx;
ALTER TABLE attrib ADD CONSTRAINT UNIQUE KEY type_val_idx (attrib_type_id, value(80));

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_77_78_d.sql|Extend the index type_val_idx in the table attrib');

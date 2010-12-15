# change the p_value datatype from varchar to double in the variation_annotation table

ALTER TABLE variation_annotation CHANGE p_value p_value DOUBLE;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_k.sql|change the p_value datatype from varchar to double in the variation_annotation table');

# add indexes in the table structural_variation_association
##################

ALTER TABLE structural_variation_association ADD INDEX structural_variation_idx (structural_variation_id);
ALTER TABLE structural_variation_association ADD INDEX supporting_structural_variation_idx (supporting_structural_variation_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_66_67_c.sql|add indexes in the table structural_variation_association');

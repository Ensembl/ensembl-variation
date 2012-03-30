# update the structural_variation schema
##################

ALTER TABLE structural_variation_association ADD INDEX structural_variation_idx (structural_variation_id);
ALTER TABLE structural_variation_association ADD INDEX supporting_structural_variation_idx (supporting_structural_variation_id);

ALTER TABLE structural_variation ADD COLUMN somatic TINYINT(1) NOT NULL DEFAULT 0;
ALTER TABLE structural_variation_feature ADD COLUMN somatic TINYINT(1) NOT NULL DEFAULT 0;
ALTER TABLE structural_variation_feature ADD COLUMN breakpoint_order TINYINT(4) DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_66_67_c.sql|update the structural_variation schema');

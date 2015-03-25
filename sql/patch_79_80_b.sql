

# create a unique key for the variation_name column in the table structural_variation
##################
ALTER TABLE structural_variation DROP KEY name_idx;
ALTER TABLE structural_variation ADD CONSTRAINT UNIQUE KEY (`variation_name`);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_79_80_b.sql|create a unique key for the variation_name column in the table structural_variation');

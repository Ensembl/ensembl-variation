## add pos_idx key to phenotype_feature
ALTER TABLE phenotype_feature ADD KEY `pos_idx` (`seq_region_id`,`seq_region_start`,`seq_region_end`)

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_j.sql|add pos_idx key to phenotype_feature');

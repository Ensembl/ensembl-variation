## Update indexes for the phenotype table
ALTER TABLE phenotype DROP KEY name_idx;
ALTER TABLE phenotype ADD KEY name_idx (name);
ALTER TABLE phenotype ADD CONSTRAINT desc_idx UNIQUE (description);

#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_e.sql|Update indexes for the phenotype table');
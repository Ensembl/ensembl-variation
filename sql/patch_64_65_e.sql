# add meta entry for ploidy
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'ploidy', 2);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_64_65_e.sql|add meta entry for ploidy');

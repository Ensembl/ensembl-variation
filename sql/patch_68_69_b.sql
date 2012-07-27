# update the variation_feature table
##################


ALTER TABLE variation_feature ADD COLUMN alignment_quality double  DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_68_69_b.sql|update the variation_feature table');

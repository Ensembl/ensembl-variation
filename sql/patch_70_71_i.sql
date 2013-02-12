## drop old annotation tables
DROP TABLE variation_annotation;
DROP TABLE structural_variation_annotation;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_i.sql|drop old annotation tables');

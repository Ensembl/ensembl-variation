# remove unused tables
DROP TABLE variation_group;
DROP TABLE variation_group_variation;
DROP TABLE variation_group_feature;
DROP TABLE allele_group;
DROP TABLE allele_group_allele;
DROP TABLE httag;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_65_66_b.sql|remove unused tables');


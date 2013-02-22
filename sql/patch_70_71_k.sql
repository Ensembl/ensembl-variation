## add pos_idx key to phenotype_feature
ALTER TABLE `phenotype_feature_attrib` ADD INDEX type_value_idx (`attrib_type_id`, `value`);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_k.sql|add type_value_idx key to phenotype_feature_attrib');

ALTER TABLE `variation` CHANGE `evidence_attribs` `evidence_attribs` SET('367','368','369','370','371','372','418');

ALTER TABLE `variation_feature` CHANGE `evidence_attribs` `evidence_attribs` SET('367','368','369','370','371','372','418');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_79_80_f.sql|add Phenotype or Disease evidence_attribs');

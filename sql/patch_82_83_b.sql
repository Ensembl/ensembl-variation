# Add the evidence "ExAC" in variation and variation_feature
ALTER TABLE `variation` CHANGE `evidence_attribs` `evidence_attribs` SET('367','368','369','370','371','372','418','421') DEFAULT NULL;

ALTER TABLE `variation_feature` CHANGE `evidence_attribs` `evidence_attribs` SET('367','368','369','370','371','372','418','421') DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_82_83_b.sql|Add the evidence ExAC in variation and variation_feature');


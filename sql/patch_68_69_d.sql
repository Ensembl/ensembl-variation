# additional variant fail class- fail if >1 match to reference rather than >3 as previous
##################

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_68_69_d.sql|additional variant fail class');


INSERT INTO failed_description (failed_description_id,description) VALUES (19,'Variation maps to more than one genomic location');

UPDATE failed_variation SET failed_description_id = 19 WHERE failed_description_id=1;

INSERT INTO failed_variation (variation_id, failed_description_id)  ( SELECT DISTINCT variation_id, '19' FROM variation_feature WHERE map_weight >1 AND map_weight < 3);

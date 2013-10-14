## Add clinical_significance to variation_feature (copied from variation)
ALTER TABLE variation_feature
ADD COLUMN clinical_significance
SET('drug-response','histocompatibility','non-pathogenic','other','pathogenic','probable-non-pathogenic','probable-pathogenic''unknown','untested') DEFAULT NULL;

UPDATE variation_feature vf, variation v
SET vf.clinical_significance = v.clinical_significance
WHERE v.variation_id = vf.variation_id;

#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_c.sql|Add clinical_significance to variation_feature table');

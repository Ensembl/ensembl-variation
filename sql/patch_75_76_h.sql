# Replace the column clinical_significance_attrib_id by clinical_significance in structural_variation

ALTER TABLE structural_variation ADD COLUMN clinical_significance SET('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective') DEFAULT NULL AFTER clinical_significance_attrib_id;

ALTER TABLE structural_variation ADD COLUMN tmp_clin_sign VARCHAR(255) DEFAULT NULL;

UPDATE structural_variation sv, attrib a SET sv.tmp_clin_sign=LOWER(a.value) WHERE sv.clinical_significance_attrib_id=a.attrib_id AND sv.clinical_significance_attrib_id is NOT NULL;

UPDATE structural_variation SET clinical_significance=tmp_clin_sign WHERE tmp_clin_sign is NOT NULL;

ALTER TABLE structural_variation DROP clinical_significance_attrib_id, DROP tmp_clin_sign;

#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_h.sql|Replace the column clinical_significance_attrib_id by clinical_significance in structural_variation');


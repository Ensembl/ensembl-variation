
## update clinical significance terms to match ClinVar website rather than dbSNP


ALTER TABLE variation MODIFY column clinical_significance set ('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective');

ALTER TABLE variation_feature MODIFY column clinical_significance set ('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective');


#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_g.sql|update variation and variation_feature to use the same clinical significance terms as ClinVar');

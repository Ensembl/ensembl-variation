-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.


# Create a table structural_variation_sample. Move the structural variation data from phenotype_feature_attrib to structural_variation_sample and structural_variation.

# Create the table structural_variation_sample

CREATE TABLE structural_variation_sample (
	structural_variation_sample_id int(10) unsigned NOT NULL auto_increment,
	structural_variation_id int(10) unsigned NOT NULL,
	individual_id int(10) unsigned DEFAULT NULL,
	strain_id int(10) unsigned DEFAULT NULL,
	
	primary key (structural_variation_sample_id),
	key structural_variation_idx(structural_variation_id),
	key individual_idx(individual_id),
	key strain_idx(strain_id)
);


# Fill the table structural_variation_sample

INSERT IGNORE INTO structural_variation_sample (structural_variation_id,individual_id,strain_id)
SELECT DISTINCT sv.structural_variation_id,pfa1.value,pfa2.value
FROM structural_variation sv,phenotype_feature pf
	LEFT JOIN phenotype_feature_attrib pfa1 ON 
		(pfa1.phenotype_feature_id=pf.phenotype_feature_id AND pfa1.attrib_type_id IN 
		  (SELECT at1.attrib_type_id FROM attrib_type at1 WHERE at1.code='sample_id')
		)
	LEFT JOIN phenotype_feature_attrib pfa2 ON 
		(pfa2.phenotype_feature_id=pf.phenotype_feature_id AND pfa2.attrib_type_id IN 
		  (SELECT at2.attrib_type_id FROM attrib_type at2 WHERE at2.code='strain_id')
		)
WHERE sv.variation_name=pf.object_id AND (pf.type='StructuralVariation' OR pf.type='SupportingStructuralVariation');


# Move the clinical significance to the table structural_variation

ALTER TABLE structural_variation ADD COLUMN tmp_clin_sign VARCHAR(255) DEFAULT NULL;

UPDATE structural_variation sv, phenotype_feature pf, phenotype_feature_attrib pfa, attrib_type at SET sv.tmp_clin_sign=pfa.value 
WHERE 
  sv.variation_name=pf.object_id AND
	pf.phenotype_feature_id=pfa.phenotype_feature_id AND
	at.attrib_type_id=pfa.attrib_type_id AND
	at.code='dgva_clin_sig' AND
	(pf.type='StructuralVariation' OR pf.type='SupportingStructuralVariation')
;
UPDATE structural_variation sv, attrib a, attrib_type at SET sv.clinical_significance_attrib_id=a.attrib_id 
WHERE a.value=sv.tmp_clin_sign AND a.attrib_type_id=at.attrib_type_id AND at.code='dgva_clin_sig' AND sv.tmp_clin_sign is not null;

ALTER TABLE structural_variation DROP COLUMN tmp_clin_sign;


# Delete entries in phenotype_feature without phenotype_id

DELETE FROM phenotype_feature_attrib WHERE phenotype_feature_id IN 
(SELECT phenotype_feature_id FROM phenotype_feature WHERE type='StructuralVariation' OR type='SupportingStructuralVariation');

DELETE FROM phenotype_feature WHERE phenotype_id is null AND (type='StructuralVariation' OR type='SupportingStructuralVariation');


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_71_72_h.sql|Create a table structural_variation_sample. Move the structural variation data from phenotype_feature_attrib to structural_variation_sample and structural_variation.');

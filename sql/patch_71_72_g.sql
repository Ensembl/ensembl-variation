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


# Add new columns to the tables structural_variation and structural_variation_feature.

ALTER TABLE structural_variation ADD COLUMN alias VARCHAR(255) DEFAULT NULL AFTER variation_name;
ALTER TABLE structural_variation ADD COLUMN clinical_significance_attrib_id int(10) unsigned DEFAULT NULL AFTER class_attrib_id;
ALTER TABLE structural_variation ADD KEY clinical_attrib_idx (clinical_significance_attrib_id);

ALTER TABLE structural_variation_feature ADD COLUMN length INT(10) DEFAULT NULL;
ALTER TABLE structural_variation_feature ADD COLUMN study_id INT(10) unsigned DEFAULT NULL AFTER source_id; 

UPDATE structural_variation_feature svf, structural_variation sv SET svf.study_id=sv.study_id 
WHERE svf.structural_variation_id=sv.structural_variation_id and sv.study_id is not NULL;

ALTER TABLE structural_variation_feature ADD KEY study_idx (study_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_71_72_g.sql|Add new columns to the tables structural_variation and structural_variation_feature.');

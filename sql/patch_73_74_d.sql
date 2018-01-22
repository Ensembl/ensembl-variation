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


## Add data_types to source table
CREATE TABLE tmp_source_type (
  source_id INT,
  source_type VARCHAR(255)
);

INSERT INTO tmp_source_type SELECT s.source_id, 'variation' FROM source s, variation v WHERE s.source_id = v.source_id GROUP BY s.source_id;
INSERT INTO tmp_source_type SELECT s.source_id, 'variation_synonym' FROM source s, variation_synonym v WHERE s.source_id = v.source_id GROUP BY s.source_id;
INSERT INTO tmp_source_type select s.source_id, 'structural_variation' FROM source s, structural_variation v WHERE s.source_id = v.source_id GROUP BY s.source_id;
INSERT INTO tmp_source_type SELECT s.source_id, 'phenotype_feature' FROM source s, phenotype_feature v WHERE s.source_id = v.source_id GROUP BY s.source_id;
INSERT INTO tmp_source_type SELECT s.source_id, 'study' FROM source s, study v WHERE s.source_id = v.source_id GROUP BY s.source_id;

CREATE TABLE tmp_source_type_grouped SELECT source_id, group_concat(source_type) AS data_types FROM tmp_source_type GROUP BY source_id;
ALTER TABLE source ADD data_types SET('variation','variation_synonym','structural_variation','phenotype_feature','study')  NULL  DEFAULT NULL;
UPDATE source s, tmp_source_type_grouped t SET s.data_types = t.data_types WHERE s.source_id = t.source_id;

DROP TABLE tmp_source_type;
DROP TABLE tmp_source_type_grouped;

#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_d.sql|Add data_types to source table');

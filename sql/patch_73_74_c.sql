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


## Add clinical_significance to variation_feature (copied from variation)
ALTER TABLE variation 
MODIFY COLUMN clinical_significance 
SET('drug-response','histocompatibility','non-pathogenic','other','pathogenic','probable-non-pathogenic','probable-pathogenic','unknown','untested');

ALTER TABLE variation_feature
ADD COLUMN clinical_significance
SET('drug-response','histocompatibility','non-pathogenic','other','pathogenic','probable-non-pathogenic','probable-pathogenic','unknown','untested') DEFAULT NULL;

UPDATE variation_feature vf, variation v
SET vf.clinical_significance = v.clinical_significance
WHERE v.variation_id = vf.variation_id;

#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_c.sql|Add clinical_significance to variation_feature table');

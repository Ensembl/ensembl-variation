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


# update the structural_variation schema
##################

ALTER TABLE structural_variation_association ADD INDEX structural_variation_idx (structural_variation_id);
ALTER TABLE structural_variation_association ADD INDEX supporting_structural_variation_idx (supporting_structural_variation_id);

ALTER TABLE structural_variation ADD COLUMN somatic TINYINT(1) NOT NULL DEFAULT 0;
ALTER TABLE structural_variation_feature ADD COLUMN somatic TINYINT(1) NOT NULL DEFAULT 0;
ALTER TABLE structural_variation_feature ADD COLUMN breakpoint_order TINYINT(4) DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_66_67_c.sql|update the structural_variation schema');

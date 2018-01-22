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


# add a study_id column in the tables structural_variation and variation_annotation

ALTER TABLE structural_variation ADD COLUMN study_id int(10);
ALTER TABLE structural_variation ADD INDEX study_idx (study_id);

ALTER TABLE variation_annotation ADD COLUMN study_id int(10);
ALTER TABLE variation_annotation ADD INDEX study_idx (study_id);

ALTER TABLE variation_annotation DROP COLUMN local_stable_id;
ALTER TABLE variation_annotation DROP COLUMN study;
ALTER TABLE variation_annotation DROP COLUMN study_type;
ALTER TABLE variation_annotation DROP COLUMN source_id;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_f.sql|add a study_id column in the tables structural_variation and variation_annotation');

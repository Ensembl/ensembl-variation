-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2025] EMBL-European Bioinformatics Institute
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

#  The UNIQUE KEY `name` (`name`,`source_id`) 
#  only allows one synonym for a source
#  Modifying to allow only one synonym for a source for a variation
ALTER TABLE variation_synonym DROP INDEX name, ADD UNIQUE INDEX name_idx (name, source_id, variation_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_95_96_b.sql|modify index on variation_synonym');

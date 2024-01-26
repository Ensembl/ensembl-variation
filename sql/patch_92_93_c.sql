-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

# add new evidence attribs the variation and variation_feature tables

ALTER TABLE `variation` CHANGE `evidence_attribs` `evidence_attribs` SET('367','368','369','370','371','372','418','421','573','585') DEFAULT NULL;

ALTER TABLE `variation_feature` CHANGE `evidence_attribs` `evidence_attribs` SET('367','368','369','370','371','372','418','421','573','585') DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_c.sql|Add new evidence attribs to the variation and variation_feature tables');

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


# change id to auto_increment
ALTER TABLE `failed_description` CHANGE `failed_description_id` `failed_description_id` int(10) UNSIGNED NOT NULL  auto_increment;

# add new failed description entry
INSERT IGNORE INTO failed_description (description) VALUES ('Variation has no genotypes');
INSERT IGNORE INTO failed_description (description) VALUES ('Genotype frequencies do not add up to 1');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_e.sql|add new failed_description entry');

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


# change structure of failed_variation table to allow
# multiple entries per variation and subsnp_id entry
ALTER TABLE failed_variation DROP PRIMARY KEY;
ALTER TABLE failed_variation ADD failed_variation_id int(11) NULL auto_increment PRIMARY KEY FIRST;
ALTER TABLE failed_variation ADD subsnp_id int(10) UNSIGNED NULL DEFAULT NULL  AFTER variation_id;
ALTER TABLE failed_variation ADD INDEX variation_idx (variation_id);
ALTER TABLE failed_variation ADD INDEX subsnp_idx (subsnp_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_b.sql|change structure of failed_variation table');

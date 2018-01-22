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


# populate the columns is_evidence and variation_set_id into the structural variation feature table
ALTER TABLE structural_variation_feature ADD COLUMN is_evidence tinyint(1) DEFAULT 0 NOT NULL;


# Add a column for variation_set_ids to strucutral_variation_feature
ALTER TABLE `structural_variation_feature` ADD `variation_set_id` SET 
('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64') NOT NULL DEFAULT '';
# Add an index on the variation_set_id column
ALTER TABLE `structural_variation_feature` ADD INDEX `variation_set_idx` (`variation_set_id`);


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_65_66_c.sql|populate the columns is_evidence and variation_set_id into the structural variation feature table');

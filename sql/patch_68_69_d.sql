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


# additional variant fail class- fail if >1 match to reference rather than >3 as previous
##################

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_68_69_d.sql|additional variant fail class');


INSERT INTO failed_description (failed_description_id,description) VALUES (19,'Variation maps to more than one genomic location');

UPDATE failed_variation SET failed_description_id = 19 WHERE failed_description_id=1;

INSERT INTO failed_variation (variation_id, failed_description_id)  ( SELECT DISTINCT(variation_id), '19' FROM variation_feature WHERE map_weight >1 AND map_weight < 3);

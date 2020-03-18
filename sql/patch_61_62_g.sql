-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2020] EMBL-European Bioinformatics Institute
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


# move the somatic flag from the source to variation and variation_feature tables

ALTER TABLE source DROP COLUMN somatic;
ALTER TABLE variation ADD COLUMN somatic tinyint(1) DEFAULT 0 NOT NULL;
ALTER TABLE variation_feature ADD COLUMN somatic tinyint(1) DEFAULT 0 NOT NULL;

# patch identifier

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_g.sql|move somatic flag from source to variation and variation_feature tables');

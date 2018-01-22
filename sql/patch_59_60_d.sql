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


# drop and rebuild index on variation_feature and structural_variation
alter table variation_feature drop index pos_idx;
alter table structural_variation drop index pos_idx;

create index pos_idx on variation_feature (seq_region_id, seq_region_start, seq_region_end);
create index pos_idx on structural_variation (seq_region_id, seq_region_start, seq_region_end);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_59_60_d.sql|modify index on variation_feature and structural_variation');

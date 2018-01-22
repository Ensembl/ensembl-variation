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


###################################################
# add a new foreign key in variation and variaiton_synonym table, to allow get variations by source
###################################################
ALTER TABLE variation ADD KEY source_idx (source_id);

ALTER TABLE variation_synonym ADD KEY source_idx (source_id);

INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_47_48_b.sql|source_id key');

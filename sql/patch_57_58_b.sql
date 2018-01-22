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


# update add two entries cds_start, cds_end in the transcript_variation table
##################
alter table transcript_variation add column cds_start int(11) after cdna_end;
alter table transcript_variation add column cds_end int(11) after cds_start;
ALTER TABLE transcript_variation CHANGE transcript_id transcript_stable_id VARCHAR(128) NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_57_58_b.sql|add cds_start, cds_end, change transcript_id to stable_id in transcript_variation');

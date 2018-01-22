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


# replace variation.clinical_significance_attrib_id with variation.clinical_significance set to allow different statuses for different traits

alter table variation add column clinical_significance SET('drug-response','histocompatibility','non-pathogenic','other','pathogenic','probable-non-pathogenic','probable-pathogenic''unknown','untested');


alter table variation drop column clinical_significance_attrib_id;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_71_72_e.sql|change variation clinical_significance column');

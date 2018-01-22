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


# add some new columns to the variation table to store new dbSNP 134 stuff

ALTER TABLE variation 
    ADD COLUMN minor_allele char(1) DEFAULT NULL,
    ADD COLUMN minor_allele_freq float DEFAULT NULL,
    ADD COLUMN minor_allele_count int(10) unsigned DEFAULT NULL, 
    ADD COLUMN clinical_significance_attrib_id int(10) unsigned DEFAULT NULL
;

# and add a new failed_description

INSERT INTO failed_description (failed_description_id,description) VALUES (16,'Flagged as suspect by dbSNP');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_64_65_i.sql|add support for new data types from dbSNP');


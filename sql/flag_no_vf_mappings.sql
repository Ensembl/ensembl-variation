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


#sql patch to find out variations that don't have variation_feature mappings
#first, if there is no table to store the failed variations, create it
CREATE TABLE IF NOT EXISTS failed_variation(
    variation_id int(10) unsigned not null,
    failed_description_id int(10) unsigned not null,

    PRIMARY KEY(variation_id)
);

#then, copy the information from the variation_feature table
INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) 
   SELECT v.variation_id, 5
   FROM  variation v LEFT JOIN variation_feature vf
   ON v.variation_id=vf.variation_id
   WHERE vf.variation_id is null;

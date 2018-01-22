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


#sql patch to find out variations that contain more than 3 alleles
#first, if there is no table to store the failed variations, create it
CREATE TABLE IF NOT EXISTS failed_variation(
    variation_id int(10) unsigned not null,
    failed_description_id int(10) unsigned not null,

    PRIMARY KEY(variation_id)
);

#then, copy the information from the variation_feature table
INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) 
   SELECT variation_id, 3
   FROM  variation_feature
   WHERE length(allele_string) - length(REPLACE(allele_string,'/','')) > 2
   AND allele_string not like '%-%';

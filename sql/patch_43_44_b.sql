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


####################
# change the schema of the database to include the failed_variation table
####################

CREATE TABLE IF NOT EXISTS failed_variation(
    variation_id int(10) unsigned not null,
    failed_description_id int(10) unsigned not null,

    PRIMARY KEY(variation_id)
);

#import the data in the failed_variation table

INSERT INTO failed_variation (variation_id, failed_description_id)
  SELECT variation_id, failed_description_id
  FROM variation
  WHERE failed_description_id <> 0;

#and remove the column from the variation table

ALTER TABLE variation DROP failed_description_id;

#and update the meta table with the patch applied
INSERT INTO meta (meta_key,meta_value) VALUES ('patch','patch_43_44_b.sql|update database schema');	

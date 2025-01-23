-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2025] EMBL-European Bioinformatics Institute
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

# create table to store allele synonyms for a variant

CREATE TABLE allele_synonym (
  allele_synonym_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  variation_id      int(10) unsigned NOT NULL,
  hgvs_genomic      varchar(600) NOT NULL,
  name              varchar(255) NOT NULL,
  PRIMARY KEY (allele_synonym_id),
  UNIQUE KEY variation_name_idx (variation_id, name),
  KEY name_idx (name)
) ENGINE = MyISAM;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_94_95_b.sql|create table to store allele synonyms');

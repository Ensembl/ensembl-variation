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


# Replace sample_synonym by individual_synonym and population_synonym

CREATE TABLE individual_synonym (
  synonym_id int(10) unsigned not null auto_increment,
  individual_id int(10) unsigned not null,
  source_id int(10) unsigned not null,
  name varchar(255),

  primary key(synonym_id),
  key individual_idx (individual_id),
  key (name, source_id)
);

INSERT INTO individual_synonym(individual_id, source_id, name)
SELECT i.individual_id, s.source_id, s.name 
FROM individual i, sample_synonym s
WHERE i.individual_id = s.sample_id; 

CREATE TABLE population_synonym (
  synonym_id int(10) unsigned not null auto_increment,
  population_id int(10) unsigned not null,
  source_id int(10) unsigned not null,
  name varchar(255),

  primary key(synonym_id),
  key population_idx (population_id),
  key (name, source_id)
);

INSERT INTO population_synonym(population_id, source_id, name)
SELECT p.population_id, s.source_id, s.name
FROM population p, sample_synonym s
WHERE p.population_id = s.sample_id;

DROP TABLE sample_synonym;
DROP TABLE sample;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_71_72_d.sql|Replace sample_synonym by individual_synonym and population_synonym');

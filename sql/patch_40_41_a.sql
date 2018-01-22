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


#this patch is suposed to change the sample, individual and population schema

#delete from the Allele table the alleles that belong to strains
DELETE a FROM allele a, population p where p.sample_id = a.sample_id and p.is_strain = 1;
#delete from Individual_population table the relation between population->individual
DELETE ip FROM individual_population ip, population p where p.sample_id = ip.population_sample_id and p.is_strain = 1;
#and delete from the Sample table all strains
DELETE s FROM sample s, population p where s.sample_id = p.sample_id and p.is_strain = 1;
#and the strains from the population table
DELETE p FROM population p where p.is_strain = 1;
#remove the is_strain column from the population table
ALTER TABLE population drop column is_strain;
#create new individual_type table
create table individual_type(
  individual_type_id int(0) unsigned not null auto_increment,
  name varchar(255) not null,
  description text,
  
  primary key (individual_type_id)
);

#this table will always contain the same values

INSERT INTO individual_type (name,description) VALUES ('fully_inbred','multiple organisms have the same genome sequence');
INSERT INTO individual_type (name,description) VALUES ('partly_inbred','single organisms have reduced genome variability due to human intervention');
INSERT INTO individual_type (name,description) VALUES ('outbred','a single organism which breeds freely');
INSERT INTO individual_type (name,description) VALUES ('mutant','a single or multiple organisms with the same genome sequence that have a natural or experimentally induced mutation');

#modify table Individual with a new column
ALTER TABLE individual add column individual_type_id int(10) unsigned not null;
#we have to set the value of individual_type_id in the individual table depending on the specie
UPDATE individual set individual_type_id = 1 where (database() like 'mus%');
UPDATE individual set individual_type_id = 2 where (database() like 'canis%' or database() like 'danio%' or database() like 'gallus%' or database() like 'rattus%');
UPDATE individual set individual_type_id = 3 where (database() like 'homo%' or database() like 'anopheles%');

UPDATE meta set meta_key = 'individual.default_strain' where meta_key = 'population.default_strain';
UPDATE meta set meta_key = 'individual.reference_strain' where meta_key = 'population.reference_strain';
UPDATE meta set meta_key = 'individual.display_strain' where meta_key = 'population.display_strain';
INSERT INTO meta (meta_key,meta_value) VALUES ('patch','patch_40_41_a.sql|change population/individual schema');

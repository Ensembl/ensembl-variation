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


## Schema 31-32
## added a new concept in the database, sample
## that tries to merge into a base class the 
## individual and population concept. There is
## the addition of a new table, Sample, main reason
## to maintain the unique sample_id identifier

CREATE TABLE sample(
	sample_id int not null auto_increment,
	name varchar(255) not null,
	size int,
	description text,

	primary key( sample_id ),
	key name_idx( name )
);
## in order to keep the relation between old population_id and individual_id
## with sample_id, need to add the field in the tables

ALTER TABLE sample ADD (individual_id int), ADD (population_id int);

## need to move to Sample table Populations
INSERT INTO sample (name, size, description, population_id)
	SELECT name, size, description, population_id
	FROM population;

## and update all previous population_id field with the new sample_id
## tables needed the update: population_structure, population_synonym, allele,
## individual_population, allele_group, tagged_variation_feature, pairwise_ld,
## population_genotype

ALTER TABLE population ADD sample_id int not null;
UPDATE population p, sample s
SET p.sample_id = s.sample_id
WHERE p.population_id = s.population_id;
ALTER TABLE population DROP population_id, DROP name, DROP size, DROP description, ADD KEY sample_idx (sample_id);


ALTER TABLE population_structure ADD super_population_sample_id int not null;
UPDATE population_structure ps, sample s
SET ps.super_population_sample_id = s.sample_id
WHERE ps.super_population_id = s.population_id;


ALTER TABLE population_structure ADD sub_population_sample_id int not null;
UPDATE population_structure ps, sample s
SET ps.sub_population_sample_id = s.sample_id
WHERE ps.sub_population_id = s.population_id;
#and create the indexes
ALTER TABLE population_structure DROP super_population_id, DROP sub_population_id,
	ADD UNIQUE super_population_sample_id (super_population_sample_id, sub_population_sample_id), 
	ADD INDEX sub_pop_sample_idx (sub_population_sample_id, super_population_sample_id);

ALTER TABLE population_synonym RENAME sample_synonym; #change table name

ALTER TABLE sample_synonym CHANGE population_synonym_id sample_synonym_id int not null auto_increment; #change field name
ALTER TABLE sample_synonym ADD sample_id int not null;
UPDATE sample_synonym ss, sample s
SET ss.sample_id = s.sample_id
WHERE ss.population_id = s.population_id;
ALTER TABLE sample_synonym DROP population_id, ADD INDEX sample_idx (sample_id), DROP INDEX name, ADD KEY name (name,source_id);


ALTER TABLE allele ADD sample_id int;
UPDATE allele a, sample s
SET a.sample_id = s.sample_id
WHERE a.population_id = s.population_id;
ALTER TABLE allele DROP population_id;

ALTER TABLE individual_population ADD population_sample_id int not null;
UPDATE individual_population ip, sample s
SET ip.population_sample_id = s.sample_id
WHERE ip.population_id = s.population_id;
ALTER TABLE individual_population DROP population_id;

ALTER TABLE allele_group ADD sample_id int;
UPDATE allele_group ag, sample s
SET ag.sample_id = s.sample_id
WHERE ag.population_id = s.population_id;
ALTER TABLE allele_group DROP population_id;

#ALTER TABLE tagged_variation_feature ADD sample_id int not null;
#UPDATE tagged_variation_feature tg, sample s
#SET tg.sample_id = s.sample_id
#WHERE tg.population_id = s.population_id;
#ALTER TABLE tagged_variation_feature DROP PRIMARY KEY, DROP population_id, ADD PRIMARY KEY (variation_feature_id, sample_id);

ALTER TABLE pairwise_ld ADD sample_id int not null;
UPDATE pairwise_ld p, sample s
SET p.sample_id = s.sample_id
WHERE p.population_id = s.population_id;
ALTER TABLE pairwise_ld DROP population_id;

ALTER TABLE population_genotype ADD sample_id int not null;
UPDATE population_genotype p, sample s
SET p.sample_id = s.sample_id
WHERE p.population_id = s.population_id;
ALTER TABLE population_genotype DROP population_id, ADD INDEX sample_idx (sample_id);

## move now to Sample table Individual 
INSERT INTO sample (name, description, individual_id)
	SELECT name, description, individual_id
	FROM individual;

## and update the individual_id field in the tables: individual,
## individual_population, individual_genotypes tables

ALTER TABLE individual ADD sample_id int not null, ADD father_individual_sample_id int, ADD mother_individual_sample_id int;

UPDATE individual i, sample s
SET i.sample_id = s.sample_id
WHERE i.individual_id = s.individual_id;

#update the father_individual_id
UPDATE individual i, sample s
SET i.father_individual_sample_id = s.sample_id
WHERE i.father_individual_id = s.individual_id;

#update the mother_individual_id
UPDATE individual i, sample s
SET i.mother_individual_sample_id = s.sample_id
WHERE i.mother_individual_id = s.individual_id;

## and remove extra fields
ALTER TABLE individual DROP individual_id, DROP name, DROP description, ADD PRIMARY KEY (sample_id), 
	DROP father_individual_id, DROP mother_individual_id;


ALTER TABLE individual_population ADD individual_sample_id int not null;
UPDATE individual_population ip, sample s
SET ip.individual_sample_id = s.sample_id
WHERE ip.individual_id = s.individual_id;
ALTER TABLE individual_population DROP individual_id, ADD KEY individual_sample_idx (individual_sample_id), 
	ADD KEY population_sample_idx (population_sample_id);

ALTER TABLE individual_genotype_single_bp ADD sample_id int not null;
UPDATE individual_genotype_single_bp i, sample s
SET i.sample_id = s.sample_id
WHERE i.individual_id = s.individual_id;
ALTER TABLE individual_genotype_single_bp DROP individual_id, ADD KEY sample_idx (sample_id);

ALTER TABLE individual_genotype_multiple_bp ADD sample_id int not null;
UPDATE individual_genotype_multiple_bp i, sample s
SET i.sample_id = s.sample_id
WHERE i.individual_id = s.individual_id;
ALTER TABLE individual_genotype_multiple_bp DROP individual_id, ADD KEY sample_idx (sample_id);

##and update the sample_synonym table with the new synonyms from dbSNP
INSERT INTO sample_synonym (sample_id, source_id, name)
	SELECT ss.sample_id, s.source_id, ss.individual_id
	FROM source s, sample ss
	WHERE s.name = 'dbSNP'
	AND ss.individual_id is not null;

##and finally drop the extra columns in population and individual
ALTER TABLE sample DROP COLUMN population_id, DROP COLUMN individual_id;

##
## New table to store read coverage data
##

CREATE TABLE read_coverage (
   seq_region_id int not null,
   seq_region_start int not null,
   seq_region_end int not null,
   level tinyint not null,
   sample_id int not null,
		  
   key seq_region_idx(seq_region_id,seq_region_start)   
);

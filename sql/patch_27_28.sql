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


## Schema 27-28
## due to space restrictions, we have
## splitted the individual_genotype table
## in 2 different tables: one will contain
## genotypes that are single base paired, and 
## the other contains the rest of genotypes
## with the change in the column type for allele_1
## and allele_2 plus the removal of the individual_genotype_id
## column we save more than 4 Gb

create table individual_genotype_single_bp (
  variation_id int not null,
  allele_1 char,
  allele_2 char,
  individual_id int,

  key variation_idx(variation_id),
  key individual_idx(individual_id)
) MAX_ROWS = 100000000;

create table individual_genotype_multiple_bp (
  variation_id int not null,
  allele_1 varchar(255),
  allele_2 varchar(255),
  individual_id int,

  key variation_idx(variation_id),
  key individual_idx(individual_id)
);

## and now split the data from individual_genotype
## first populate the individual_genotype_single_bp

INSERT INTO individual_genotype_single_bp 
   SELECT individual_id, variation_id, allele_1, allele_2 FROM individual_genotype WHERE length(allele_1) = 1 AND length(allele_2) = 1;

## and then do the same for the multiple_bp table
INSERT INTO individual_genotype_multiple_bp 
   SELECT individual_id, variation_id, allele_1, allele_2 FROM individual_genotype WHERE length(allele_1) > 1 OR length(allele_2) > 1;

## and drop the individual_genotype table
DROP TABLE individual_genotype;

## another significant change is the relation between individual
## and population. Before the relation was 1 -> N (one individual
## only belonged to 1 population). In this schema, the relation is
## N -> N (one individual can belong to different populations)
## to implement this change, we have to add a new table, called
## individual_population

create table individual_population (
  individual_id int not null,
  population_id int not null,

  key individual_idx(individual_id),
  key population_idx(population_id)

);

## populate the new table

INSERT INTO individual_population (individual_id, population_id)
   SELECT individual_id, population_id
   FROM individual;

## and remove the population_id column from the table

ALTER TABLE individual DROP population_id;

##another minor change is the addition of a new column
## in the population table to deal with strains

ALTER TABLE population ADD (is_strain int(1) default 0 NOT NULL);

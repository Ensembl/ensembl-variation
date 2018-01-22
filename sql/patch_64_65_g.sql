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


# create allele_new table
CREATE TABLE IF NOT EXISTS `allele_new` (
  `allele_id` int(11) NOT NULL AUTO_INCREMENT,
  `variation_id` int(11) unsigned NOT NULL,
  `subsnp_id` int(11) unsigned DEFAULT NULL,
  `allele_code_id` int(11) unsigned NOT NULL,
  `sample_id` int(11) unsigned DEFAULT NULL,
  `frequency` float unsigned DEFAULT NULL,
  `count` int(11) unsigned DEFAULT NULL,
  PRIMARY KEY (`allele_id`),
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`),
  KEY `sample_idx` (`sample_id`)
);

# create population_genotype_new table
CREATE TABLE IF NOT EXISTS `population_genotype_new` (
  `population_genotype_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `variation_id` int(11) unsigned NOT NULL,
  `subsnp_id` int(11) unsigned DEFAULT NULL,
  `genotype_code_id` int(11) DEFAULT NULL,
  `frequency` float DEFAULT NULL,
  `sample_id` int(10) unsigned DEFAULT NULL,
  `count` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`population_genotype_id`),
  KEY `sample_idx` (`sample_id`),
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`)
);


# populate allele_new
alter table allele_new disable keys;

insert ignore into allele_new select a.allele_id, a.variation_id, a.subsnp_id, ac.allele_code_id, a.sample_id, a.frequency, a.count from allele a, allele_code ac where a.allele = ac.allele;

alter table allele_new enable keys;

# populate population_genotype_new
alter table population_genotype_new disable keys;

insert into population_genotype_new select pg.population_genotype_id, pg.variation_id, pg.subsnp_id, gc.genotype_code_id, pg.frequency, pg.sample_id, pg.count from population_genotype pg, genotype_code_tmp gc where pg.allele_1 = gc.allele_1 and pg.allele_2 = gc.allele_2;

alter table population_genotype_new enable keys;

# drop genotype_code_tmp
DROP TABLE genotype_code_tmp;

# rename old tables for Mart
RENAME TABLE allele TO MTMP_allele;
RENAME TABLE population_genotype TO MTMP_population_genotype;

# rename new tables
RENAME TABLE allele_new TO allele;
RENAME TABLE population_genotype_new TO population_genotype;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_64_65_g.sql|creates and populates new allele and population_genotype tables');

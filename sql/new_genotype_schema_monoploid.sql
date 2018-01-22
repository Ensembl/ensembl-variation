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


# create tables

SELECT "Creating tables", CURDATE(), CURTIME();

CREATE TABLE `subsnp_proxy` (
  `subsnp_proxy_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `variation_id` int(11) NOT NULL,
  `subsnp_id` int(11) DEFAULT NULL,
  PRIMARY KEY (`subsnp_proxy_id`),
  UNIQUE KEY `variation_idx` (`variation_id`,`subsnp_id`),
  KEY `subsnp_idx` (`subsnp_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS `allele_code` (
  `allele_code_id` int(11) NOT NULL AUTO_INCREMENT,
  `allele` varchar(60000) DEFAULT NULL,
  PRIMARY KEY (`allele_code_id`),
  UNIQUE KEY `allele_idx` (`allele`(1000))
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS `genotype_code_tmp` (
  `genotype_code_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `allele_1` varchar(60000) NOT NULL,
  PRIMARY KEY (`genotype_code_id`),
  UNIQUE KEY `genotype_idx` (`allele_1`(1000))
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS `genotype_code` (
  `genotype_code_id` int(11) unsigned NOT NULL,
  `allele_code_id` int(11) unsigned NOT NULL,
  `haplotype_id` tinyint(2) unsigned NOT NULL,
  KEY `genotype_code_id` (`genotype_code_id`),
  KEY `allele_code_id` (`allele_code_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS `allele_proxy` (
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
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS `population_genotype_proxy` (
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
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS `compressed_genotype_var` (
  `variation_id` int(11) unsigned NOT NULL,
  `subsnp_id` int(11) unsigned DEFAULT NULL,
  `genotypes` blob,
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS `compressed_genotype_region` (
  `sample_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(11) NOT NULL,
  `seq_region_end` int(11) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `genotypes` blob,
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


SELECT "Populating code tables", CURDATE(), CURTIME();

# add basic alleles to allele_code
INSERT INTO `allele_code` (`allele_code_id`,`allele`)
VALUES
	(1, 'T'),
	(2, 'A'),
	(3, 'G'),
	(4, 'C'),
	(5, '-'),
    (6, 'N');


# add basic genotypes to genotype_code_tmp
INSERT IGNORE INTO genotype_code_tmp (allele_1) SELECT allele FROM allele_code;

# add genotypes from genotype tables to genotype_code_tmp
INSERT IGNORE INTO genotype_code_tmp(allele_1) SELECT DISTINCT(allele_1) FROM individual_genotype_multiple_bp;
INSERT IGNORE INTO genotype_code_tmp(allele_1) SELECT DISTINCT(allele_1) FROM population_genotype;

# add alleles from genotype_code_tmp and allele
INSERT IGNORE INTO allele_code(allele) SELECT allele_1 FROM genotype_code_tmp;
INSERT IGNORE INTO allele_code(allele) SELECT DISTINCT(a.allele) FROM allele a LEFT JOIN allele_code ac ON a.allele = ac.allele WHERE ac.allele IS NULL;


SELECT "Populating genotype_code", CURDATE(), CURTIME();

# populate genotype_code
insert into genotype_code select t.genotype_code_id, ac.allele_code_id, 1 from genotype_code_tmp t, allele_code ac where t.allele_1 = ac.allele;
alter table genotype_code order by genotype_code_id, haplotype_id asc;


SELECT "Populating proxy tables", CURDATE(), CURTIME();


# populate allele_proxy
alter table allele_proxy disable keys;

insert ignore into allele_proxy select a.allele_id, a.variation_id, a.subsnp_id, ac.allele_code_id, a.sample_id, a.frequency, a.count from allele a, allele_code ac where a.allele = ac.allele;

alter table allele_proxy enable keys;


# populate population_genotype_proxy
alter table population_genotype_proxy disable keys;

insert into population_genotype_proxy select pg.population_genotype_id, pg.variation_id, pg.subsnp_id, gc.genotype_code_id, pg.frequency, pg.sample_id, pg.count from population_genotype pg, genotype_code_tmp gc where pg.allele_1 = gc.allele_1;

alter table population_genotype_proxy enable keys;

SELECT "DONE", CURDATE(), CURTIME();
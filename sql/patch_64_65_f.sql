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


# add allele_code table
CREATE TABLE IF NOT EXISTS `allele_code` (
  `allele_code_id` int(11) NOT NULL AUTO_INCREMENT,
  `allele` varchar(60000) DEFAULT NULL,
  PRIMARY KEY (`allele_code_id`),
  UNIQUE KEY `allele_idx` (`allele`(1000))
);

# add genotype_code_tmp table
CREATE TABLE IF NOT EXISTS `genotype_code_tmp` (
  `genotype_code_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `allele_1` varchar(30000) NOT NULL,
  `allele_2` varchar(30000) NOT NULL,
  PRIMARY KEY (`genotype_code_id`),
  UNIQUE KEY `genotype_idx` (`allele_1`(500),`allele_2`(500))
);

# add genotype_code table
CREATE TABLE IF NOT EXISTS `genotype_code` (
  `genotype_code_id` int(11) unsigned NOT NULL,
  `allele_code_id` int(11) unsigned NOT NULL,
  `haplotype_id` tinyint(2) unsigned NOT NULL,
  KEY `genotype_code_id` (`genotype_code_id`),
  KEY `allele_code_id` (`allele_code_id`)
);

# populate genotype_code_tmp from genotype_code if there are already entries
TRUNCATE genotype_code_tmp;

INSERT IGNORE INTO genotype_code_tmp (genotype_code_id, allele_1, allele_2)
SELECT gc1.genotype_code_id, ac1.allele, ac2.allele
FROM genotype_code gc1, genotype_code gc2, allele_code ac1, allele_code ac2
WHERE gc1.allele_code_id = ac1.allele_code_id
AND gc2.allele_code_id = ac2.allele_code_id
AND gc1.haplotype_id < gc2.haplotype_id;

# populate allele_code entries
INSERT IGNORE INTO `allele_code` (`allele_code_id`,`allele`)
VALUES
	(1, 'T'),
	(2, 'A'),
	(3, 'G'),
	(4, 'C'),
	(5, '-'),
    (6, 'N');
    
# add basic genotypes to genotype_code_tmp
INSERT IGNORE INTO genotype_code_tmp (allele_1, allele_2) SELECT ac1.allele, ac2.allele FROM allele_code ac1, allele_code ac2 ORDER BY ac1.allele_code_id, ac1.allele_code_id + ac2.allele_code_id;

# add genotypes from genotype tables to genotype_code_tmp
INSERT IGNORE INTO genotype_code_tmp(allele_1, allele_2) SELECT allele_1, allele_2 FROM individual_genotype_multiple_bp;
INSERT IGNORE INTO genotype_code_tmp(allele_1, allele_2) SELECT allele_1, allele_2 FROM population_genotype;

# add alleles from genotype_code_tmp and allele
INSERT IGNORE INTO allele_code(allele) SELECT allele_1 FROM genotype_code_tmp;
INSERT IGNORE INTO allele_code(allele) SELECT allele_2 FROM genotype_code_tmp;
INSERT IGNORE INTO allele_code(allele) SELECT a.allele FROM allele a LEFT JOIN allele_code ac ON a.allele = ac.allele WHERE ac.allele IS NULL;

# truncate and re-populate genotype_code
TRUNCATE genotype_code;

INSERT INTO genotype_code
SELECT t.genotype_code_id, ac.allele_code_id, 1
FROM genotype_code_tmp t, allele_code ac
WHERE t.allele_1 = ac.allele;

INSERT INTO genotype_code
SELECT t.genotype_code_id, ac.allele_code_id, 2
FROM genotype_code_tmp t, allele_code ac
WHERE t.allele_2 = ac.allele;

ALTER TABLE genotype_code ORDER BY genotype_code_id, haplotype_id ASC;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_64_65_f.sql|creates and populates allele_code and genotype_code');

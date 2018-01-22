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


# create new genotype tables
CREATE TABLE IF NOT EXISTS `compressed_genotype_var` (
  `variation_id` int(11) unsigned NOT NULL,
  `subsnp_id` int(11) unsigned DEFAULT NULL,
  `genotypes` blob,
  KEY `variation_idx` (`variation_id`),
  KEY `subsnp_idx` (`subsnp_id`)
);

CREATE TABLE IF NOT EXISTS `compressed_genotype_region` (
  `sample_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(11) NOT NULL,
  `seq_region_end` int(11) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `genotypes` blob,
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`),
  KEY `sample_idx` (`sample_id`)
);

# drop old compressed table
DROP TABLE compressed_genotype_single_bp;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_64_65_h.sql|creates new genotype tables');

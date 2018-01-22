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


# create table if it doesn't exist first
CREATE TABLE IF NOT EXISTS tmp_individual_genotype_single_bp (
	variation_id int(10) not null,
	subsnp_id int(15) unsigned,   
	allele_1 char(1),allele_2 char(1),sample_id int,
	key variation_idx(variation_id),
	key subsnp_idx(subsnp_id),
	key sample_idx(sample_id)
) MAX_ROWS = 100000000;

# make allele fields char(1) instead of varchar(255)
ALTER TABLE `tmp_individual_genotype_single_bp` CHANGE `allele_1` `allele_1` char(1) NULL DEFAULT NULL;
ALTER TABLE `tmp_individual_genotype_single_bp` CHANGE `allele_2` `allele_2` char(1) NULL DEFAULT NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_c.sql|change allele columns in tmp_individual_genotype_single_bp to char(1)');

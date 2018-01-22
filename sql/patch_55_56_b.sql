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


# added subsnp_id in the  variation_synonym, allele, population_genotype, tmp_individual_genotype_single_bp, individual_genotype_multiple_bp table
##################
alter table variation_synonym add subsnp_id int(15) unsigned after variation_id, add key subsnp_idx(subsnp_id);
alter table allele add subsnp_id int(15) unsigned after variation_id, add key subsnp_idx(subsnp_id);
alter table population_genotype add subsnp_id int(15) unsigned after variation_id, add key subsnp_idx(subsnp_id);
alter table tmp_individual_genotype_single_bp add subsnp_id int(15) unsigned after variation_id, add key subsnp_idx(subsnp_id);
alter table individual_genotype_multiple_bp add subsnp_id int(15) unsigned after variation_id, add key subsnp_idx(subsnp_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_55_56_b.sql|add subsnp_id');

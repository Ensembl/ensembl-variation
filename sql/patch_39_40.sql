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


####################
#  All ids are going to be int(10) unsigned for consistency and speed purposes
###################

alter table variation modify variation_id int(10) unsigned not null auto_increment, modify source_id int(10) unsigned not null;
alter table variation_synonym modify variation_synonym_id int(10) unsigned not null auto_increment, modify variation_id int(10) unsigned not null, modify source_id int(10) unsigned not null;
alter table sample_synonym modify sample_synonym_id int(10) unsigned not null auto_increment, modify sample_id int(10) unsigned not null, modify source_id int(10) unsigned not null;
alter table allele modify allele_id int(10) unsigned not null auto_increment, modify variation_id int(10) unsigned not null, modify sample_id int(10) unsigned;
alter table sample modify sample_id int(10) unsigned not null auto_increment;
alter table population modify sample_id int(10) unsigned not null;
alter table population_structure modify super_population_sample_id int(10) unsigned not null, modify sub_population_sample_id int(10) unsigned not null;
alter table individual modify sample_id int(10) unsigned not null, modify father_individual_sample_id int(10) unsigned, modify mother_individual_sample_id int(10) unsigned;
alter table variation_feature modify variation_feature_id int(10) unsigned not null auto_increment, modify seq_region_id int(10) unsigned not null, modify variation_id int(10) unsigned not null, modify source_id int(10) unsigned not null;
alter table variation_group modify variation_group_id int(10) unsigned not null auto_increment, modify source_id int(10) unsigned not null;
alter table variation_group_variation modify variation_id int(10) unsigned not null, modify variation_group_id int(10) unsigned not null;
alter table variation_group_feature modify variation_group_feature_id int(10) unsigned not null auto_increment, modify seq_region_id int(10) unsigned not null, modify variation_group_id int(10) unsigned not null;
alter table transcript_variation modify transcript_variation_id int(10) unsigned not null auto_increment, modify transcript_id int(10) unsigned not null, modify variation_feature_id int(10) unsigned not null;
alter table allele_group modify allele_group_id int(10) not null auto_increment, modify variation_group_id int(10) unsigned not null, modify sample_id int(10) unsigned, modify source_id int(10) unsigned;
alter table allele_group_allele modify allele_group_id int(10) unsigned not null, modify variation_id int(10) unsigned not null;
alter table flanking_sequence modify variation_id int(10) unsigned not null, modify seq_region_id int(10) unsigned;
alter table httag modify httag_id int(10) unsigned not null auto_increment, modify variation_group_id int(10) unsigned not null, modify source_id int(10) unsigned;
alter table source modify source_id int(10) unsigned not null auto_increment;
alter table population_genotype modify population_genotype_id int(10) unsigned not null auto_increment, modify variation_id int(10) unsigned not null, modify sample_id int(10) unsigned;
alter table individual_population modify individual_sample_id int(10) unsigned not null, modify population_sample_id int(10) unsigned not null;
alter table individual_genotype_multiple_bp modify variation_id int(10) unsigned not null, modify sample_id int(10) unsigned;
alter table meta_coord modify coord_system_id int(10) unsigned not null;
alter table meta modify meta_id int(10) unsigned not null auto_increment;
alter table tagged_variation_feature modify variation_feature_id int(10) unsigned not null, modify sample_id int(10) unsigned not null;
alter table read_coverage modify seq_region_id int(10) unsigned not null, modify sample_id int(10) unsigned not null;
alter table compressed_genotype_single_bp modify sample_id int(10) unsigned not null, modify seq_region_id int(10) unsigned not null;


###################
# update the schema_version entry in the meta table
##################
update meta set meta_value = 40 where meta_key = 'schema_version';

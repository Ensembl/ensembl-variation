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


# patch_52_53_b.sql
#
# title: change variation_annotation table
#
# description:
# Table containing annotation associated with the variation
# such as GWAS

alter table variation_annotation 
add column study varchar(30) default NULL after source_id,
add column associated_gene varchar(255) default NULL,
add column associated_variant_risk_allele varchar(255) default NULL,
add column variation_names varchar(255) default NULL,
add column risk_allele_freq_in_controls varchar(30) default NULL,
add column p_value varchar(20) default NULL;




INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_52_53_b.sql|add columns in variation_annotation table');

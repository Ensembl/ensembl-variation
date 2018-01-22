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


# patch_51_52_b.sql
#
# title: add variation_annotation and phenotype tables
#
# description:
# Table containing annotation associated with the variation
# such as GWAS

create table variation_annotation (
        variation_annotation_id int(10) unsigned not null auto_increment,
        variation_id int(10) unsigned not null,
        phenotype_id int(10) unsigned not null,
        source_id int(10) unsigned not null,
        study_type set('GWAS'),
        local_stable_id varchar(255),
        primary key (variation_annotation_id),
        key variation_idx(variation_id),
        key phenotype_idx(phenotype_id),
        key source_idx(source_id)
);

create table phenotype (
        phenotype_id int(10) unsigned not null auto_increment,
        name varchar(50),
        description varchar(255),

        primary key (phenotype_id),
        unique key name_idx(name)
);

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_51_52_c.sql|add variation_annotation/phenotype table');

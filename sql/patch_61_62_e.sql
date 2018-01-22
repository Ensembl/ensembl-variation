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


# add the tables study and supporting_structural_variation

create table study (
	study_id int(10) unsigned not null auto_increment,
	source_id int(10) unsigned not null,
	name varchar(255),
	description varchar(255),
	url varchar(255),
	external_reference varchar(255),
	study_type set('GWAS'),
	
	primary key( study_id ),
	key source_idx (source_id)
);


create table supporting_structural_variation (
	supporting_structural_variation_id int(10) unsigned not null auto_increment,
	name varchar(255),
	structural_variation_id int(10) unsigned not null,
	
	primary key( supporting_structural_variation_id ),
	unique key name_idx(name),
	key structural_variation_idx (structural_variation_id)
);


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_e.sql|add the tables study and supporting_structural_variation');

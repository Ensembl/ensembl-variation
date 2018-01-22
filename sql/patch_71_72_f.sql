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


# store limited details on publications and link to variations

CREATE TABLE publication(
publication_id int(10) unsigned not null auto_increment, 
title          varchar(255),
authors        varchar(255),
pmid           int(10),
pmcid          varchar(255),
primary key( publication_id ),
key pmid_idx (pmid)
);

CREATE TABLE variation_citation (
   variation_id int(10) unsigned not null,
   publication_id int(10) unsigned not null,
   PRIMARY KEY variation_citation_idx (variation_id, publication_id)
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_71_72_f.sql|new tables for citations');

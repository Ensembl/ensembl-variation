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


# add top level menus
insert into meta(species_id, meta_key, meta_value)
values
(NULL, 'web_config', "menu#dbSNP#dbsnp#"),
(NULL, 'web_config', "menu#1000 Genomes & HapMap#1kg_hapmap#"),
(NULL, 'web_config', "menu#Phenotype and curated variants#phenotype#"),
(NULL, 'web_config', "menu#Individual genomes#ind_genomes#"),
(NULL, 'web_config', "menu#Arrays and other#var_other#"),
(NULL, 'web_config', "menu#Failed variants#failed#");

# add all variants
insert into meta(species_id, meta_key, meta_value)
values
(NULL, 'web_config', "source#Sequence variants (dbSNP and all other sources)#variation_feature_variation#dbsnp");

# add dbSNP
insert into meta(species_id, meta_key, meta_value)
values
(NULL, 'web_config', "source#dbSNP variants#variation_feature_variation_dbSNP#dbsnp");


# add HapMap sets
insert into meta(species_id, meta_key, meta_value)
values
(NULL, 'web_config', "menu#HapMap#hapmap#1kg_hapmap");

insert into meta(species_id, meta_key, meta_value)
select NULL, "web_config",  concat("set#", vs.name, "#variation_set_", a.value, "#hapmap")
from variation_set vs, attrib a
where vs.short_name_attrib_id = a.attrib_id
and a.value like 'hapmap%';

# add 1000 genomes sets
insert into meta(species_id, meta_key, meta_value)
values
(NULL, 'web_config', "menu#1000 Genomes#1kg#1kg_hapmap");

insert into meta(species_id, meta_key, meta_value)
select NULL, "web_config",  concat("set#", vs.name, "#variation_set_", a.value, "#1kg")
from variation_set vs, attrib a
where vs.short_name_attrib_id = a.attrib_id
and a.value like '1kg%'
and a.value != '1kg_hq'
order by vs.name asc;

# add phenotype sets
insert into meta(species_id, meta_key, meta_value)
select NULL, "web_config",  concat("set#", vs.name, "#variation_set_", a.value, "#phenotype")
from variation_set vs, attrib a
where vs.short_name_attrib_id = a.attrib_id
and a.value like 'ph%'
and vs.name not like '%cosmic%';

# add phenotype sources
insert into meta(species_id, meta_key, meta_value)
values
(NULL, 'web_config', "menu#LSDB-associated variants#lsdb#phenotype");

insert into meta(species_id, meta_key, meta_value)
select NULL, "web_config",  concat("source#", name, " variants", "#variation_feature_variation_", replace(name, " ", "_"), "#lsdb")
from source
where name like 'LSDB%';

# add clinical set (phenotype)
insert into meta(species_id, meta_key, meta_value)
select NULL, "web_config",  concat("set#", vs.name, "#variation_set_", a.value, "#phenotype")
from variation_set vs, attrib a
where vs.short_name_attrib_id = a.attrib_id
and a.value like 'precious%';

# add individual genotype sets
insert into meta(species_id, meta_key, meta_value)
select NULL, "web_config",  concat("set#", vs.name, "#variation_set_", a.value, "#ind_genomes")
from variation_set vs, attrib a
where vs.short_name_attrib_id = a.attrib_id
and a.value like 'ind%'
order by vs.name asc;

# add failed sets
insert into meta(species_id, meta_key, meta_value)
select NULL, "web_config",  concat("set#", vs.name, "#variation_set_", a.value, "#failed")
from variation_set vs, attrib a
where vs.short_name_attrib_id = a.attrib_id
and a.value like 'fail%'
order by vs.variation_set_id asc;

# add chips and uniprot
insert into meta(species_id, meta_key, meta_value)
select NULL, "web_config",  concat("source#", name, " variants", "#variation_feature_variation_",
replace(name, " ", "_"), "#var_other")
from source
where name like 'Affy%'
or name like 'Illumina%'
or name like 'Uniprot%'
order by name asc;

# add structural variation sets
insert into meta(species_id, meta_key, meta_value)
values (NULL, 'web_config', "sv_set#1000 Genomes - Low coverage#sv_set_1kg_lc");

insert into meta(species_id, meta_key, meta_value)
values (NULL, 'web_config', "sv_set#1000 Genomes - High coverage - Trios#sv_set_1kg_hct");

insert into meta(species_id, meta_key, meta_value)
values (NULL, 'web_config', "sv_set#1000 Genomes - All#sv_set_1kg");

insert into meta(species_id, meta_key, meta_value)
values (NULL, 'web_config', "sv_set#1000 Genomes - EUR#sv_set_1kg_eur");

insert into meta(species_id, meta_key, meta_value)
values (NULL, 'web_config', "sv_set#1000 Genomes - ASN#sv_set_1kg_asn");

insert into meta(species_id, meta_key, meta_value)
values (NULL, 'web_config', "sv_set#1000 Genomes - AMR#sv_set_1kg_amr");

insert into meta(species_id, meta_key, meta_value)
values (NULL, 'web_config', "sv_set#1000 Genomes - AFR#sv_set_1kg_afr");

insert into meta(species_id, meta_key, meta_value) 
values (NULL, 'web_config', "sv_set#1000 Genomes - High quality#sv_set_1kg_hq");

# add structural variation studies (selection of)
insert into meta(species_id, meta_key, meta_value) values (NULL, 'web_config', "sv_study#Study estd199#estd199");
insert into meta(species_id, meta_key, meta_value) values (NULL, 'web_config', "sv_study#Study nstd37#nstd37");
insert into meta(species_id, meta_key, meta_value) values (NULL, 'web_config', "sv_study#Study estd203#estd203");
insert into meta(species_id, meta_key, meta_value) values (NULL, 'web_config', "sv_study#Study nstd71#nstd71");

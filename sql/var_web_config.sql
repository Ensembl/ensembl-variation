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

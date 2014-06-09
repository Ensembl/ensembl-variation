
## Add tables used in website index building to schema

create table if not exists variation_hgvs(
variation_id int(10) unsigned not null,
hgvs_name varchar(255) not null,
primary key(variation_id, hgvs_name));

create table if not exists variation_genename (
variation_id int(10) unsigned not null, 
gene_name varchar(255) not null, 
primary key(variation_id, gene_name));

#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_i.sql|Add tables required for HGVS index creation');

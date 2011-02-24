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

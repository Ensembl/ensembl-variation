# add the table associate_study

create table associate_study (
	study1_id int(10) unsigned not null,
	study2_id int(10) unsigned not null,
	
	primary key( study1_id,study2_id )
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_n.sql|add the table associate_study');

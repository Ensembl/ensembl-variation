# add the table associate_study

create table associate_study (
	study_id_1 int(10) unsigned not null,
	study_id_2 int(10) unsigned not null,
	
	primary key( study_id_1,study_id_2 ),
	key study_idx (study_id_1,study_id_2)
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_n.sql|add the table associate_study');

## create link between variation and study to allow PubMed ids to be stored against variants


create table study_variation (
   variation_id int(10) unsigned not null,
   study_id int(10) unsigned not null,
   primary key study_variation_idx (variation_id, study_id)
);


INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_e.sql|create study_variation table');

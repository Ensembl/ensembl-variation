
## group populations for display on PopulationGenetics page
## 

create table display_group(
  display_group_id int(10) unsigned not null auto_increment ,
  display_priority int(10) unsigned not null, 
  display_name     varchar(255) not null,

	primary key( display_group_id ),
	unique ( display_name ),
	unique ( display_priority )

 );

ALTER TABLE population ADD column display_group_id tinyint(1) default null;



#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_c.sql|Add new table and extra column to population table to specify if population is to be displayed seperately on the PopulationGenetics page and if so with what priority');

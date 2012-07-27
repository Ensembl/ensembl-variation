# add more support for dbSNP submitter handles
##################


ALTER TABLE allele MODIFY COLUMN frequency_submitter_handle int(10);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_68_69_c.sql|update for dbSNP handle information');

create table submitter_handle (
  handle_id int(10) unsigned not null auto_increment, # PK 
  handle varchar(25),
 primary key( handle_id ),
        unique ( handle )
);
 

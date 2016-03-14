CREATE TABLE sample_synonym (
  synonym_id int(10) unsigned not null auto_increment,
  sample_id int(10) unsigned not null,
  source_id int(10) unsigned not null,
  name varchar(255),
  primary key(synonym_id),
  key sample_idx (sample_id),
  key (name, source_id)
);

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_b.sql|create sample_synonym');

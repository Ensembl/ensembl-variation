# Replace sample_synonym by individual_synonym and population_synonym

CREATE TABLE individual_synonym (
  synonym_id int(10) unsigned not null auto_increment,
  individual_id int(10) unsigned not null,
  source_id int(10) unsigned not null,
  name varchar(255),

  primary key(synonym_id),
  key individual_idx (individual_id),
  key (name, source_id)
);

INSERT INTO individual_synonym(individual_id, source_id, name)
SELECT i.individual_id, s.source_id, s.name 
FROM individual i, sample_synonym s
WHERE i.individual_id = s.sample_id; 

CREATE TABLE population_synonym (
  synonym_id int(10) unsigned not null auto_increment,
  population_id int(10) unsigned not null,
  source_id int(10) unsigned not null,
  name varchar(255),

  primary key(synonym_id),
  key population_idx (population_id),
  key (name, source_id)
);

INSERT INTO population_synonym(population_id, source_id, name)
SELECT p.population_id, s.source_id, s.name
FROM population p, sample_synonym s
WHERE p.population_id = s.sample_id;

DROP TABLE sample_synonym;
DROP TABLE sample;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_71_72_d.sql|Replace sample_synonym by individual_synonym and population_synonym');

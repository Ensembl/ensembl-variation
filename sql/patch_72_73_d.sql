## revert population synonym and individual synonym to varchars

alter table individual_synonym modify column  name varchar(255);

alter table population_synonym modify column  name varchar(255);

##patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_72_73_d.sql| revert population_synonym.name and individual_synonym.name to varchars');

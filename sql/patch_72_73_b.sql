## Add year to publication table

alter table publication add (year  int(10) unsigned);


##patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_72_73_b.sql|Add year to publication table');

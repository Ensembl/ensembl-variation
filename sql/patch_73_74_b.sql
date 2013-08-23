## Add doi to publication table


alter table publication add column doi varchar(50) DEFAULT NULL;
CREATE INDEX doi_idx on publication (doi); 

#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_b.sql|Add doi to publication table');

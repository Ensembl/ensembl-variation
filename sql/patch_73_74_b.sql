## Add doi and UCSC id to publication table


alter table publication add column doi varchar(50) DEFAULT NULL;
CREATE INDEX doi_idx on publication (doi); 

alter table publication add column ucsc_id varchar(50) DEFAULT NULL;

#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_b.sql|Add doi and UCSC id to publication table');

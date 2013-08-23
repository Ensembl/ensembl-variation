## Add phenotype_citation table


CREATE TABLE phenotype_citation (
   phenotype_id int(11) unsigned not null,
   publication_id int(10) unsigned not null,
   PRIMARY KEY phenotype_citation_idx (phenotype_id, publication_id)
);



#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_c.sql|Add phenotype_citation table');

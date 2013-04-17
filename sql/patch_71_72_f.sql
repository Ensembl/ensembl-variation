# store limited details on publications and link to variations

CREATE TABLE publication(
publication_id int(10) unsigned not null auto_increment, 
title          varchar(255),
authors        varchar(255),
pmid           int(10),
pmcid          varchar(255),
primary key( publication_id ),
key pmid_idx (pmid)
);

CREATE TABLE variation_citation (
   variation_id int(10) unsigned not null,
   publication_id int(10) unsigned not null,
   PRIMARY KEY variation_citation_idx (variation_id, publication_id)
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_71_72_f.sql|new tables for citations');

# extend the characters supported in the publication.authors column

alter table publication modify authors varchar(255) CHARACTER SET latin2 ;  

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_90_91_b.sql|extend the characters supported in the publication.authors column'); 

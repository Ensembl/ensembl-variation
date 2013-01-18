## increase size of minor_allele in variation & variation_feature to allow storage of frequency information for small indels

alter table variation modify column minor_allele varchar(10);

alter table variation_feature modify column minor_allele varchar(10);

# patch identifier
#INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_c.sql|increase size of minor_allele');

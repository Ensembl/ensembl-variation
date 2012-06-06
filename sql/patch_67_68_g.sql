# update to description text in the failed_description table
##################


update  failed_description  set description = 'Additional allele data from dbSNP does not match the refSNP alleles' where failed_description_id = 11 ;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_67_68_g.sql|update to description text in the failed_description table');



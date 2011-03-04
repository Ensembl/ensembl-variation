# add table for failed_alleles
CREATE TABLE failed_allele (
  failed_allele_id int(11) NOT NULL AUTO_INCREMENT,
  allele_id int(10) unsigned NOT NULL,
  failed_description_id int(10) unsigned NOT NULL,
  PRIMARY KEY (failed_allele_id),
  UNIQUE KEY allele_idx (allele_id,failed_description_id)
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_c.sql|add failed_allele table');

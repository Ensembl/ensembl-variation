CREATE TABLE ancestral_alleles (
  ancestral_allele_id int(11) NOT NULL AUTO_INCREMENT,
  variation_id int(10) unsigned NOT NULL,
  ancestral_allele varchar(255) DEFAULT NULL,
  PRIMARY KEY (ancestral_allele_id)
);

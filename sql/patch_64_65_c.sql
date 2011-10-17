# adds a failed table for structural variation

CREATE TABLE failed_structural_variation (
  failed_structural_variation_id int(11) NOT NULL AUTO_INCREMENT,
  structural_variation_id int(10) unsigned NOT NULL,
  failed_description_id int(10) unsigned NOT NULL,
	
  PRIMARY KEY (failed_structural_variation_id),
  UNIQUE KEY structural_variation_idx (structural_variation_id,failed_description_id)
);

INSERT INTO failed_description (failed_description_id,description) VALUES (17,'Variation can not be re-mapped to the current assembly');
INSERT INTO failed_description (failed_description_id,description) VALUES (18,'Supporting evidence can not be re-mapped to the current assembly');


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_64_65_c.sql|adds a failed table for structural variation');

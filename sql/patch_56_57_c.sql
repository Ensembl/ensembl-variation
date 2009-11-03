# add a table to hold structural variations
CREATE TABLE structural_variation (
  structural_variation_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  seq_region_id int(10) unsigned NOT NULL,
  seq_region_start int(11) NOT NULL,
  seq_region_end int(11) NOT NULL,
  seq_region_strand tinyint(4) NOT NULL,
  variation_name varchar(255) DEFAULT NULL,
  source_id int(10) unsigned NOT NULL,
  class varchar(255) DEFAULT NULL,
  bound_start int(11) DEFAULT NULL,
  bound_end int(11) DEFAULT NULL,
  PRIMARY KEY (structural_variation_id),
  KEY pos_idx (seq_region_id,seq_region_start)
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_56_57_c.sql|add structural variation');

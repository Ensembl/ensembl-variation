# add a table to hold structural variation features
CREATE TABLE structural_variation_feature (
  structural_variation_feature_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  seq_region_id int(10) unsigned NOT NULL,
  seq_region_start int(11) NOT NULL,
  seq_region_end int(11) NOT NULL,
  seq_region_strand tinyint(4) NOT NULL,
  variation_id int(10) unsigned NOT NULL,
  variation_name varchar(255) DEFAULT NULL,
  map_weight int(11) NOT NULL,
  source_id int(10) unsigned NOT NULL,
  class varchar(255) DEFAULT NULL,
  bound_start int(11) DEFAULT NULL,
  bound_end int(11) DEFAULT NULL,
  PRIMARY KEY (structural_variation_feature_id),
  KEY pos_idx (seq_region_id,seq_region_start),
  KEY variation_idx (variation_id)
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_56_57_b.sql|schema version');

# create variation_attrib table
CREATE TABLE variation_attrib (
  variation_id INT(11) UNSIGNED NOT NULL,
  attrib_type_id INT(11) DEFAULT NULL,
  value VARCHAR(255) DEFAULT NULL,
  KEY variation_idx (variation_id),
  KEY type_value_idx (attrib_type_id,value)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_76_77_e.sql|add variation_attrib table');

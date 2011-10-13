# adds a variation set table for structural variation

CREATE TABLE IF NOT EXISTS variation_set_structural_variation (
	structural_variation_id int(10) unsigned NOT NULL,
	variation_set_id int(10) unsigned NOT NULL,
	
	PRIMARY KEY (structural_variation_id,variation_set_id)
);


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_64_65_d.sql|adds a variation set table for structural variation');

# structural variation schema changes

DROP TABLE supporting_structural_variation;


ALTER TABLE structural_variation ADD COLUMN is_evidence TINYINT(4) DEFAULT 0;

ALTER TABLE structural_variation ADD INDEX source_idx (source_id);
ALTER TABLE structural_variation_feature ADD INDEX source_idx (source_id);

CREATE TABLE structural_variation_association (
  structural_variation_id int(10) unsigned NOT NULL,
  supporting_structural_variation_id int(10) unsigned NOT NULL,
	
  PRIMARY KEY (structural_variation_id, supporting_structural_variation_id)
);


CREATE TABLE structural_variation_annotation (
	structural_variation_annotation_id int(10) unsigned NOT NULL auto_increment,
	structural_variation_id int(10) unsigned NOT NULL,
	clinical_attrib_id int(10) unsigned DEFAULT NULL,
	phenotype_id int(10) unsigned DEFAULT NULL,
	sample_id int(10) unsigned DEFAULT NULL,
	strain_id int(10) unsigned DEFAULT NULL,
	
	primary key (structural_variation_annotation_id),
	key structural_variation_idx(structural_variation_id),
	key clinical_attrib_idx(clinical_attrib_id),
	key phenotype_idx(phenotype_id),
	key sample_idx(sample_id),
	key strain_idx(strain_id)
);


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_64_65_b.sql|structural variation schema changes');

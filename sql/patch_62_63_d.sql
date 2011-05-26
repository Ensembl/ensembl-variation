# alter some structural_variation columns and a supporting_structural_variation index

ALTER TABLE structural_variation CHANGE bound_start inner_start int(11) DEFAULT NULL;
ALTER TABLE structural_variation CHANGE bound_end inner_end int(11) DEFAULT NULL;
ALTER TABLE structural_variation ADD COLUMN class_attrib_id int(10) unsigned not null default 0;
ALTER TABLE structural_variation DROP COLUMN class;
ALTER TABLE structural_variation ADD COLUMN validation_status ENUM('validated','not validated','high quality');
ALTER TABLE structural_variation ADD INDEX attrib_idx (class_attrib_id);

ALTER TABLE supporting_structural_variation DROP INDEX name_idx;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_62_63_d.sql|alter some structural_variation columns and a supporting_structural_variation index');

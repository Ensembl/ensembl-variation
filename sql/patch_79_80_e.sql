

# Update the attrib tables by changing the default values
##################
ALTER TABLE attrib_type CHANGE code code VARCHAR(20) NOT NULL DEFAULT '';
ALTER TABLE attrib_type CHANGE name name VARCHAR(255) NOT NULL DEFAULT '';

ALTER TABLE attrib CHANGE attrib_type_id attrib_type_id SMALLINT(5) UNSIGNED NOT NULL DEFAULT 0;
ALTER TABLE attrib CHANGE value value TEXT NOT NULL;

ALTER TABLE attrib_set CHANGE attrib_set_id attrib_set_id INT(11) UNSIGNED NOT NULL DEFAULT 0;
ALTER TABLE attrib_set CHANGE attrib_id attrib_id INT(11) UNSIGNED NOT NULL DEFAULT 0;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_79_80_e.sql|update the attrib tables by changing the default values');

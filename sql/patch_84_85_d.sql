-- Making attrib_id auto_increment
ALTER TABLE attrib MODIFY COLUMN attrib_id int(11) unsigned auto_increment;

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_d.sql|Making attrib_id auto_increment');

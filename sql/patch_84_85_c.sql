-- drop column moltype from variation_synonym
ALTER TABLE variation_synonym DROP COLUMN moltype;

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_c.sql|drop column moltype from variation_synonym');

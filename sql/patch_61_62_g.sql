# move the somatic flag from the source to variation table

ALTER TABLE source DROP COLUMN somatic;
ALTER TABLE variation ADD COLUMN somatic tinyint(1);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_g.sql|move somatic flag from source to variation table');

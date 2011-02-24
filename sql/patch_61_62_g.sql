# move the somatic flag from the source to variation and variation_feature tables

ALTER TABLE source DROP COLUMN somatic;
ALTER TABLE variation ADD COLUMN somatic tinyint(1) DEFAULT 0 NOT NULL;
ALTER TABLE variation_feature ADD COLUMN somatic tinyint(1) DEFAULT 0 NOT NULL;

# patch identifier

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_g.sql|move somatic flag from source to variation and variation_feature tables');

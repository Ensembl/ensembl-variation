# convert column of type text to varchar in motif_feature_variation table

ALTER TABLE motif_feature_variation MODIFY motif_name VARCHAR(60);

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_80_81_d.sql|Update type of motif_name to varchar.');

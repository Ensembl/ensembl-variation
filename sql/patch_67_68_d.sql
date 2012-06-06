# add distance_to_transcript field to transcript_variation
ALTER TABLE transcript_variation ADD COLUMN distance_to_transcript INT(11) UNSIGNED NULL DEFAULT NULL AFTER translation_end;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_67_68_d.sql|add distance_to_transcript field to transcript_variation');
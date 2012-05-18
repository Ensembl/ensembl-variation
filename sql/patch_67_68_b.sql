# update the transcript_variation table
##################


ALTER TABLE transcript_variation CHANGE COLUMN hgvs_coding  hgvs_transcript text;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_67_68_b.sql|update the transcript_variation table');

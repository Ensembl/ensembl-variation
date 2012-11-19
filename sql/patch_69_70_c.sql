# add freqs_from_gts column to sample table
ALTER TABLE sample ADD freqs_from_gts BOOLEAN NOT NULL DEFAULT 0 AFTER display;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_69_70_c.sql|add freqs_from_gts column to sample table');

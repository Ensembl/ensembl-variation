# Add new column to transcript_variation table to flag whether variants should be displayed or not

ALTER TABLE transcript_variation ADD COLUMN display int(1) DEFAULT 1;

UPDATE transcript_variation SET display = 0 WHERE variation_feature_id in (SELECT variation_feature_id FROM variation_feature where display =0 );


#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_77_78_f.sql|Add new column to the transcript_variation table to flag whether variants should be displayed or not'
);

# Add new column to Variation_feature table to flag whether variants should be displayed or not

ALTER TABLE variation_feature ADD COLUMN display int(1) DEFAULT 1;

UPDATE variation_feature SET display = 0 WHERE variation_id in (SELECT failed_variation.variation_id FROM failed_variation
left outer join variation_citation on failed_variation.variation_id = variation_citation.variation_id 
where variation_citation.publication_id is null );


#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_76_77_g.sql|Add new column to Variation_feature table to flag whether variants should be displayed or not'
);

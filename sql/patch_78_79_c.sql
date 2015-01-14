
## store more detailed Sift information

CREATE TABLE protein_function_predictions_attrib (
    translation_md5_id int(11) unsigned NOT NULL,
    analysis_attrib_id int(11) unsigned NOT NULL,
    attrib_type_id     int(11) unsigned NOT NULL,
    position_values    blob,
    
    PRIMARY KEY (translation_md5_id, analysis_attrib_id,attrib_type_id )
);


alter table transcript_variation MODIFY COLUMN sift_prediction enum('tolerated', 'deleterious','tolerated - low confidence', 'deleterious - low confidence' ) DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_78_79_c.sql|Store more detailed Sift information');
# update to the new protein function prediction schema

DROP TABLE protein_info;
DROP TABLE protein_position;
DROP TABLE polyphen_prediction;
DROP TABLE sift_prediction;

CREATE TABLE protein_function_predictions (

    translation_stable_id   VARCHAR(128) NOT NULL,
    transcript_stable_id    VARCHAR(128) NOT NULL,
    translation_md5         CHAR(32) NOT NULL,
    polyphen_predictions    MEDIUMBLOB,
    sift_predictions        MEDIUMBLOB,
    
    PRIMARY KEY (translation_stable_id),
    KEY transcript_idx (transcript_stable_id)
);

# add sift and polyphen scores to the transcript_variation table for biomart

ALTER TABLE transcript_variation ADD COLUMN polyphen_score FLOAT DEFAULT NULL; 
ALTER TABLE transcript_variation ADD COLUMN sift_score FLOAT DEFAULT NULL; 

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_63_64_c.sql|update to new protein function prediction schema');

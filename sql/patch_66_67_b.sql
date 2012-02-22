DROP TABLE protein_function_predictions;

CREATE TABLE protein_function_predictions (

    translation_md5             CHAR(32) NOT NULL,
    sift_predictions            MEDIUMBLOB,
    polyphen_humdiv_predictions MEDIUMBLOB,
    polyphen_humvar_predictions MEDIUMBLOB,
    
    PRIMARY KEY (translation_md5)
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_66_67_b.sql|update protein function predictions table schema');

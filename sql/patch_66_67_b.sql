DROP TABLE protein_function_predictions;

CREATE TABLE protein_function_predictions (
    translation_md5_id int(11) unsigned NOT NULL,
    analysis_attrib_id int(11) unsigned NOT NULL,
    prediction_matrix mediumblob,
    
    PRIMARY KEY (translation_md5_id, analysis_attrib_id)
);

CREATE TABLE translation_md5 (
    translation_md5_id int(11) NOT NULL AUTO_INCREMENT,
    translation_md5 char(32) NOT NULL,

    PRIMARY KEY (translation_md5_id),
    UNIQUE KEY md5_idx (translation_md5)
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_66_67_b.sql|update protein function predictions table schema');

# change so_id to so_accession and change the class_so_id column to just class

DROP TABLE IF EXISTS variation_class;

CREATE TABLE variation_class (
    so_accession        VARCHAR(128) NOT NULL,
    so_term             VARCHAR(128),
    display_term        VARCHAR(128),
    
    PRIMARY KEY (so_accession)
);

INSERT INTO variation_class (
        so_accession,
        so_term,
        display_term
    )
    VALUES (
        'SO:0001483',
        'SNV',
        'SNP'
    ), (
        'SO:0001019',
        'copy_number_variaton',
        'CNV'
    ), (
        'SO:0000667',
        'insertion',
        'insertion'
    ), (
        'SO:0000159',
        'deletion',
        'deletion'
    ), (
        'SO:1000032',
        'indel',
        'indel'
    ), (
        'SO:0001059',
        'sequence_alteration',
        'sequence_alteration'
    ), (
        'SO:1000002',
        'substitution',
        'substitution'
    ), (
        'SO:0000705',
        'tandem_repeat',
        'tandem_repeat'
    )
;

ALTER TABLE variation CHANGE class_so_id class ENUM(
    'SO:0001483', # SNV
    'SO:1000002', # substitution
    'SO:0000667', # insertion
    'SO:0000159', # deletion
    'SO:0000705', # tandem_repeat
    'SO:1000032', # indel
    'SO:0001059', # sequence_alteration
    'SO:0001019'  # copy_number_variation
) DEFAULT 'SO:0001059';

ALTER TABLE variation_feature CHANGE class_so_id class ENUM(
    'SO:0001483', # SNV
    'SO:1000002', # substitution
    'SO:0000667', # insertion
    'SO:0000159', # deletion
    'SO:0000705', # tandem_repeat
    'SO:1000032', # indel
    'SO:0001059', # sequence_alteration
    'SO:0001019'  # copy_number_variation
) DEFAULT 'SO:0001059';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_h.sql|change so_id to so_accession and change the class_so_id column to just class');

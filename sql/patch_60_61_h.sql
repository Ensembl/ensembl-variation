# add class_so_id columns to the variation and variation_feature tables

ALTER TABLE variation ADD COLUMN class_so_id ENUM(
    'SO:0001483', # SNV
    'SO:1000002', # substitution
    'SO:0000667', # insertion
    'SO:0000159', # deletion
    'SO:0000705', # tandem_repeat
    'SO:1000032', # indel
    'SO:0001059', # sequence_alteration
    'SO:0001019'  # copy_number_variation
) DEFAULT 'SO:0001059';

ALTER TABLE variation_feature ADD COLUMN class_so_id ENUM(
    'SO:0001483', # SNV
    'SO:1000002', # substitution
    'SO:0000667', # insertion
    'SO:0000159', # deletion
    'SO:0000705', # tandem_repeat
    'SO:1000032', # indel
    'SO:0001059', # sequence_alteration
    'SO:0001019'  # copy_number_variation
) DEFAULT 'SO:0001059';


# create a table that maps the SO variation class ID to the ensembl display term and SO term

CREATE TABLE variation_class (
    so_id               VARCHAR(128) NOT NULL,
    so_term             VARCHAR(128),
    display_term        VARCHAR(128),
    
    PRIMARY KEY (so_id)
);

INSERT INTO variation_class (
        so_id,
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
        'indel'
    ), (
        'SO:0000159',
        'deletion',
        'indel'
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


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_h.sql|add support for storing the class of a variation in the database using the SO id');


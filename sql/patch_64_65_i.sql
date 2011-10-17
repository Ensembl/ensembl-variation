# add some new columns to the variation table to store new dbSNP 134 stuff

ALTER TABLE variation 
    ADD COLUMN minor_allele char(1) DEFAULT NULL,
    ADD COLUMN minor_allele_freq float DEFAULT NULL,
    ADD COLUMN minor_allele_count int(10) unsigned DEFAULT NULL, 
    ADD COLUMN clinical_significance_attrib_id int(10) unsigned DEFAULT NULL
;

# and add a new failed_description

INSERT INTO failed_description (failed_description_id,description) VALUES (16,'Flagged as suspect by dbSNP');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_64_65_i.sql|add support for new data types from dbSNP');


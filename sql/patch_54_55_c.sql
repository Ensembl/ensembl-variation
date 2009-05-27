# update the description entry in the failed_description table
# better wording
##################
truncate table failed_description;

#possible values in the failed_description table
INSERT INTO failed_description (failed_description_id,description) VALUES (1,'Variation maps to more than 3 different locations');
INSERT INTO failed_description (failed_description_id,description) VALUES (2,'None of the variant alleles match the reference allele');
INSERT INTO failed_description (failed_description_id,description) VALUES (3,'Variation has more than 3 different alleles');
INSERT INTO failed_description (failed_description_id,description) VALUES (4,'Loci with no observed variant alleles in dbSNP');
INSERT INTO failed_description (failed_description_id,description) VALUES (5,'Variation does not map to the genome'); 

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_54_55_c.sql|changed description for failed_description table');


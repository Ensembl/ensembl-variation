# add some failed descriptions

INSERT INTO failed_description (failed_description_id,description) VALUES (9,'Variation submission has been withdrawn by the 1000 genomes project due to high false positive rate');
INSERT INTO failed_description (failed_description_id,description) VALUES (11,'Allele string obtained from dbSNP does not agree with the submitted alleles'); 
INSERT INTO failed_description (failed_description_id,description) VALUES (12,'Variation has more than 3 different submitted alleles');         
INSERT INTO failed_description (failed_description_id,description) VALUES (13,'Alleles contain non-nucleotide characters');  
INSERT INTO failed_description (failed_description_id,description) VALUES (14,'Alleles contain ambiguity codes');  
INSERT INTO failed_description (failed_description_id,description) VALUES (15,'Mapped position is not compatible with reported alleles');

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch','patch_61_62_p.sql|add some failed descriptions');


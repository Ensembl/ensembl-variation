#this patch changes database schema version from 42 -> 43


#first, add a new column in the Variation table

ALTER Variation ADD failed_description_id int(10) unsigned not null;

#
# create the failed_description table
#
# Contains reasons for removing some variations from the Variation database
#
# failed_description_id  - primary key, internal identifier
# description - text containing the reason why the Variation information has been removed from the 
#               Variation databse except in the Variation table
#

CREATE TABLE failed_description(

 failed_description_id int(10) unsigned not null,
 description  text not null,

 PRIMARY KEY (failed_description_id)
);

#possible values in the failed_description table
INSERT INTO failed_description (failed_description_id,description) VALUES (1,'Variation has more than 3 different locations');
INSERT INTO failed_description (failed_description_id,description) VALUES (2,'Reference allele not present in the alleles of the variation');
INSERT INTO failed_description (failed_description_id,description) VALUES (3,'Variation containing more than 3 alleles');
INSERT INTO failed_description (failed_description_id,description) VALUES (4,'Variation with \'NoVariation\' alleles');

#and update the meta table with the patch applied
INSERT INTO meta (meta_key,meta_value) VALUES ('patch','patch_42_43_a.sql|update database schema');	

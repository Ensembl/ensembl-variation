#sql patch to find out variations that contain more than 3 alleles
#first, if there is no table to store the failed variations, create it
CREATE TABLE IF NOT EXISTS failed_variation(
    variation_id int(10) unsigned not null,
    failed_description_id int(10) unsigned not null,

    PRIMARY KEY(variation_id)
);

#then, copy the information from the variation_feature table
INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) 
   SELECT variation_id, 3
   FROM  variation_feature
   WHERE length(allele_string) - length(REPLACE(allele_string,'/','')) > 2
   AND allele_string not like '%-%';

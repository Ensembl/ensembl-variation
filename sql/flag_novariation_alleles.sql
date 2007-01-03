#sql patch to find out variations that contain alleles with "NoVariation" allele
#first, if there is no table to store the failed variations, create it
CREATE TABLE IF NOT EXISTS failed_variation(
    variation_id int(10) unsigned not null,
    failed_description_id int(10) unsigned not null,

    PRIMARY KEY(variation_id)
);

#then, copy the information from the Allele table
INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) 
   SELECT variation_id,4
   FROM  allele
   WHERE allele = 'NoVariation';

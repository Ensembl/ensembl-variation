#sql patch to find out variations that don't have variation_feature mappings
#first, if there is no table to store the failed variations, create it
CREATE TABLE IF NOT EXISTS failed_variation(
    variation_id int(10) unsigned not null,
    failed_description_id int(10) unsigned not null,

    PRIMARY KEY(variation_id)
);

#then, copy the information from the variation_feature table
INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) 
   SELECT v.variation_id, 5
   FROM  variation v LEFT JOIN variation_feature vf
   ON v.variation_id=vf.variation_id
   WHERE vf.variation_id is null;

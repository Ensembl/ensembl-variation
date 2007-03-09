####################
# change the schema of the database to include the failed_variation table
####################

CREATE TABLE IF NOT EXISTS failed_variation(
    variation_id int(10) unsigned not null,
    failed_description_id int(10) unsigned not null,

    PRIMARY KEY(variation_id)
);

#import the data in the failed_variation table

INSERT INTO failed_variation (variation_id, failed_description_id)
  SELECT variation_id, failed_description_id
  FROM variation
  WHERE failed_description_id <> 0;

#and remove the column from the variation table

ALTER TABLE variation DROP failed_description_id;

#and update the meta table with the patch applied
INSERT INTO meta (meta_key,meta_value) VALUES ('patch','patch_43_44_b.sql|update database schema');	

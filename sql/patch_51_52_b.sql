# patch_51_52_b.sql
#
# title: add seq_region table
#
# description:
# add toplevel seq_region table so people know which seq_region_id is corresponding to what chromosome name
# also to make sure seq_region_id in variation database is in sync with core db

CREATE TABLE seq_region (

  seq_region_id               INT(10) UNSIGNED NOT NULL,
  name                        VARCHAR(40) NOT NULL,

  PRIMARY KEY (seq_region_id),
  UNIQUE KEY name_idx (name)

) ;

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_51_52_b.sql|add seq_region table');

# added  a empty table subsnp_handle table
##################
CREATE TABLE `subsnp_handle` (
  `subsnp_id` int(11) unsigned NOT NULL,
  `handle` varchar(20) DEFAULT NULL,
  PRIMARY KEY (`subsnp_id`)
);
# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_55_56_c.sql|add subsnp_handle');

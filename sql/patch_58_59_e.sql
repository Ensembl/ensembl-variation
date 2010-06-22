# Add schema_type entry to meta table
ALTER TABLE `variation` CHANGE `validation_status` `validation_status` SET('cluster','freq','submitter','doublehit','hapmap','1000Genome','failed','precious');

ALTER TABLE `variation_feature` CHANGE `validation_status` `validation_status` SET('cluster','freq','submitter','doublehit','hapmap','1000Genome','failed','precious');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_58_59_e.sql|add 1000Genome validation_status');

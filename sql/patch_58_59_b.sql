# Add 'precious' validation_status to variation
ALTER TABLE `variation` CHANGE `validation_status` `validation_status` SET('cluster','freq','submitter','doublehit','hapmap','failed','precious');
# Add 'precious' validation_status to variation_feature
ALTER TABLE `variation_feature` CHANGE `validation_status` `validation_status` SET('cluster','freq','submitter','doublehit','hapmap','precious');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_58_59_b.sql|new validation status: precious');

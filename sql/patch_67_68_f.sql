# update the transcript_variation table
##################


ALTER TABLE allele ADD COLUMN  frequency_submitter_handle varchar(20);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_67_68_c.sql|update to store frequency submitter handle in allele table');

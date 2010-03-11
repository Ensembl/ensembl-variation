# update add two entries cds_start, cds_end in the transcript_variation table
##################
alter table transcript_variation add column cds_start int(11) after cdna_end;
alter table transcript_variation add column cds_end int(11) after cds_start;
ALTER TABLE transcript_variation CHANGE transcript_id transcript_stable_id VARCHAR(128) NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_57_58_b.sql|add cds_start, cds_end, change transcript_id to stable_id in transcript_variation');

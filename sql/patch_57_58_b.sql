# update add two entries cds_start, cds_end in the transcript_variation table
##################
alter table transcript_variation add column cds_start int(11) after cdna_end;
alter table transcript_variation add column cds_end int(11) after cds_start;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_57_58_b.sql|add cds_start, cds_end in transcript_variation');

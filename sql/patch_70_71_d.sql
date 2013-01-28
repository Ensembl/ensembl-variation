## drop unnecessary indices from consequence tables
alter table transcript_variation drop index feature_idx;
alter table motif_feature_variation drop index feature_idx;
alter table regulatory_feature_variation drop index feature_idx;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_d.sql|drop feature_idx from consequence tables');

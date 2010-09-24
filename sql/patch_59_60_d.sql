# drop and rebuild index on variation_feature and structural_variation
alter table variation_feature drop index pos_idx;
alter table structural_variation drop index pos_idx;

create index pos_idx on variation_feature (seq_region_id, seq_region_start, seq_region_end);
create index pos_idx on structural_variation (seq_region_id, seq_region_start, seq_region_end);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_59_60_d.sql|modify index on variation_feature and structural_variation');

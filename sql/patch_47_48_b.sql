###################################################
# add a new foreign key in variation and variaiton_synonym table, to allow get variations by source
###################################################
ALTER TABLE variation ADD KEY source_idx (source_id);

ALTER TABLE variation_synonym ADD KEY source_idx (source_id);

INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_47_48_b.sql|source_id key');

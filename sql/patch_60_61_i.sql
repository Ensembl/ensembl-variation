# add a failed description for variations that have no associated sequence

INSERT INTO failed_description (description) VALUES ('Variation has no associated sequence');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_i.sql|add a failed description for variations with no associated sequence');


# change the display name for insertions and deletions

UPDATE variation_class SET display_term = 'insertion' WHERE so_term = 'insertion';
UPDATE variation_class SET display_term = 'deletion' WHERE so_term = 'deletion';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_j.sql|change the display term for insertions and deletions');


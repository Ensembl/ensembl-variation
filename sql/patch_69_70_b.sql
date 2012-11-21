# drop table flanking_sequence
DROP TABLE IF EXISTS flanking_sequence;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_69_70_b.sql|drop table flanking_sequence');

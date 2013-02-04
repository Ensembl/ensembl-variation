# add ensembl evidence status

alter table variation add column evidence SET('Multiple_observations', 'Frequency','HapMap','1000Genomes', 'Cited');
alter table variation_feature add column evidence SET('Multiple_observations', 'Frequency','HapMap','1000Genomes', 'Cited');


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_b.sql|add evidence column');

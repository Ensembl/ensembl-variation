## Add ESP to evidence set in variation and variation_feature tables

ALTER TABLE variation CHANGE evidence evidence SET('Multiple_observations','Frequency','HapMap','1000Genomes','Cited', 'ESP');
ALTER TABLE variation_feature CHANGE evidence evidence SET('Multiple_observations','Frequency','HapMap','1000Genomes','Cited', 'ESP');

##patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_72_73_e.sql|Add ESP to varition set in variation and variation_feature tables');

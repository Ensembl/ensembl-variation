
## update evidence storage to use attribs rather than hard-coded values in the set for EnsemblGenomes 
## Note - data is not copied to new structure


ALTER TABLE variation ADD COLUMN evidence_attribs   SET('367','368','369','370','371','372') DEFAULT NULL;

ALTER TABLE variation_feature ADD COLUMN  evidence_attribs  SET('367','368','369','370','371','372') DEFAULT NULL;

ALTER TABLE variation DROP COLUMN  evidence;

ALTER TABLE variation_feature DROP COLUMN evidence;


#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_b.sql|Change evidence storage in Variation & Variation_feature table to attribs');
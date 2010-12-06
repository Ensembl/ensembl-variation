# change risk_allele_freq_in_controls column to double
ALTER TABLE `variation_annotation` CHANGE `risk_allele_freq_in_controls` `risk_allele_freq_in_controls` double NULL DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_60_61_d.sql|change type of column in variation_annotation to double');

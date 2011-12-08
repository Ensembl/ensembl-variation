# fix ploidy meta entry
UPDATE meta SET species_id = 1 WHERE meta_key = 'ploidy';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_65_66_a.sql|fix species_id in ploidy entry');


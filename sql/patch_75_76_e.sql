# add has_coverage flag to individual table
ALTER TABLE individual ADD has_coverage TINYINT(1) UNSIGNED NOT NULL DEFAULT 0 AFTER `display`;

#patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_e.sql|Add has_coverage flag to individual table');

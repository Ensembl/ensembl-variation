ALTER TABLE variation CHANGE ancestral_allele ancestral_allele varchar(255);

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch','patch_61_62_q.sql|revert ancestral_allele to varchar 255');


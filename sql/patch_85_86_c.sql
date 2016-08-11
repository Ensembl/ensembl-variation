-- index on study.external_reference

ALTER TABLE study add index external_reference_idx(external_reference);

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_85_86_c.sql|add index on study.external_reference');

# add a somatic_status column to source

ALTER TABLE source ADD COLUMN somatic_status ENUM ('germline','somatic','mixed') DEFAULT 'germline';

UPDATE source SET somatic_status = 'somatic' WHERE name = 'COSMIC';
UPDATE source SET somatic_status = 'mixed' WHERE name = 'dbSNP';

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch','patch_61_62_o.sql|add a somatic_status columns to source');


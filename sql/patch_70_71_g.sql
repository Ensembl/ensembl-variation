## add stable_id field to phenotype
ALTER TABLE phenotype ADD COLUMN `stable_id` VARCHAR(255) DEFAULT NULL AFTER phenotype_id;
CREATE INDEX stable_idx on phenotype(stable_id);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_g.sql|add stable_id column to phenotype');

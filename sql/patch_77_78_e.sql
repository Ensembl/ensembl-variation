# add indexes on father_ and mother_individual_id
CREATE INDEX father_individual_idx ON individual(father_individual_id);
CREATE INDEX mother_individual_idx ON individual(mother_individual_id);

# add index on population name
CREATE INDEX name_idx ON population(name);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_77_78_e.sql|add indexes on father_ and mother_individual_id and population name');

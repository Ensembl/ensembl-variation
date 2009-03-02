# ADD DISPLAY COLUMN TO SAMPLE

alter table sample
add column display enum('REFERENCE','DEFAULT','DISPLAYABLE','UNDISPLAYABLE') default 'UNDISPLAYABLE';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_53_54_c.sql|add display column to sample table');

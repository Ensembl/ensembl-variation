# ADD DESCRIPTION COLUMN TO SOURCE TABLE

alter table sample
add column display enum('REFERENCE','DEFAULT','DISPLAYABLE','UNDISPLAYABLE') default 'DISPLAYABLE';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_53_54_c.sql|schema version');

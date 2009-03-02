# ADD DESCRIPTION COLUMN TO SOURCE TABLE

alter table source
add column description varchar(255) default NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_53_54_b.sql|schema version');

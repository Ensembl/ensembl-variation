#Limit source.name column to be a maximum of 24 characters to be compatible with BioMart
##################

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_68_69_e.sql|length change for source.name ');


ALTER TABLE source MODIFY COLUMN name  varchar(24) NOT NULL;

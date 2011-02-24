# patch a MARTDISPLAYABLE option into the display field in sample
ALTER TABLE sample CHANGE display display enum('REFERENCE','DEFAULT','DISPLAYABLE','UNDISPLAYABLE','LD','MARTDISPLAYABLE') NULL DEFAULT 'UNDISPLAYABLE';


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_k.sql|add MARTDISPLAYABLE option to sample display column');
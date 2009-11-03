# update the sample table to include an 'LD' option for populations
ALTER TABLE sample MODIFY `display` ENUM('REFERENCE','DEFAULT','DISPLAYABLE','UNDISPLAYABLE','LD') DEFAULT 'UNDISPLAYABLE';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_56_57_b.sql|add LD option to display column in sample');

# Add coord_system_id column into the seq_region table
ALTER TABLE seq_region ADD COLUMN coord_system_id INT(10) UNSIGNED NOT NULL;
ALTER TABLE seq_region ADD KEY cs_idx (coord_system_id);
# Update unique key from name_idx to name_cs_idx ? 
# Check with Will
#ALTER TABLE seq_region DROP KEY name_idx;
#ALTER TABLE seq_region ADD UNIQUE KEY name_cs_idx(name, coord_system_id);


# Add coord_system table 

CREATE TABLE coord_system (

  coord_system_id             INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  species_id                  INT(10) UNSIGNED NOT NULL DEFAULT 1,
  name                        VARCHAR(40) NOT NULL,
  version                     VARCHAR(255) DEFAULT NULL,
  rank                        INT NOT NULL,
  attrib                      SET('default_version', 'sequence_level'),

  PRIMARY   KEY (coord_system_id),
  UNIQUE    KEY rank_idx (rank, species_id),
  UNIQUE    KEY name_idx (name, version, species_id),
            KEY species_idx (species_id)

);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_65_66_d.sql|Added coord_system table and update seq_region table in include coord_system_id in order to support multiple-species variation databases');

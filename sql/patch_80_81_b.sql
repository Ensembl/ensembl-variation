# create new table sample
CREATE TABLE IF NOT EXISTS sample (
  sample_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  individual_id int(10) unsigned NOT NULL,
  name varchar(255) DEFAULT NULL,
  description text,
  study_id int(10) unsigned DEFAULT NULL,
  display enum('REFERENCE','DEFAULT','DISPLAYABLE','UNDISPLAYABLE','LD','MARTDISPLAYABLE') DEFAULT 'UNDISPLAYABLE',
  has_coverage tinyint(1) unsigned NOT NULL DEFAULT '0',
  variation_set_id set('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64') DEFAULT NULL,
  PRIMARY KEY (sample_id),
  KEY individual_idx (individual_id),
  KEY study_idx (study_id)
);

RENAME TABLE individual TO individual_old;

# create new individual table
CREATE TABLE individual (
  individual_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  name varchar(255) DEFAULT NULL,
  description text,
  gender enum('Male','Female','Unknown') NOT NULL DEFAULT 'Unknown',
  father_individual_id int(10) unsigned DEFAULT NULL,
  mother_individual_id int(10) unsigned DEFAULT NULL,
  individual_type_id int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (individual_id),
  KEY father_individual_idx (father_individual_id),
  KEY mother_individual_idx (mother_individual_id)
);


INSERT INTO sample(sample_id, individual_id, name, description, display, has_coverage, variation_set_id) SELECT individual_id, individual_id, name, description, display, has_coverage, variation_set_id FROM individual_old;
INSERT INTO individual(individual_id, name, description, gender, father_individual_id, mother_individual_id, individual_type_id) SELECT individual_id, name, description, gender, father_individual_id, mother_individual_id, individual_type_id FROM individual_old;

DROP TABLE individual_old;


INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_80_81_b.sql|Create new sample table and update individual table. Copy individual data into new sample table.');

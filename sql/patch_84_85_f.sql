# add new phenotype_onology_accession table

CREATE TABLE `phenotype_ontology_accession` (
  `phenotype_id` int(11) unsigned NOT NULL,
  `accession` varchar(255) NOT NULL,
  `linked_by_attrib` set('437','438','439','440','441','442','443','444') DEFAULT NULL,
  PRIMARY KEY (`phenotype_id`,`accession`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_f.sql|add phenotype_ontology_accession');

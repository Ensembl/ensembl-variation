## create phenotype_feature and phenotype_feature_attrib tables
CREATE TABLE IF NOT EXISTS `phenotype_feature` (
  `phenotype_feature_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `phenotype_id` int(11) unsigned DEFAULT NULL,
  `source_id` int(11) unsigned DEFAULT NULL,
  `study_id` int(11) unsigned DEFAULT NULL,
  `type` enum('Gene','Variation','StructuralVariation','SupportingStructuralVariation','QTL','RegulatoryFeature') DEFAULT NULL,
  `object_id` varchar(255) DEFAULT NULL,
  `is_significant` tinyint(1) unsigned DEFAULT '1',
  `seq_region_id` int(11) unsigned DEFAULT NULL,
  `seq_region_start` int(11) unsigned DEFAULT NULL,
  `seq_region_end` int(11) unsigned DEFAULT NULL,
  `seq_region_strand` tinyint(4) DEFAULT NULL,
  PRIMARY KEY (`phenotype_feature_id`),
  KEY `phenotype_idx` (`phenotype_id`),
  KEY `object_idx` (`object_id`,`type`),
  KEY `type_idx` (`type`)
);

CREATE TABLE IF NOT EXISTS `phenotype_feature_attrib` (
  `phenotype_feature_id` int(11) unsigned NOT NULL,
  `attrib_type_id` int(11) DEFAULT NULL,
  `value` varchar(255) DEFAULT NULL,
  KEY `phenotype_feature_idx` (`phenotype_feature_id`)
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_f.sql|add phenotype_feature and phenotype_feature_attrib');

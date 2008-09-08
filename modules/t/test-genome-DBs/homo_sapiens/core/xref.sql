CREATE TABLE `xref` (
  `xref_id` int(10) unsigned NOT NULL auto_increment,
  `external_db_id` int(11) NOT NULL default '0',
  `dbprimary_acc` varchar(40) collate latin1_bin NOT NULL default '',
  `display_label` varchar(128) collate latin1_bin NOT NULL default '',
  `version` varchar(10) collate latin1_bin NOT NULL default '',
  `description` varchar(255) collate latin1_bin default NULL,
  `info_type` ENUM('PROJECTION', 'MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH','INFERRED_PAIR','PROBE','UNMAPPED'),
  `info_text` VARCHAR(255),
  `priority` INT,

  PRIMARY KEY  (`xref_id`),
  UNIQUE KEY `id_index` (`dbprimary_acc`,`external_db_id`),
  KEY `display_index` (`display_label`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;


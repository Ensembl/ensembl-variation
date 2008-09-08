-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `sample_synonym`
--

CREATE TABLE `sample_synonym` (
  `sample_synonym_id` int(11) NOT NULL auto_increment,
  `source_id` int(11) NOT NULL default '0',
  `name` varchar(255) default NULL,
  `sample_id` int(11) NOT NULL default '0',
  PRIMARY KEY  (`sample_synonym_id`),
  KEY `sample_idx` (`sample_id`),
  KEY `name` (`name`,`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




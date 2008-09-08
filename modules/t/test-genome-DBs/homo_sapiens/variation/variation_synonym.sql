-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `variation_synonym`
--

CREATE TABLE `variation_synonym` (
  `variation_synonym_id` int(11) NOT NULL auto_increment,
  `variation_id` int(11) NOT NULL default '0',
  `source_id` int(11) NOT NULL default '0',
  `name` varchar(255) default NULL,
  `moltype` varchar(50) default NULL,
  PRIMARY KEY  (`variation_synonym_id`),
  UNIQUE KEY `name` (`name`,`source_id`),
  KEY `variation_idx` (`variation_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




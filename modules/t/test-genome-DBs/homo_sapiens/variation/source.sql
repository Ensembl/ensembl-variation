-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `source`
--

CREATE TABLE `source` (
  `source_id` int(11) NOT NULL auto_increment,
  `name` varchar(255) default NULL,
  `version` int(11) default NULL,
  PRIMARY KEY  (`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




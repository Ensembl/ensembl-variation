-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `sample`
--

CREATE TABLE `sample` (
  `sample_id` int(11) NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `size` int(11) default NULL,
  `description` text,
  PRIMARY KEY  (`sample_id`),
  KEY `name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `httag`
--

CREATE TABLE `httag` (
  `httag_id` int(11) NOT NULL auto_increment,
  `variation_group_id` int(11) NOT NULL default '0',
  `name` varchar(255) default NULL,
  `source_id` int(11) NOT NULL default '0',
  PRIMARY KEY  (`httag_id`),
  KEY `variation_group_idx` (`variation_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




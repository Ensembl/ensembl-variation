-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `variation_group`
--

CREATE TABLE `variation_group` (
  `variation_group_id` int(11) NOT NULL auto_increment,
  `name` varchar(255) default NULL,
  `source_id` int(11) NOT NULL default '0',
  `type` enum('haplotype','tag') default NULL,
  PRIMARY KEY  (`variation_group_id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




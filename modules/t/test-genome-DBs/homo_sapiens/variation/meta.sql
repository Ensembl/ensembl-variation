-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `meta`
--

CREATE TABLE `meta` (
  `meta_id` int(11) NOT NULL auto_increment,
  `meta_key` varchar(40) NOT NULL default '',
  `meta_value` varchar(255) NOT NULL default '',
  PRIMARY KEY  (`meta_id`),
  KEY `meta_key_index` (`meta_key`),
  KEY `meta_value_index` (`meta_value`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




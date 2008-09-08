-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `variation_group_feature`
--

CREATE TABLE `variation_group_feature` (
  `variation_group_feature_id` int(11) NOT NULL auto_increment,
  `seq_region_id` int(11) NOT NULL default '0',
  `seq_region_start` int(11) NOT NULL default '0',
  `seq_region_end` int(11) NOT NULL default '0',
  `seq_region_strand` tinyint(4) NOT NULL default '0',
  `variation_group_id` int(11) NOT NULL default '0',
  `variation_group_name` varchar(255) default NULL,
  PRIMARY KEY  (`variation_group_feature_id`),
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`),
  KEY `variation_group_idx` (`variation_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




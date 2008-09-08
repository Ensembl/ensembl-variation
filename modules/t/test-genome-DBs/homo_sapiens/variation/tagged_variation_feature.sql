-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `tagged_variation_feature`
--

CREATE TABLE `tagged_variation_feature` (
  `variation_feature_id` int(11) NOT NULL default '0',
  `sample_id` int(11) NOT NULL default '0',
  PRIMARY KEY  (`variation_feature_id`,`sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `variation_group_variation`
--

CREATE TABLE `variation_group_variation` (
  `variation_id` int(11) NOT NULL default '0',
  `variation_group_id` int(11) NOT NULL default '0',
  UNIQUE KEY `variation_group_id` (`variation_group_id`,`variation_id`),
  KEY `variation_idx` (`variation_id`,`variation_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `population`
--

CREATE TABLE `population` (
  `is_strain` int(1) NOT NULL default '0',
  `sample_id` int(11) NOT NULL default '0',
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `allele_group_allele`
--

CREATE TABLE `allele_group_allele` (
  `allele_group_id` int(11) NOT NULL default '0',
  `allele` varchar(255) NOT NULL default '',
  `variation_id` int(11) NOT NULL default '0',
  UNIQUE KEY `allele_group_id` (`allele_group_id`,`variation_id`),
  KEY `allele_idx` (`variation_id`,`allele_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




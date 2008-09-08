-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `individual_genotype_multiple_bp`
--

CREATE TABLE `individual_genotype_multiple_bp` (
  `variation_id` int(11) NOT NULL default '0',
  `allele_1` varchar(255) default NULL,
  `allele_2` varchar(255) default NULL,
  `sample_id` int(11) NOT NULL default '0',
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




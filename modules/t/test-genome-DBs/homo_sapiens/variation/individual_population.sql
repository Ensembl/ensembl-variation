-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `individual_population`
--

CREATE TABLE `individual_population` (
  `population_sample_id` int(11) NOT NULL default '0',
  `individual_sample_id` int(11) NOT NULL default '0',
  KEY `individual_sample_idx` (`individual_sample_id`),
  KEY `population_sample_idx` (`population_sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `population_structure`
--

CREATE TABLE `population_structure` (
  `super_population_sample_id` int(11) NOT NULL default '0',
  `sub_population_sample_id` int(11) NOT NULL default '0',
  UNIQUE KEY `super_population_sample_id` (`super_population_sample_id`,`sub_population_sample_id`),
  KEY `sub_pop_sample_idx` (`sub_population_sample_id`,`super_population_sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `population_genotype`
--

CREATE TABLE `population_genotype` (
  `population_genotype_id` int(11) NOT NULL auto_increment,
  `variation_id` int(11) NOT NULL default '0',
  `allele_1` varchar(255) default NULL,
  `allele_2` varchar(255) default NULL,
  `frequency` float default NULL,
  `sample_id` int(11) NOT NULL default '0',
  PRIMARY KEY  (`population_genotype_id`),
  KEY `variation_idx` (`variation_id`),
  KEY `sample_idx` (`sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




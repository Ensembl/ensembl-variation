-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `allele`
--

CREATE TABLE `allele` (
  `allele_id` int(11) NOT NULL auto_increment,
  `variation_id` int(11) NOT NULL default '0',
  `allele` text,
  `frequency` float default NULL,
  `sample_id` int(11) default NULL,
  PRIMARY KEY  (`allele_id`),
  KEY `variation_idx` (`variation_id`),
  KEY `allele_idx` (`allele_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




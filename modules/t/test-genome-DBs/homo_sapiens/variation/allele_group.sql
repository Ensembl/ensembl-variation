-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `allele_group`
--

CREATE TABLE `allele_group` (
  `allele_group_id` int(11) NOT NULL auto_increment,
  `variation_group_id` int(11) NOT NULL default '0',
  `name` varchar(255) default NULL,
  `source_id` int(11) default NULL,
  `frequency` float default NULL,
  `sample_id` int(11) default NULL,
  PRIMARY KEY  (`allele_group_id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




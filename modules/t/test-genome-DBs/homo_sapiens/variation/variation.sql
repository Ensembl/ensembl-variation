-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `variation`
--

CREATE TABLE `variation` (
  `variation_id` int(11) NOT NULL auto_increment,
  `source_id` int(11) NOT NULL default '0',
  `name` varchar(255) default NULL,
  `validation_status` set('cluster','freq','submitter','doublehit','hapmap') default NULL,
  `ancestral_allele` text,
  PRIMARY KEY  (`variation_id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




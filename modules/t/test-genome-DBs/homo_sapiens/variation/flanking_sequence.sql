-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `flanking_sequence`
--

CREATE TABLE `flanking_sequence` (
  `variation_id` int(11) NOT NULL default '0',
  `up_seq` text,
  `down_seq` text,
  `up_seq_region_start` int(11) default NULL,
  `up_seq_region_end` int(11) default NULL,
  `down_seq_region_start` int(11) default NULL,
  `down_seq_region_end` int(11) default NULL,
  `seq_region_id` int(11) default NULL,
  `seq_region_strand` tinyint(4) default NULL,
  PRIMARY KEY  (`variation_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;




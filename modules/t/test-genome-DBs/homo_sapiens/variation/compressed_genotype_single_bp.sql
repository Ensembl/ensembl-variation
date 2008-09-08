-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `compressed_genotype_single_bp`
--

CREATE TABLE `compressed_genotype_single_bp` (
  `sample_id` int(11) NOT NULL default '0',
  `seq_region_id` int(11) NOT NULL default '0',
  `seq_region_start` int(11) NOT NULL default '0',
  `seq_region_end` int(11) NOT NULL default '0',
  `seq_region_strand` tinyint(4) NOT NULL default '0',
  `genotypes` blob,
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




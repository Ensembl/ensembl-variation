-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `variation_feature`
--

CREATE TABLE `variation_feature` (
  `variation_feature_id` int(11) NOT NULL auto_increment,
  `seq_region_id` int(11) NOT NULL default '0',
  `seq_region_start` int(11) NOT NULL default '0',
  `seq_region_end` int(11) NOT NULL default '0',
  `seq_region_strand` tinyint(4) NOT NULL default '0',
  `variation_id` int(11) NOT NULL default '0',
  `allele_string` text,
  `variation_name` varchar(255) default NULL,
  `map_weight` int(11) NOT NULL default '0',
  `flags` set('genotyped') default NULL,
  `source_id` int(1) NOT NULL default '1',
  `validation_status` set('cluster','freq','submitter','doublehit','hapmap') default NULL,
  `consequence_type` enum('INTRONIC','UPSTREAM','DOWNSTREAM','SYNONYMOUS_CODING','NON_SYNONYMOUS_CODING','FRAMESHIFT_CODING','5PRIME_UTR','3PRIME_UTR','INTERGENIC') NOT NULL default 'INTERGENIC',
  PRIMARY KEY  (`variation_feature_id`),
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`),
  KEY `variation_idx` (`variation_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




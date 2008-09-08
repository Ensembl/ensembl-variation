-- MySQL dump 10.9
--
-- Host: ecs2    Database: _test_db_homo_sapiens_variation_dr2_25_11_12432
-- ------------------------------------------------------
-- Server version	4.1.12-log


--
-- Table structure for table `transcript_variation`
--

CREATE TABLE `transcript_variation` (
  `transcript_variation_id` int(11) NOT NULL auto_increment,
  `transcript_id` int(11) NOT NULL default '0',
  `variation_feature_id` int(11) NOT NULL default '0',
  `cdna_start` int(11) default NULL,
  `cdna_end` int(11) default NULL,
  `translation_start` int(11) default NULL,
  `translation_end` int(11) default NULL,
  `peptide_allele_string` varchar(255) default NULL,
  `consequence_type` enum('INTRONIC','UPSTREAM','DOWNSTREAM','SYNONYMOUS_CODING','NON_SYNONYMOUS_CODING','FRAMESHIFT_CODING','5PRIME_UTR','3PRIME_UTR') NOT NULL default 'INTRONIC',
  PRIMARY KEY  (`transcript_variation_id`),
  KEY `variation_idx` (`variation_feature_id`),
  KEY `transcript_idx` (`transcript_id`),
  KEY `consequence_type_idx` (`consequence_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




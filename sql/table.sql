-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.


/**
Use MyISAM storage engine
*/
SET default_storage_engine=MYISAM;

/**
@header  Variation tables
@desc    These tables define the central variation data.
@colour  #000
*/


/**
@table variation

@desc This is the schema's generic representation of a variation, defined as a genetic feature that varies between individuals of the same species. 
      The most common type is the single nucleotide variation (SNP) though the schema also accommodates copy number variations (CNVs) and structural variations (SVs).<br />
			In Ensembl, a variation is defined by its flanking sequence rather than its mapped location on a chromosome; a variation may in fact have multiple mappings across a genome, 
			although this fails our <a href="/info/genome/variation/data_description.html#quality_control">Quality Control</a>.<br /> 
      This table stores a variation's name (commonly an ID of the form e.g. rs123456, assigned by dbSNP), along with an ancestral (or reference) allele.

@column variation_id		       Primary key, internal identifier.
@column source_id			         Foreign key references to the @link source table.
@column name				           Name of the variation. e.g. "rs1333049".
@column ancestral_allele	     Ancestral allele for the variant.
@column flipped				         This is set to 1 if the variant is flipped from the negative to the positive strand during import.
@column class_attrib_id		     Class of the variation, key into the @link attrib table.<br /> The list of variation classes is available <a href="/info/genome/variation/data_description.html#classes">here</a>.
@column somatic                Flags whether this variation is known to be somatic or not
@column minor_allele           The minor allele of this variant. The minor allele is the second most frequent allele.
@column minor_allele_freq      The 'global' frequency of the minor allele of this variant, as reported by dbSNP. The minor allele frequency is the frequency of the second most frequent allele.
@column minor_allele_count     The number of samples the minor allele of this variant is found in. The minor allele is the second most frequent allele.
@column clinical_significance  A set of clinical significance classes assigned to the variant.<br /> 
                               The list of clinical significances is available <a href="/info/genome/variation/data_description.html#clin_significance">here</a>.
@column evidence_attribs            A summary of the evidence supporting a variant as a guide to its potential reliability. See the evidence descriptions <a href="/info/genome/variation/data_description.html#evidence_status">here</a>.
@column display                Flags whether this variation should be displayed in browser tracks and returned by default by the API

@see variation_synonym
@see failed_variation
@see variation_feature
@see allele
@see sample_genotype_multiple_bp
@see compressed_genotype_var
@see attrib
*/

create table variation (
  variation_id int(10) unsigned not null auto_increment, # PK
  source_id int(10) unsigned not null, 
  name varchar(255),
  ancestral_allele varchar(255) DEFAULT NULL,
  flipped tinyint(1) unsigned NULL DEFAULT NULL,
  class_attrib_id int(10) unsigned default 0,
  somatic tinyint(1) DEFAULT 0 NOT NULL,
  minor_allele varchar(50) DEFAULT NULL,
  minor_allele_freq float DEFAULT NULL,
  minor_allele_count int(10) unsigned DEFAULT NULL,
  clinical_significance SET('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective'),
  evidence_attribs   SET('367','368','369','370','371','372','418','421') DEFAULT NULL,
  display int(1) DEFAULT 1,

	primary key( variation_id ),
	unique ( name ),
	key source_idx (source_id)
);


/**
@table variation_attrib

@desc This table stores miscellaneous attributes associated with a variation entry.

@column variation_id     Foreign key references @link variation table
@column attrib_id			   Foreign key references @link attrib table, describes the attribute
@column value            Attribute value

@see variation
@see attrib
*/

CREATE TABLE variation_attrib (
  variation_id INT(11) UNSIGNED NOT NULL,
  attrib_id INT(11) DEFAULT NULL,
  value VARCHAR(255) DEFAULT NULL,
  KEY variation_idx (variation_id),
  KEY attrib_value_idx (attrib_id,value)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table variation_feature

@desc This table represents mappings of variations to genomic locations. It stores an allele string representing the different possible alleles that are found at that locus e.g. "A/T" for a SNP, as well as a "worst case" consequence of the mutation. It also acts as part of the relationship between variations and transcripts.

@column variation_feature_id	 Primary key, internal identifier.
@column seq_region_id			     Foreign key references @link seq_region in core db. Refers to the seq_region which this variant is on, which may be a chromosome, a clone, etc...
@column seq_region_start		   The start position of the variation on the @link seq_region.
@column seq_region_end			   The end position of the variation on the @link seq_region.
@column seq_region_strand		   The orientation of the variation on the @link seq_region.
@column variation_id				   Foreign key references to the @link variation table.
@column allele_string			     This is a denormalised string taken from the alleles in the allele table associated with this variation. The reference allele (i.e. one on the reference genome comes first).
@column variation_name			   A denormalisation taken from the variation table. This is the name or identifier that is used for displaying the feature.
@column map_weight				     The number of times that this variation has mapped to the genome. This is a denormalisation as this particular feature is one example of a mapped location. This can be used to limit the the features that come back from a query.
@column flags						       Flag to filter the selection of variations.
@column source_id					     Foreign key references to the source table.
@column consequence_types		   The SO term(s) of all unique observed consequence types of this variation feature.<br /> The list of consequence descriptions is available <a href="/info/genome/variation/predicted_data.html#consequences">here</a>.
@column variation_set_id		   The variation feature can belong to a @link variation_set.
@column class_attrib_id			   Class of the variation, key in the @link attrib table.<br /> The list of variation classes is available <a href="/info/genome/variation/data_description.html#classes">here</a>.
@column somatic                Flags whether this variation_feature is somatic or germline
@column minor_allele           The minor allele of this variant. The minor allele is the second most frequent allele.
@column minor_allele_freq      The 'global' frequency of the minor allele of this variant, as reported by dbSNP. The minor allele frequency is the frequency of the second most frequent allele.
@column minor_allele_count     The number of samples the minor allele of this variant is found in. The minor allele is the second most frequent allele.
@column alignment_quality      Quality of alignment for variants mapped by flanks rather than position justified.
@column evidence_attribs       A summary of the evidence supporting a variant as a guide to its potential reliability. See the evidence descriptions <a href="/info/genome/variation/data_description.html#evidence_status">here</a>.
@column clinical_significance  A set of clinical significance classes assigned to the variant.<br /> 
                               The list of clinical significances is available <a href="/info/genome/variation/data_description.html#clin_significance">here</a>.
@column display                Flags whether this variation should be displayed in browser tracks and returned by default by the API

@see variation
@see transcript_variation
@see seq_region
@see attrib
*/

create table variation_feature(
  variation_feature_id int(10) unsigned not null auto_increment,
  seq_region_id int(10) unsigned not null,
  seq_region_start int not null,
  seq_region_end int not null,
  seq_region_strand tinyint not null,
  variation_id int(10) unsigned not null,
  allele_string varchar(50000),
  variation_name varchar(255),
  map_weight int not null,
  flags SET('genotyped'),
  source_id int(10) unsigned not null,
  consequence_types SET (
        'intergenic_variant',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'stop_lost',
        'coding_sequence_variant',
        'missense_variant',
        'stop_gained',
        'synonymous_variant',
        'frameshift_variant',
        'non_coding_transcript_variant',
        'non_coding_transcript_exon_variant',
        'mature_miRNA_variant',
        'NMD_transcript_variant',
        '5_prime_UTR_variant',
        '3_prime_UTR_variant',
        'incomplete_terminal_codon_variant',
        'intron_variant',
        'splice_region_variant',
        'downstream_gene_variant',
        'upstream_gene_variant',
        'start_lost',
        'stop_retained_variant',
        'inframe_insertion',
        'inframe_deletion',
        'transcript_ablation',
        'transcript_fusion',
        'transcript_amplification',
        'transcript_translocation',
        'TFBS_ablation',
        'TFBS_fusion',
        'TFBS_amplification',
        'TFBS_translocation',
        'regulatory_region_ablation',
        'regulatory_region_fusion',
        'regulatory_region_amplification',
        'regulatory_region_translocation',
        'feature_elongation',
        'feature_truncation',
        'regulatory_region_variant',
        'TF_binding_site_variant',
        'protein_altering_variant',
        'start_retained_variant'
    ) DEFAULT 'intergenic_variant' NOT NULL,
    variation_set_id SET (
            '1','2','3','4','5','6','7','8',
            '9','10','11','12','13','14','15','16',
            '17','18','19','20','21','22','23','24',
            '25','26','27','28','29','30','31','32',
            '33','34','35','36','37','38','39','40',
            '41','42','43','44','45','46','47','48',
            '49','50','51','52','53','54','55','56',
            '57','58','59','60','61','62','63','64'
    ) NOT NULL DEFAULT '',
    class_attrib_id int(10) unsigned default 0,
    somatic tinyint(1) DEFAULT 0 NOT NULL,
    minor_allele varchar(50) DEFAULT NULL,
    minor_allele_freq float DEFAULT NULL,
    minor_allele_count int(10) unsigned DEFAULT NULL,
    alignment_quality double  DEFAULT NULL,
    evidence_attribs   SET('367','368','369','370','371','372','418','421') DEFAULT NULL,    
    clinical_significance SET('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective') DEFAULT NULL,
    display int(1) DEFAULT 1,

   	primary key( variation_feature_id ),
	  key pos_idx( seq_region_id, seq_region_start, seq_region_end ),
	  key variation_idx( variation_id ),
    key variation_set_idx ( variation_set_id ),
    key consequence_type_idx (consequence_types),
    key source_idx (source_id)
);


/**
@table variation_synonym

@desc This table allows for a variation to have multiple IDs, generally given by multiple sources.

@column variation_synonym_id	Primary key, internal identifier.
@column variation_id					Foreign key references to the variation table.
@column subsnp_id							Foreign key references to the subsnp_handle table.
@column source_id							Foreign key references to the source table.
@column name									Name of the synonym variation. e.g. 'rs1333049'.

@see source
@see variation
@see subsnp_handle
*/

create table variation_synonym (
  variation_synonym_id int(10) unsigned not null auto_increment,
  variation_id int(10) unsigned not null,
  subsnp_id int(15) unsigned ,
  source_id int(10) unsigned not null,
  name varchar(255),

  primary key(variation_synonym_id),
  key variation_idx (variation_id),
  key subsnp_idx(subsnp_id),
  unique (name, source_id),
  key source_idx (source_id)
);


/**
@table allele

@desc This table stores information about each of a variation's alleles, along with population frequencies.

@column allele_id		   Primary key, internal identifier.
@column variation_id	 	   Foreign key references to the @link variation table.
@column subsnp_id		   Foreign key references to the @link subsnp_handle table.
@column allele_code_id 		   Foreign key reference to @link allele_code table.
@column population_id		   Foreign key references to the @link population table.
@column frequency		   Frequency of this allele in the population.
@column count			   Number of individuals/samples in the population where this allele is found.
@column frequency_submitter_handle dbSNP handle for submitter of frequency data [may be different to submitter of observed variant]

@see variation
@see subsnp_handle
@see allele_code
@see population
@see submitter_handle
*/

CREATE TABLE allele (
  allele_id int(11) NOT NULL AUTO_INCREMENT,
  variation_id int(11) unsigned NOT NULL,
  subsnp_id int(11) unsigned DEFAULT NULL,
  allele_code_id int(11) unsigned NOT NULL,
  population_id int(11) unsigned DEFAULT NULL,
  frequency float unsigned DEFAULT NULL,
  count int(11) unsigned DEFAULT NULL,
  frequency_submitter_handle int(10) DEFAULT NULL,

  PRIMARY KEY (allele_id),
  KEY variation_idx (variation_id),
  KEY subsnp_idx (subsnp_id),
  KEY population_idx (population_id)
);




/**
@header  Phenotype tables
@desc    These tables store information linking entities (variants, genes, QTLs) with phenotypes and other annotations.
@colour  #22949b
*/


/**
@table phenotype_feature

@desc This table stores information linking entities (variants, genes, QTLs) and phenotypes.

@column phenotype_feature_id	  Primary key, internal identifier.
@column phenotype_id			  Foreign key references to the @link phenotype table.
@column source_id				  Foreign key references to the @link source table.
@column study_id				  Foreign key references to the @link study table.
@column type					  Type of object associated.
@column object_id	              Stable identifier for associated object.
@column is_significant			  Flag indicating if the association is statistically significant in the given study.
@column seq_region_id			  Foreign key references @link seq_region in core db. Refers to the seq_region which this feature is on, which may be a chromosome, a clone, etc...
@column seq_region_start		  The start position of the feature on the @link seq_region.
@column seq_region_end			  The end position of the feature on the @link seq_region.
@column seq_region_strand		  The orientation of the feature on the @link seq_region.

@see variation
@see phenotype
@see source
@see study
*/

CREATE TABLE IF NOT EXISTS `phenotype_feature` (
  `phenotype_feature_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `phenotype_id` int(11) unsigned DEFAULT NULL,
  `source_id` int(11) unsigned DEFAULT NULL,
  `study_id` int(11) unsigned DEFAULT NULL,
  `type` enum('Gene','Variation','StructuralVariation','SupportingStructuralVariation','QTL','RegulatoryFeature') DEFAULT NULL,
  `object_id` varchar(255) DEFAULT NULL,
  `is_significant` tinyint(1) unsigned DEFAULT '1',
  `seq_region_id` int(11) unsigned DEFAULT NULL,
  `seq_region_start` int(11) unsigned DEFAULT NULL,
  `seq_region_end` int(11) unsigned DEFAULT NULL,
  `seq_region_strand` tinyint(4) DEFAULT NULL,
  PRIMARY KEY (`phenotype_feature_id`),
  KEY `phenotype_idx` (`phenotype_id`),
  KEY `object_idx` (`object_id`,`type`),
  KEY `type_idx` (`type`),
  KEY `pos_idx` (`seq_region_id`,`seq_region_start`,`seq_region_end`),
  KEY `source_idx` (`source_id`)
);


/**
@table phenotype_feature_attrib

@desc This table stores additional information on a given phenotype/object association. It is styled as an attrib table to allow for a variety of fields to be populated across different object types.

@column phenotype_feature_id	  Foreign key, references to the @link phenotype_feature table.
@column attrib_type_id			  Foreign key references to the @link attrib_type table.
@column value	    			  The value of the attribute.

@see phenotype_feature
@see attrib_type
*/

CREATE TABLE IF NOT EXISTS `phenotype_feature_attrib` (
  `phenotype_feature_id` int(11) unsigned NOT NULL,
  `attrib_type_id` int(11) DEFAULT NULL,
  `value` varchar(255) DEFAULT NULL,
  KEY `phenotype_feature_idx` (`phenotype_feature_id`),
  KEY `type_value_idx` (`attrib_type_id`,`value`)
);


/**
@table phenotype

@desc This table stores details of the phenotypes associated with phenotype_features.

@column phenotype_id         Primary key, internal identifier.
@column stable_id            Ensembl stable identifier for the phenotype
@column name                 Phenotype short name. e.g. "CAD".
@column description varchar  Phenotype long name. e.g. "Coronary Artery Disease".

@see phenotype_feature
*/

CREATE TABLE `phenotype` (
  `phenotype_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `stable_id` varchar(255) DEFAULT NULL,
  `name` varchar(50) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`phenotype_id`),
  KEY `name_idx` (`name`),
  UNIQUE KEY `desc_idx` (`description`),
  KEY `stable_idx` (`stable_id`)
);

/**
@table phenotype_ontology_accession

@desc This table stores accessions of phenotype ontology terms which have been linked to phenotype.descriptions

@column phenotype_id         Foreign key, references to the @link phenotype table.
@column accession            The accession of an ontology term held in the ontology database (eg. EFO:0000378) 
@column mapped_by_attrib     The method used to annotate the phenotype.description with the ontology term
@column mapping_type         The relation defining the association between the ontology term and the phenotype.description

@see phenotype
*/
CREATE TABLE `phenotype_ontology_accession` (
  `phenotype_id` int(11) unsigned NOT NULL,
  `accession` varchar(255) NOT NULL,
  `mapped_by_attrib` set('437','438','439','440','441','442','443','444') DEFAULT NULL,
  `mapping_type` enum('is','involves') DEFAULT NULL,
  PRIMARY KEY (`phenotype_id`,`accession`),
  KEY `accession_idx` (`accession`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@header  Other tables
@desc    These tables define the other data associated with a variation.
@colour  #000
*/


/**
@table subsnp_handle

@desc This table contains the SubSNP(ss) ID and the name of the submitter handle of dbSNP.

@column subsnp_id	Primary key. It corresponds to the subsnp identifier (ssID) from dbSNP.<br />This ssID is stored in this table without the "ss" prefix. e.g. "120258606" instead of "ss120258606".
@column handle		The name of the dbSNP handler who submitted the ssID.<br />Name of the synonym (a different <b>sample_id</b>).

@see allele
@see failed_variation
@see population_genotype
@see sample
@see variation_synonym
*/

create table subsnp_handle (
  subsnp_id int(11) unsigned not null,
  handle varchar(20),

  primary key(subsnp_id)
);


/**
@table submitter_handle

@desc This table holds a short string to distinguish data submitters

@column handle_id	Primary key, internal identifier.
@column handle      	Short string assigned to the data submitter. 

@see allele 
*/

create table submitter_handle (
  handle_id int(10) unsigned not null auto_increment,
  handle varchar(25),
 primary key( handle_id ),
        unique ( handle )
);


/**
@table allele_code

@desc This table stores the relationship between the internal allele identifiers and the alleles themselves.

@column allele_code_id	Primary key, internal identifier.
@column allele      	String representing the allele. Has a unique constraint on the first 1000 characters (max allowed by MySQL).

@example See below the first entries in the allele_code table:
         @sql SELECT * FROM allele_code LIMIT 6;

@see allele
@see genotype_code
*/

CREATE TABLE allele_code (
  allele_code_id int(11) NOT NULL AUTO_INCREMENT,
  allele varchar(60000) DEFAULT NULL,
  
  PRIMARY KEY (allele_code_id),
  UNIQUE KEY allele_idx (allele(1000))
);


/**
@table genotype_code

@desc This table stores genotype codes as multiple rows of allele_code identifiers, linked by genotype_code_id and ordered by haplotype_id.

@column genotype_code_id	Internal identifier.
@column allele_code_id 	    Foreign key reference to @link allele_code table.
@column haplotype_id        Sorting order of the genotype's alleles.
@column phased              Indicates if this genotype is phased

@see allele_code
@see population_genotype
*/

CREATE TABLE genotype_code (
  genotype_code_id int(11) unsigned NOT NULL,
  allele_code_id int(11) unsigned NOT NULL,
  haplotype_id tinyint(2) unsigned NOT NULL,
  phased tinyint(2) unsigned DEFAULT NULL,
  
  KEY genotype_code_id (genotype_code_id),
  KEY allele_code_id (allele_code_id)
);


/**
@table seq_region

@desc This table stores the relationship between Ensembl's internal coordinate system identifiers and traditional chromosome names.

@column seq_region_id	   Primary key. Foreign key references seq_region in core db. Refers to the seq_region which this variant is on, which may be a chromosome, a clone, etc...
@column name				     The name of this sequence region.
@column coord_system_id  Foreign key references to the @link coord_system table.

@see variation_feature
@see compressed_genotype_region
@see read_coverage
*/

CREATE TABLE seq_region (

  seq_region_id               INT(10) UNSIGNED NOT NULL,
  name                        VARCHAR(255) NOT NULL,
  coord_system_id             INT(10) UNSIGNED NOT NULL,

  PRIMARY KEY (seq_region_id),
  # Which one, check with Will
  #UNIQUE KEY name_idx (name),
  UNIQUE KEY name_cs_idx (name, coord_system_id),
  KEY cs_idx (coord_system_id)

) ;


/**
@table coord_system
@desc Stores information about the available co-ordinate systems for the species identified through the species_id field.
Note that for each species, there must be one co-ordinate system that has the attribute "top_level" and one that has the attribute "sequence_level".

@column coord_system_id      Primary key, internal identifier.
@column species_id           Identifies the species for multi-species databases.
@column name                 Co-oridinate system name, e.g. 'chromosome', 'contig', 'scaffold' etc.
@column version              Assembly.
@column rank                 Co-oridinate system rank.
@column attrib               Co-oridinate system attrib (e.g. "top_level", "sequence_level").

@see seq_region
@see meta_coord
@see meta

*/

CREATE TABLE coord_system (

  coord_system_id             INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  species_id                  INT(10) UNSIGNED NOT NULL DEFAULT 1,
  name                        VARCHAR(40) NOT NULL,
  version                     VARCHAR(255) DEFAULT NULL,
  rank                        INT NOT NULL,
  attrib                      SET('default_version', 'sequence_level'),

  PRIMARY   KEY (coord_system_id),
  UNIQUE    KEY rank_idx (rank, species_id),
  UNIQUE    KEY name_idx (name, version, species_id),
            KEY species_idx (species_id)

);



/**
@header  Sample tables
@desc    These tables define the sample, individual and population information.
@colour  #FF8500
*/


/**
@table population

@colour #FF8500
@desc Stores information about a population. A population may be an ethnic group (e.g. Caucasian, Hispanic), assay group (e.g. 24 Europeans), phenotypic group (e.g. blue eyed, diabetes) etc. Populations may be composed of other populations by defining relationships in the population_structure table.

@column population_id     Primary key, internal identifier.
@column name              Name of the population.
@column size              Size of the population.
@column description       Description of the population.
@column collection        Flag indicating if the population is defined based on geography (0) or a collection of individuals/samples with respect to some other criteria (1).
@column freqs_from_gts    Flag indicating if the population frequencies can be retrieved from the allele table (0) or from the individual/sample genotypes (1).
@column display           Information used by BioMart.
@column display_group_id  Used to group population for display on the Population Genetics page

@see population_synonym
@see sample_population
@see population_structure
@see population_genotype
@see allele
@see display_group
*/

CREATE TABLE population(
    population_id int(10) unsigned not null auto_increment,
    name varchar(255),
    size int(10),
    description text,
    collection tinyint(1) default 0,
    freqs_from_gts tinyint(1),
    display enum('LD', 'MARTDISPLAYABLE', 'UNDISPLAYABLE') default 'UNDISPLAYABLE',
    display_group_id tinyint(1) ,

    primary key(population_id),
    key name_idx (name)
);


/**
@table population_structure

@colour #FF8500
@desc This table stores hierarchical relationships between populations by relating them as populations and sub-populations.

@column super_population_id    Foreign key references to the population table.
@column sub_population_id      Foreign key references to the population table.

@see population
*/

create table population_structure (
  super_population_id int(10) unsigned not null,
  sub_population_id int(10) unsigned not null,

  unique super_population_idx (super_population_id, sub_population_id),
  key sub_population_idx (sub_population_id)
);

/**
@table individual

@colour #FF8500
@desc Stores information about an identifiable individual, including gender and the identifiers of the individual's parents (if known).

@column individual_id           Primary key, internal identifier.
@column name                    Name of the individual.
@column description             Description of the individual.
@column gender                  The sex of this individual.
@column father_individual_id    Self referential ID, the father of this individual if known.
@column mother_individual_id    Self referential ID, the mother of this individual if known.
@column individual_type_id      Foreign key references to the @link individual_type table.

@see individual_synonym
@see individual_type
@see sample
*/

create table individual(
  individual_id int(10) unsigned not null auto_increment,
  name varchar(255),
  description text,
  gender enum('Male', 'Female', 'Unknown') default 'Unknown' NOT NULL,
  father_individual_id int(10) unsigned,
  mother_individual_id int(10) unsigned,
  individual_type_id int(10) unsigned NOT NULL DEFAULT 0,
  primary key(individual_id),
  key father_individual_idx (father_individual_id),
  key mother_individual_idx (mother_individual_id)
);

/**
@table individual_type

@colour #FF8500
@desc This table gives a detailed description for each of the possible individual types: fully_inbred, partly_inbred, outbred, mutant  
@column individual_type_id	Primary key, internal identifier.
@column name				Short name of the individual type. e.g. "fully_inbred","mutant".
@column description			Long name of the individual type.

@example See below the list of individual types:
         @sql SELECT * FROM individual_type;

@see individual
*/

create table individual_type(
  individual_type_id int(0) unsigned not null auto_increment,
  name varchar(255) not null,
  description text,
  
  primary key (individual_type_id)
);

#this table will always contain the same values

INSERT INTO individual_type (name,description) VALUES ('fully_inbred','multiple organisms have the same genome sequence');
INSERT INTO individual_type (name,description) VALUES ('partly_inbred','single organisms have reduced genome variability due to human intervention');
INSERT INTO individual_type (name,description) VALUES ('outbred','a single organism which breeds freely');
INSERT INTO individual_type (name,description) VALUES ('mutant','a single or multiple organisms with the same genome sequence that have a natural or experimentally induced mutation');

/**
@table sample

@colour #FF8500
@desc Stores information about a sample. A sample belongs to an individual. An individual can have multiple samples. A sample can belong only to one individual. A sample can be associated
with a study. 

@column sample_id               Primary key, internal identifier.
@column individual_id           Foreign key references to the @link individual table.
@column name                    Name of the sample.
@column description             Description of the sample.
@column study_id                Foreign key references to the @link study table.
@column display                 Information used by the website: samples with little information are filtered from some web displays.
@column has_coverage            Indicate if the sample has coverage data populated in the read coverage table
@column variation_set_id        Indicates the variation sets for which a sample has genotypes

@see individual
@see study
@see sample_population
@see sample_synonym
@see sample_genotype_multiple_bp
@see read_coverage
@see compressed_genotype_region
@see compressed_genotype_var
@see variation_set
*/

CREATE TABLE sample(
  sample_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  individual_id int(10) unsigned NOT NULL,
  name varchar(255) DEFAULT NULL,
  description text,
  study_id int(10) unsigned DEFAULT NULL,
  display enum('REFERENCE','DEFAULT','DISPLAYABLE','UNDISPLAYABLE','LD','MARTDISPLAYABLE') DEFAULT 'UNDISPLAYABLE',
  has_coverage tinyint(1) unsigned NOT NULL DEFAULT '0',
  variation_set_id set('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64') DEFAULT NULL,
  PRIMARY KEY (sample_id),
  KEY individual_idx (individual_id),
  KEY study_idx (study_id)
);

/**
@table sample_synonym

@colour #FF8500
@desc Used to store alternative names for samples when data comes from multiple sources.

@column synonym_id      Primary key, internal identifier.
@column sample_id       Foreign key references to the @link sample table.
@column source_id       Foreign key references to the @link source table.
@column name            Name of the synonym.

@see sample
@see source
*/

CREATE TABLE sample_synonym (
  synonym_id int(10) unsigned not null auto_increment,
  sample_id int(10) unsigned not null,
  source_id int(10) unsigned not null,
  name varchar(255),
  primary key(synonym_id),
  key sample_idx (sample_id),
  key (name, source_id)
);
  
/**
@table sample_population

@colour #FF8500
@desc This table resolves the many-to-many relationship between the sample and population tables; i.e. samples may belong to more than one population. Hence it is composed of rows of sample and population identifiers.

@column sample_id	Foreign key references to the @link sample table.
@column population_id	Foreign key references to the @link population table.

@see sample
@see population
*/

create table sample_population (
  sample_id int(10) unsigned not null,
  population_id int(10) unsigned not null,

  key sample_idx(sample_id),
  key population_idx(population_id)

);

/**
@table individual_synonym

@colour #FF8500
@desc Used to store alternative names for individuals when data comes from multiple sources.

@column synonym_id       Primary key, internal identifier.
@column individual_id    Foreign key references to the @link individual table.
@column source_id        Foreign key references to the @link source table.
@column name             Name of the synonym.

@see individual
@see source
*/

CREATE TABLE individual_synonym (
  synonym_id int(10) unsigned not null auto_increment,
  individual_id int(10) unsigned not null,
  source_id int(10) unsigned not null,
  name varchar(255),

  primary key(synonym_id),
  key individual_idx (individual_id),
  key (name, source_id)
);

/**
@table population_synonym

@colour #FF8500
@desc Used to store alternative names for populations when data comes from multiple sources.

@column synonym_id       Primary key, internal identifier.
@column population_id    Foreign key references to the @link population table.
@column source_id        Foreign key references to the @link source table.
@column name             Name of the synonym.

@see population
@see source
*/

CREATE TABLE population_synonym (
  synonym_id int(10) unsigned not null auto_increment,
  population_id int(10) unsigned not null,
  source_id int(10) unsigned not null,
  name varchar(255),

  primary key(synonym_id),
  key population_idx (population_id),
  key (name, source_id)
);

/**
@table display_group

@colour #FF8500
@desc Used to store groups of populations displayed separately on the Population Genetics page

@column display_group_id     Primary key, internal identifier.
@column display_priority     Priority level for group (smallest number is highest on page) 
@column display_name         Name of the group to be displayed as the table header.

@see population
*/
CREATE TABLE display_group(
  display_group_id int(10) unsigned not null auto_increment ,
  display_priority int(10) unsigned not null, 
  display_name     varchar(255) not null,

	primary key( display_group_id ),
	unique ( display_name ),
	unique ( display_priority )
 );


/**
@header  Genotype tables
@desc    These tables define the genotype data at the sample and population levels.
@colour  #FF8500
*/


/**
@table population_genotype

@colour #FF8500
@desc This table stores genotypes and frequencies for variations in given populations.

@column population_genotype_id    Primary key, internal identifier.
@column variation_id              Foreign key references to the @link variation table.
@column subsnp_id                 Foreign key references to the subsnp_handle table.
@column genotype_code_id          Foreign key reference to the @link genotype_code table.
@column frequency                 Frequency of the genotype in the population.
@column population_id             Foreign key references to the @link population table.
@column count                     Number of individuals/samples who have this genotype, in this population.

@see population
@see variation
@see subsnp_handle
@see genotype_code
*/


CREATE TABLE population_genotype (
  population_genotype_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  variation_id int(11) unsigned NOT NULL,
  subsnp_id int(11) unsigned DEFAULT NULL,
  genotype_code_id int(11) DEFAULT NULL,
  frequency float DEFAULT NULL,
  population_id int(10) unsigned DEFAULT NULL,
  count int(10) unsigned DEFAULT NULL,
  
  PRIMARY KEY (population_genotype_id),
  KEY population_idx (population_id),
  KEY variation_idx (variation_id),
  KEY subsnp_idx (subsnp_id)
);


/**
@table tmp_sample_genotype_single_bp

@colour #FF8500
@desc This table is only needed to create master schema when run healthcheck system. Needed for other species, but human, so keep it.

@column variation_id     Foreign key references to the @link variation table.
@column subsnp_id        Foreign key references to the @link subsnp_handle table.
@column allele_1         One of the alleles of the genotype, e.g. "TAG".
@column allele_2         The other allele of the genotype.
@column sample_id        Foreign key references to the @link sample table.

@see sample
@see variation
@see subsnp_handle
*/

CREATE TABLE tmp_sample_genotype_single_bp (
	variation_id int(10) not null,
	subsnp_id int(15) unsigned,   
	allele_1 char(1),
	allele_2 char(1),
	sample_id int(10) unsigned not null,

	key variation_idx(variation_id),
    key subsnp_idx(subsnp_id),
    key sample_idx(sample_id)
) MAX_ROWS = 100000000;


/**
@table sample_genotype_multiple_bp

@colour #FF8500
@desc This table holds uncompressed genotypes for given variations.

@column variation_id     Primary key. Foreign key references to the @link variation table.
@column subsnp_id        Foreign key references to the @link subsnp_handle table.
@column allele_1         One of the alleles of the genotype, e.g. "TAG".
@column allele_2         The other allele of the genotype.
@column sample_id        Foreign key references to the @link sample table.

@see sample
@see variation
@see subsnp_handle
*/

create table sample_genotype_multiple_bp (
  variation_id int(10) unsigned not null,
  subsnp_id int(15) unsigned,	
  allele_1 varchar(25000),
  allele_2 varchar(25000),
  sample_id int(10) unsigned,

  key variation_idx(variation_id),
  key subsnp_idx(subsnp_id),
  key sample_idx(sample_id)
);


/**
@table compressed_genotype_region

@colour #FF8500
@desc This table holds genotypes compressed using the pack() method in Perl. These genotypes are mapped to particular genomic locations rather than variation objects. The data have been compressed to reduce table size and increase the speed of the web code when retrieving strain slices and LD data. Only data from resequenced samples are used for LD calculations are included in this table

@column sample_id            Foreign key references to the @link sample table.
@column seq_region_id        Foreign key references @link seq_region in core db. ers to the seq_region which this variant is on, which may be a chromosome, a clone, etc...
@column seq_region_start     The start position of the variation on the @link seq_region.
@column seq_region_end       The end position of the variation on the @link seq_region.
@column seq_region_strand    The orientation of the variation on the @link seq_region.
@column genotypes            Encoded representation of the genotype data:<br />Each row in the compressed table stores genotypes from one individual/sample in one fixed-size region of the genome (arbitrarily defined as 100 Kb). The compressed string (using Perl's pack method) consisting of a repeating triplet of elements: a  <span style="color:#D00">distance</span> in base pairs from the previous genotype; a <span style="color:#090">variation dbID</span>; a <span style="color:#00D">genotype_code_id</span> identifier.<br />For example, a given row may have a start position of 1000, indicating the chromosomal position of the first genotype in this row. The unpacked genotypes field then may contain the following elements:<br /><b><span style="color:#D00">0</span>, <span style="color:#090">1</span>,  <span style="color:#00D">1</span>, <span style="color:#D00">20</span>, <span style="color:#090">2</span>, <span style="color:#00D">5</span>, <span style="color:#D00">35</span>, <span style="color:#090">3</span>, <span style="color:#00D">3</span>, ...</b><br />The first genotype ("<span style="color:#D00">0</span>,<span style="color:#090">1</span>,<span style="color:#00D">1</span>") has a position of 1000 + <span style="color:#D00">0</span> = 1000, and corresponds to the variation with the internal identifier <span style="color:#090">1</span> and genotype_code_id corresponding to the genotype A|G (internal ID <span style="color:#00D">1</span>).<br />The second genotype ("<span style="color:#D00">20</span>,<span style="color:#090">2</span>,<span style="color:#00D">5</span>") has a position of 1000 + <span style="color:#D00">20</span> = 1020, internal variation_id <span style="color:#090">2</span> and genotype_code_id corresponding to the genotype C|C ( internal ID <span style="color:#00D">5</span>).<br />The third genotype similarly has a position of 1055, and so on.

@see sample
@see seq_region
@see variation
@see genotype_code
*/

CREATE TABLE compressed_genotype_region (
  sample_id int(10) unsigned NOT NULL,
  seq_region_id int(10) unsigned NOT NULL,
  seq_region_start int(11) NOT NULL,
  seq_region_end int(11) NOT NULL,
  seq_region_strand tinyint(4) NOT NULL,
  genotypes blob,
  
  KEY pos_idx (seq_region_id,seq_region_start),
  KEY sample_idx (sample_id)
);

/**
@table compressed_genotype_var

@colour #FF8500
@desc This table holds genotypes compressed using the pack() method in Perl. These genotypes are mapped directly to variation objects. The data have been compressed to reduce table size. All genotypes in the database are included in this table (included duplicates of those genotypes contained in the compressed_genotype_region table). This table is optimised for retrieval from variation.

@column variation_id	Foreign key references to the @link variation table.
@column subsnp_id		  Foreign key references to the @link subsnp_handle table.
@column genotypes     Encoded representation of the genotype data:<br />Each row in the compressed table stores genotypes from one subsnp of a variation (or one variation if no subsnp is defined). The compressed string (using Perl's pack method) consisting of a repeating pair of elements: an internal sample_id corresponding to a sample; a genotype_code_id identifier.

@see sample
@see variation
@see genotype_code
*/

CREATE TABLE compressed_genotype_var (
  variation_id int(11) unsigned NOT NULL,
  subsnp_id int(11) unsigned DEFAULT NULL,
  genotypes blob,
  
  KEY variation_idx (variation_id),
  KEY subsnp_idx (subsnp_id)
);


/**
@table read_coverage

@colour #FF8500
@desc This table stores the read coverage of resequenced samples. Each row contains sample ID, chromosomal coordinates and a read coverage level.

@column seq_region_id       Foreign key references @link seq_region in core db. ers to the seq_region which this variant is on, which may be a chromosome, a clone, etc...
@column seq_region_start    The start position of the variation on the @link seq_region.
@column seq_region_end      The end position of the variation on the @link seq_region.
@column level               Minimum number of reads.
@column sample_id           Foreign key references to the @link sample table.

@see sample
@see seq_region
*/

CREATE TABLE read_coverage (
  seq_region_id int(10) unsigned not null,
  seq_region_start int not null,
  seq_region_end int not null,
  level tinyint not null,
  sample_id int(10) unsigned not null,
  KEY seq_region_idx(seq_region_id,seq_region_start),
  KEY sample_idx (sample_id)
);


/**
@header  Structural variation tables
@desc    These tables define the structural variation data.
@colour  #01D4F7
*/


/**
@table structural_variation

@colour #01D4F7
@desc This table stores information about structural variation.

@column structural_variation_id	Primary key, internal identifier.
@column variation_name					The external identifier or name of the variation. e.g. "esv9549".
@column alias                   Other structural variation name.
@column source_id								Foreign key references to the @link source table.
@column study_id								Foreign key references to the @link study table.	
@column class_attrib_id					Foreign key references to the @link attrib table. Defines the type of structural variant.<br /> 
                                The list of structural variation classes is available <a href="/info/genome/variation/data_description.html#classes">here</a>.
@column clinical_significance   A set of clinical significance classes assigned to the structural variant.<br /> 
                                The list of clinical significances is available <a href="/info/genome/variation/data_description.html#clin_significance">here</a>.
@column validation_status				Validation status of the variant.
@column is_evidence             Flag indicating if the structural variation is a supporting evidence (1) or not (0).
@column somatic                 Flags whether this structural variation is known to be somatic or not
@column copy_number             Add the copy number for the CNV supporting structural variants when available.

@see source
@see study
@see attrib
*/

CREATE TABLE structural_variation (
  structural_variation_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  variation_name varchar(255) DEFAULT NULL,
  alias varchar(255) DEFAULT NULL,
	source_id int(10) unsigned NOT NULL,
  study_id int(10) unsigned DEFAULT NULL,
	class_attrib_id int(10) unsigned NOT NULL DEFAULT 0,
	clinical_significance SET('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective') DEFAULT NULL,
  validation_status ENUM('validated','not validated','high quality'),
	is_evidence TINYINT(4) DEFAULT 0,
	somatic TINYINT(1) NOT NULL DEFAULT 0,
	copy_number TINYINT(2) DEFAULT NULL,
	
  PRIMARY KEY (structural_variation_id),
  UNIQUE (variation_name),
	KEY source_idx (source_id),
	KEY study_idx (study_id),
	KEY attrib_idx (class_attrib_id)
);


/**
@table structural_variation_association

@colour #01D4F7
@desc This table stores the associations between structural variations and their supporting evidences.

@column structural_variation_id	            Primary key. Foreign key references to the @link structural_variation table.
@column supporting_structural_variation_id	Primary key. Foreign key references to the @link structural_variation table.

@see structural_variation
*/

CREATE TABLE structural_variation_association (
  structural_variation_id int(10) unsigned NOT NULL,
  supporting_structural_variation_id int(10) unsigned NOT NULL,
	
  PRIMARY KEY (structural_variation_id, supporting_structural_variation_id),
	KEY structural_variation_idx (structural_variation_id),
	KEY supporting_structural_variation_idx (supporting_structural_variation_id)
);


/**
@table structural_variation_feature

@colour #01D4F7
@desc This table stores information about structural variation features (i.e. mappings of structural variations to genomic locations).

@column structural_variation_feature_id	 Primary key, internal identifier.
@column seq_region_id						         Foreign key references @link seq_region in core db. Refers to the seq_region which this variant is on, which may be a chromosome, a clone, etc...
@column outer_start				               The 5' outer bound position of the variation on the @link seq_region.
@column seq_region_start			 	         The start position of the variation on the @link seq_region.
@column inner_start				               The 5' inner bound position of the variation on the @link seq_region.
@column inner_end   			               The 3' inner bound position of the variation on the @link seq_region.
@column seq_region_end					         The end position of the variation on the @link seq_region.
@column outer_end				                 The 3' outer bound position of the variation on the @link seq_region.
@column seq_region_strand				         The orientation of the variation on the @link seq_region.
@column structural_variation_id	         Foreign key references to the @link structural_variation table.
@column variation_name					         A denormalisation taken from the structural_variation table. This is the name or identifier that is used for displaying the feature (e.g. "esv9549").
@column source_id								         Foreign key references to the @link source table.
@column study_id								         Foreign key references to the @link study table
@column class_attrib_id					         Foreign key references to the @link attrib table. Defines the type of structural variant.<br /> 
                                         The list of structural variation classes is available <a href="/info/genome/variation/data_description.html#classes">here</a>.
@column allele_string						         The variant allele, where known.
@column is_evidence                      Flag indicating if the structural variation is a supporting evidence (1) or not (0).
@column variation_set_id		             The structural variation feature can belong to a @link variation_set.
@column somatic                          Flags whether this structural variation is known to be somatic or not
@column breakpoint_order                 Defines the order of the breakpoints when several events/mutation occurred for a structural variation (e.g. somatic mutations)
@column length                           Length of the structural variant. Used for the variants with a class "insertion", when the size of the insertion is known.

@see structural_variation
@see source
@see study
@see seq_region
@see attrib
@see variation_set
*/

create table structural_variation_feature (
	structural_variation_feature_id int(10) unsigned NOT NULL AUTO_INCREMENT,
	seq_region_id int(10) unsigned NOT NULL,
	outer_start int,	
	seq_region_start int NOT NULL,
	inner_start int,
	inner_end int,
	seq_region_end int NOT NULL,
	outer_end int,
	seq_region_strand tinyint NOT NULL,
	structural_variation_id int(10) unsigned NOT NULL,
  variation_name varchar(255),
	source_id int(10) unsigned NOT NULL,
  study_id int(10) unsigned DEFAULT NULL,
  class_attrib_id int(10) unsigned NOT NULL DEFAULT 0,
	allele_string longtext DEFAULT NULL,
	is_evidence tinyint(1) NOT NULL DEFAULT 0,
	somatic TINYINT(1) NOT NULL DEFAULT 0,
	breakpoint_order TINYINT(4) DEFAULT NULL,
	length int(10) DEFAULT NULL,
  variation_set_id SET (
          '1','2','3','4','5','6','7','8',
          '9','10','11','12','13','14','15','16',
          '17','18','19','20','21','22','23','24',
          '25','26','27','28','29','30','31','32',
          '33','34','35','36','37','38','39','40',
          '41','42','43','44','45','46','47','48',
          '49','50','51','52','53','54','55','56',
          '57','58','59','60','61','62','63','64'
  ) NOT NULL DEFAULT '',
	
  PRIMARY KEY (structural_variation_feature_id),
	KEY pos_idx( seq_region_id, seq_region_start, seq_region_end ),
	KEY structural_variation_idx (structural_variation_id),
	KEY source_idx (source_id),
	KEY study_idx (study_id),
	KEY attrib_idx (class_attrib_id),
	KEY variation_set_idx (variation_set_id)
);


/**
@table structural_variation_sample

@colour #01D4F7
@desc This table stores sample and strain information for structural variants and their supporting evidences.

@column structural_variation_sample_id  Primary key, internal identifier.
@column structural_variation_id         Foreign key references to the @link structural_variation table.
@column sample_id		                    Foreign key references to the @link sample table. Defines the individual or sample name.
@column zygosity                        Define the numeric zygosity of the structural variant for the sample, when available.

@see structural_variation
@see sample
@see individudal
*/

CREATE TABLE structural_variation_sample (
	structural_variation_sample_id int(10) unsigned NOT NULL auto_increment,
	structural_variation_id int(10) unsigned NOT NULL,
	sample_id int(10) unsigned DEFAULT NULL,
	zygosity TINYINT(1) DEFAULT NULL,
	
	primary key (structural_variation_sample_id),
	key structural_variation_idx(structural_variation_id),
	key sample_idx(sample_id)
);

/**
@header  Variation set tables
@desc    These tables define the variation and structural variation set data. The complete list of variation sets with their descriptions is available <a href="/info/genome/variation/data_description.html#variation_sets">here</a>.
@colour  #FFD700
*/


/**
@table variation_set_variation

@colour #FFD700
@desc A table for mapping variations to variation_sets.

@column variation_id			Primary key. Foreign key references to the @link variation table.
@column variation_set_id	Primary key. Foreign key references to the @link variation_set table.

@see variation
@see variation_set
*/

CREATE TABLE IF NOT EXISTS variation_set_variation (
	variation_id int(10) unsigned NOT NULL,
	variation_set_id int(10) unsigned NOT NULL,
	PRIMARY KEY (variation_id,variation_set_id),
	KEY variation_set_idx (variation_set_id,variation_id)
);

/**
@table variation_set

@colour #FFD700
@desc This table contains the name of sets and subsets of variations stored in the database. It usually represents the name of the project or subproject where a group of variations has been identified.

@column variation_set_id			Primary key, internal identifier.
@column name									Name of the set e.g. "Phenotype-associated variations".
@column description						Description of the set.
@column short_name_attrib_id	Foreign key references to the @link attrib table. Short name used for web purpose.

@example See below the command to display the list of variation set entries, e.g. for human:
         @sql SELECT * FROM variation_set;

@see variation_set_variation
@see variation_set_structure
*/
 
CREATE TABLE IF NOT EXISTS variation_set (
	variation_set_id int(10) unsigned NOT NULL AUTO_INCREMENT,
	name VARCHAR(255),
	description TEXT,
	short_name_attrib_id INT(10) UNSIGNED DEFAULT NULL,
	PRIMARY KEY (variation_set_id),
	KEY name_idx (name)
);


/**
@table variation_set_structure

@colour #FFD700
@desc This table stores hierarchical relationships between variation sets by relating them as variation sets and variation subsets.

@column variation_set_super	Primary key. Foreign key references to the @link variation_set table.
@column variation_set_sub		Primary key. Foreign key references to the @link variation_set table.

@see variation_set
*/

CREATE TABLE IF NOT EXISTS variation_set_structure (
	variation_set_super int(10) unsigned NOT NULL,
	variation_set_sub int(10) unsigned NOT NULL,
	PRIMARY KEY (variation_set_super,variation_set_sub),
	KEY sub_idx (variation_set_sub,variation_set_super)
);


/**
@table variation_set_structural_variation

@colour #FFD700
@desc A table for mapping structural variations to variation_sets.

@column structural_variation_id  Primary key. Foreign key references to the @link structural_variation table.
@column variation_set_id	       Primary key. Foreign key references to the @link variation_set table.

@see structural_variation
@see variation_set
*/

CREATE TABLE IF NOT EXISTS variation_set_structural_variation (
	structural_variation_id int(10) unsigned NOT NULL,
	variation_set_id int(10) unsigned NOT NULL,
	PRIMARY KEY (structural_variation_id,variation_set_id)
);




/**
@header  Variation effect tables
@desc    These tables define the variation effect prediction data in different Ensembl features.
@colour  #FF4DC8
*/


/**
@table transcript_variation

@colour #FF4DC8
@desc This table relates a single allele of a variation_feature to a transcript (see Core documentation). It contains the consequence of the allele e.g. intron_variant, non_synonymous_codon, stop_lost etc, along with the change in amino acid in the resulting protein if applicable.

@column transcript_variation_id	 Primary key, internal identifier.
@column feature_stable_id		     Foreign key to core databases. Unique stable id of related transcript.
@column variation_feature_id		 Foreign key references to the @link variation_feature table.
@column allele_string            Shows the reference sequence and variant sequence of this allele
@column somatic                  Flags if the associated variation is known to be somatic
@column consequence_types			   The consequence(s) of the variant allele on this transcript.<br /> The list of consequence descriptions is available <a href="/info/genome/variation/predicted_data.html#consequences">here</a>. 
@column cds_start					       The start position of variation in cds coordinates.
@column cds_end						       The end position of variation in cds coordinates.
@column cdna_start					     The start position of variation in cdna coordinates.
@column cdna_end					       The end position of variation in cdna coordinates.
@column translation_start			   The start position of variation on peptide.
@column translation_end				   The end position of variation on peptide.
@column distance_to_transcript   Only for upstream or downstream variants, it gives the distance from the start or the end of the transcript
@column codon_allele_string      The reference and variant codons
@column pep_allele_string        The reference and variant peptides
@column hgvs_genomic             HGVS representation of this allele with respect to the genomic sequence
@column hgvs_transcript          HGVS representation of this allele with respect to the [coding or non-coding] transcript
@column hgvs_protein             HGVS representation of this allele with respect to the protein
@column polyphen_prediction      The PolyPhen prediction for the effect of this allele on the protein
@column polyphen_score           The PolyPhen score corresponding to the prediction 
@column sift_prediction          The SIFT prediction for the effect of this allele on the protein 
@column sift_score               The SIFT score corresponding to this prediction
@column display                  Flags whether this transcript_variation should be displayed in browser tracks and returned by default by the API

@see variation_feature
*/

CREATE TABLE transcript_variation (
    transcript_variation_id             int(11) unsigned NOT NULL AUTO_INCREMENT,
    variation_feature_id                int(11) unsigned NOT NULL,
    feature_stable_id                   varchar(128) DEFAULT NULL,
    allele_string                       text,
    somatic                             tinyint(1) NOT NULL DEFAULT 0,
    consequence_types                   set (
                                            'splice_acceptor_variant',
                                            'splice_donor_variant',
                                            'stop_lost',
                                            'coding_sequence_variant',
                                            'missense_variant',
                                            'stop_gained',
                                            'synonymous_variant',
                                            'frameshift_variant',
                                            'non_coding_transcript_variant',
                                            'non_coding_transcript_exon_variant',
                                            'mature_miRNA_variant',
                                            'NMD_transcript_variant',
                                            '5_prime_UTR_variant',
                                            '3_prime_UTR_variant',
                                            'incomplete_terminal_codon_variant',
                                            'intron_variant',
                                            'splice_region_variant',
                                            'downstream_gene_variant',
                                            'upstream_gene_variant',
                                            'start_lost',
                                            'stop_retained_variant',
                                            'inframe_insertion',
                                            'inframe_deletion',
                                            'transcript_ablation',
                                            'transcript_fusion',
                                            'transcript_amplification',
                                            'transcript_translocation',
                                            'feature_elongation',
                                            'feature_truncation',
                                            'protein_altering_variant',
                                            'start_retained_variant'
                                        ),
    cds_start                           int(11) unsigned,
    cds_end                             int(11) unsigned,
    cdna_start                          int(11) unsigned,
    cdna_end                            int(11) unsigned,
    translation_start                   int(11) unsigned,
    translation_end                     int(11) unsigned,
    distance_to_transcript              int(11) unsigned,
    codon_allele_string                 text,
    pep_allele_string                   text,
    hgvs_genomic                        text,
    hgvs_transcript                     text,
    hgvs_protein                        text,
    polyphen_prediction                 enum('unknown', 'benign', 'possibly damaging', 'probably damaging') DEFAULT NULL,
    polyphen_score                      float DEFAULT NULL,
    sift_prediction                     enum('tolerated', 'deleterious', 'tolerated - low confidence', 'deleterious - low confidence') DEFAULT NULL,
    sift_score                          float DEFAULT NULL,
    display                             int(1) DEFAULT 1,

    PRIMARY KEY                         (transcript_variation_id),
    KEY variation_feature_idx           (variation_feature_id),
    KEY consequence_type_idx            (consequence_types),
    KEY somatic_feature_idx             (feature_stable_id, somatic)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table variation_hgvs

@colour #FF4DC8
@desc This table is used in web index creation. It links a variation_id to all possible transcript and protein level change descriptions in HGVS annotation.
@column variation_id         Primary key, foreign key references @link variation
@column hgvs_name            Primary key, HGVS change description
*/

create table variation_hgvs(
variation_id int(10) unsigned not null,
hgvs_name varchar(255) not null,
primary key(variation_id, hgvs_name));


/**
@table variation_genename

@colour #FF4DC8
@desc This table is used in web index creation. It links a variation_id to the names of the genes the variation is within
@column variation_id         Primary key, foreign key references @link variation
@column gene_name            Primary key, display name of gene
*/

create table variation_genename (
variation_id int(10) unsigned not null, 
gene_name varchar(255) not null, 
primary key(variation_id, gene_name));


/**
@table motif_feature_variation

@colour #FF4DC8
@desc This table relates a single allele of a variation_feature to a motif feature (see Regulation documentation). It contains the consequence of the allele.

@column motif_feature_variation_id  Primary key, internal identifier.
@column variation_feature_id        Foreign key references to the @link variation_feature table.
@column feature_stable_id		        Foreign key to regulation databases. Unique stable id of related regulatory_feature.
@column motif_feature_id            Foreign key to regulation databases. Internal id of related motif_feature.
@column allele_string               Shows the reference sequence and variant sequence of this allele.
@column somatic                     Flags if the associated variation is known to be somatic.
@column consequence_types		        The consequence(s) of the variant allele on this motif_feature.<br /> The list of consequence descriptions is available <a href="/info/genome/variation/predicted_data.html#consequences">here</a>.
@column motif_name                  The display label of the motif.
@column motif_start                 The start position of the variation in the motif.
@column motif_end                   The end position of the variation in the motif.
@column motif_score_delta           The deviation from the score (that is derived from alignment software (e.g. MOODS)) caused by the variation.
@column in_informative_position     Flags if the variation is in an informative position.

@see variation_feature
*/

CREATE TABLE IF NOT EXISTS motif_feature_variation (
    motif_feature_variation_id          int(11) unsigned NOT NULL AUTO_INCREMENT,
    variation_feature_id                int(11) unsigned NOT NULL,
    feature_stable_id                   varchar(128) DEFAULT NULL,
    motif_feature_id                    int(11) unsigned NOT NULL,
    allele_string                       text,
    somatic                             tinyint(1) NOT NULL DEFAULT 0,
    consequence_types                   set (
                                          'TF_binding_site_variant',
                                          'TFBS_ablation',
                                          'TFBS_fusion',
                                          'TFBS_amplification',
                                          'TFBS_translocation'
                                        ),
    motif_name                          varchar(60) DEFAULT NULL,
    motif_start                         int(11) unsigned,
    motif_end                           int(11) unsigned,
    motif_score_delta                   float DEFAULT NULL,
    in_informative_position             tinyint(1) NOT NULL DEFAULT 0,

    PRIMARY KEY                         (motif_feature_variation_id),
    KEY variation_feature_idx           (variation_feature_id),
    KEY feature_stable_idx              (feature_stable_id),
    KEY consequence_type_idx            (consequence_types),
    KEY somatic_feature_idx             (feature_stable_id, somatic)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table regulatory_feature_variation

@colour #FF4DC8
@desc This table relates a single allele of a variation_feature to a regulatory feature (see Regulation documentation). It contains the consequence of the allele.

@column regulatory_feature_variation_id  Primary key, internal identifier.
@column variation_feature_id             Foreign key references to the @link variation_feature table.
@column feature_stable_id		             Foreign key to regulation databases. Unique stable id of related regulatory_feature.
@column feature_type                     The name of the feature type.
@column allele_string                    Shows the reference sequence and variant sequence of this allele.
@column somatic                          Flags if the associated variation is known to be somatic.
@column consequence_types		            The consequence(s) of the variant allele on this regulatory feature.<br /> The list of consequence descriptions is available <a href="/info/genome/variation/predicted_data.html#consequences">here</a>.

@see variation_feature
*/

CREATE TABLE IF NOT EXISTS regulatory_feature_variation (
    regulatory_feature_variation_id     int(11) unsigned NOT NULL AUTO_INCREMENT,
    variation_feature_id                int(11) unsigned NOT NULL,
    feature_stable_id                   varchar(128) DEFAULT NULL,
    feature_type                        text, 
    allele_string                       text,
    somatic                             tinyint(1) NOT NULL DEFAULT 0,
    consequence_types                   set (
                                          'regulatory_region_variant',
                                          'regulatory_region_ablation',
                                          'regulatory_region_fusion',
                                          'regulatory_region_amplification',
                                          'regulatory_region_translocation'
                                        ),

    PRIMARY KEY                         (regulatory_feature_variation_id),
    KEY variation_feature_idx           (variation_feature_id),
    KEY feature_stable_idx              (feature_stable_id),
    KEY consequence_type_idx            (consequence_types),
    KEY somatic_feature_idx             (feature_stable_id, somatic)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
	
	

/**
@header  Source/study tables
@desc    These tables define the variation source and study information.
@colour  #7CFC00
*/


/**
@table source

@colour #7CFC00
@desc This table contains details of the source from which a variation is derived. Most commonly this is NCBI's dbSNP; other sources include SNPs called by Ensembl.<br />
You can see the complete list, by species, <a href="/info/genome/variation/sources_documentation.html">here</a>.

@column source_id		    Primary key, internal identifier.
@column name			      Name of the source. e.g. "dbSNP"
@column version		      Version number of the source (if available). e.g. "132"
@column description	    Description of the source.
@column url				      URL of the source.
@column type			      Define the type of the source, e.g. 'chip'
@column somatic_status  Indicates if this source includes somatic or germline mutations, or a mixture
@column data_types      Indicates the type(s) of data provided by the source
 
@example See below the command listing all the data sources in the human variation database:
         @sql SELECT * FROM source ORDER BY source_id;

@see variation
@see variation_synonym
@see variation_feature
@see individual_synonym
@see phenotype_feature
@see population_synonym
@see sample_synonym
@see structural_variation
@see structural_variation_feature
@see study
*/

create table source(
	source_id int(10) unsigned not null auto_increment,
	name varchar(24) not null,
	version int,
	description varchar(400),
	url varchar(255),
	type ENUM('chip','lsdb') DEFAULT NULL,
  somatic_status ENUM ('germline','somatic','mixed') DEFAULT 'germline',
  data_types SET('variation','variation_synonym','structural_variation','phenotype_feature','study') DEFAULT NULL,
	
	primary key( source_id )
);


/**
@table study

@colour #7CFC00
@desc This table contains details of the studies.
			The studies information can come from internal studies (DGVa, EGA) or from external studies (UniProt, NHGRI, ...).

@column study_id						Primary key, internal identifier.
@column source_id						Foreign key references to the @link source table.
@column name								Name of the study. e.g. "EGAS00000000001"
@column description					Description of the study.
@column url									URL to find the study data (http or ftp).
@column external_reference	The PubMed/id or project name associated with this study.
@column study_type					Displays the type of the study (e.g. genome-wide association study, control-set, case-set, curated, ...).

@see source
@see phenotype_feature
@see structural_variation
@see structural_variation_feature
*/

CREATE TABLE study (
	study_id int(10) unsigned not null auto_increment,
	source_id int(10) unsigned not null,
	name varchar(255) DEFAULT null,
	description text DEFAULT NULL,
	url varchar(255) DEFAULT NULL,
	external_reference varchar(255) DEFAULT NULL,
	study_type varchar(255) DEFAULT NULL,
	
	primary key( study_id ),
	key source_idx (source_id),
  key external_reference_idx(external_reference)
);


/**
@table associate_study

@colour #7CFC00
@desc This table contains identifiers of associated studies (e.g. NHGRI and EGA studies with the same PubMed identifier).

@column study1_id		Primary key. Foreign key references to the @link study table.
@column study2_id		Primary key. Foreign key references to the @link study table.

@see study
*/
CREATE TABLE associate_study (
	study1_id int(10) unsigned not null,
	study2_id int(10) unsigned not null,
	
	primary key( study1_id,study2_id )
);

/**
@table submitter

@colour #7CFC00
@desc This table contains descriptions of group submitting data to public repositories such as ClinVar

@column submitter_id            Primary key
@column description             Description of data submitter

@see phenotype_feature_attrib
*/
CREATE TABLE submitter(
  submitter_id int(10) unsigned not null auto_increment,
  description          varchar(255),
  PRIMARY KEY ( submitter_id )
);


/**
@table publication
@colour #7CFC00
@desc This table contains details of publications citing variations.
      This information comes from dbSNP, UCSC and Europe PMC.

@column publication_id  Primary key, internal identifier.
@column title           Title of the publication
@column authors         Authors of the publication
@column pmid            The PubMed id for the publication if available
@column pmcid           The PubMed Central id for the publication if available
@column year            The year the publication was published
@column doi             The DOI (Digital Object Identifier) for the publication
@column ucsc_id         The external id used in the UCSC database and URL

@see variation_citation
*/

CREATE TABLE publication(
  publication_id int(10) unsigned not null auto_increment, 
  title          varchar(255),
  authors        varchar(255) CHARACTER SET utf8mb4,
  pmid           int(10),
  pmcid          varchar(255),
  year           int(10) unsigned,
  doi            varchar(50),
  ucsc_id        varchar(50),
  PRIMARY KEY ( publication_id ),
  KEY pmid_idx (pmid),
  KEY doi_idx (doi)
);

/**
@table variation_citation
@colour #7CFC00
@desc This table links a variation to a publication

@column publication_id       Primary key, internal identifier.
@column variation_id         Primary key, foreign key references @link variation

@see publication
@see variation
*/

CREATE TABLE variation_citation (
   variation_id int(10) unsigned not null,
   publication_id int(10) unsigned not null,
   PRIMARY KEY variation_citation_idx (variation_id, publication_id)
);



/**
@header  Metadata tables
@desc    These tables define some metadata information.
@colour  #BC5CEC
*/


/**
@table meta_coord

@colour #BC5CEC
@desc This table gives the coordinate system used by various tables in the database.

@column table_name			Name of the feature table, e.g. "variation_feature".
@column coord_system_id		Foreign key to core database coord_system table refers to coordinate system that features from this table can be found in.
@column max_length			Maximum length of the feature. 
*/

CREATE TABLE meta_coord (

  table_name      VARCHAR(40) NOT NULL,
  coord_system_id INT(10) UNSIGNED NOT NULL,
  max_length		  INT,

  UNIQUE(table_name, coord_system_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table meta

@colour #BC5CEC
@desc This table stores various metadata relating to the database, generally used by the Ensembl web code.

@column meta_id		Primary key, internal identifier.
@column species_id	...
@column meta_key		Name of the meta entry, e.g. "schema_version".
@column meta_value	Corresponding value of the key, e.g. "61".
*/

CREATE TABLE meta (

  meta_id 		INT(10) UNSIGNED not null auto_increment,
  species_id  INT UNSIGNED DEFAULT 1,
  meta_key    varchar( 40 ) not null,
  meta_value  varchar( 255 ) not null,

  PRIMARY KEY( meta_id ),
  UNIQUE KEY species_key_value_idx (species_id, meta_key, meta_value ),
  KEY species_value_idx (species_id, meta_value )

);


# Add schema type and schema version to the meta table.
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'schema_type', 'variation'), (NULL, 'schema_version', '92');


# Patch IDs for new release
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_91_92_a.sql|schema version');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_91_92_b.sql|Add index on feature_stable_id in regulatory_feature_variation and motif_feature_variation');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_91_92_c.sql|Change the column zygosity in the table structural_variation_sample');
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_91_92_d.sql|create table to store names of submitters to ClinVar etc');

/**
@header  Failed tables
@desc    These tables define the list of variants/alleles flagged as "failed" in the Variation pipeline.<br />
         The list of reasons for a variation being flagged as failed is available in the <a href="/info/genome/variation/data_description.html#quality_control">Quality Control documentation</a>.
@colour  #3CB371
*/


/**
@table failed_description

@colour #3CB371
@desc This table contains descriptions of reasons for a variation being flagged as failed.

@column failed_description_id	Primary key, internal identifier.
@column description				Text containing the reason why the Variation has been flagged as failed. e.g. "Variation does not map to the genome".

@example See below the list of the descriptions available in the Ensembl variation databases:
         @sql SELECT * FROM failed_description;

@see failed_variation
@see failed_allele
@see failed_structural_variation
*/

CREATE TABLE failed_description(

 failed_description_id int(10) unsigned not null AUTO_INCREMENT,
 description  text not null,

 PRIMARY KEY (failed_description_id)
);


/**
@table failed_variation

@colour #3CB371
@desc For various reasons it may be necessary to store information about a variation that has failed quality checks in the Variation pipeline. This table acts as a flag for such failures.

@column failed_variation_id		Primary key, internal identifier.
@column variation_id					Foreign key references to the @link variation table.
@column failed_description_id		Foreign key references to the @link failed_description table.

@see failed_description
@see variation
*/

CREATE TABLE failed_variation (
  failed_variation_id int(11) NOT NULL AUTO_INCREMENT,
  variation_id int(10) unsigned NOT NULL,
  failed_description_id int(10) unsigned NOT NULL,
  PRIMARY KEY (failed_variation_id),
  UNIQUE KEY variation_idx (variation_id,failed_description_id)
);


/**
@table failed_allele

@colour #3CB371
@desc Contains alleles that did not pass the Ensembl filters

@column failed_allele_id			Primary key, internal identifier.
@column allele_id						Foreign key references to the @link allele table.
@column failed_description_id		Foreign key references to the @link failed_description table.

@see failed_description
@see allele
*/

CREATE TABLE failed_allele (
  failed_allele_id int(11) NOT NULL AUTO_INCREMENT,
  allele_id int(10) unsigned NOT NULL,
  failed_description_id int(10) unsigned NOT NULL,
  PRIMARY KEY (failed_allele_id),
  UNIQUE KEY allele_idx (allele_id,failed_description_id)
);


/**
@table failed_structural_variation

@colour #3CB371
@desc For various reasons it may be necessary to store information about a structural variation that has failed quality checks (mappings) in the Structural Variation pipeline. This table acts as a flag for such failures.

@column failed_structural_variation_id		Primary key, internal identifier.
@column structural_variation_id					  Foreign key references to the @link structural_variation table.
@column failed_description_id		          Foreign key references to the @link failed_description table.

@see failed_description
@see structural_variation
*/

CREATE TABLE failed_structural_variation (
  failed_structural_variation_id int(11) NOT NULL AUTO_INCREMENT,
  structural_variation_id int(10) unsigned NOT NULL,
  failed_description_id int(10) unsigned NOT NULL,
	
  PRIMARY KEY (failed_structural_variation_id),
  UNIQUE KEY structural_variation_idx (structural_variation_id,failed_description_id)
);

/**
@header  Attributes tables
@desc    These tables define the variation attributes data.
@colour  #FF0000
*/


/**
@table  attrib_type

@colour #FF0000
@desc   Defines the set of possible attribute types used in the attrib table

@column attrib_type_id  Primary key
@column code            A short codename for this type (indexed, so should be used for lookups)
@column name            The name of this type
@column description     Longer description of this type

@example See below the command to display a subset of the attrib_type entries:
         @sql SELECT * FROM attrib_type WHERE attrib_type_id > 468 LIMIT 10;

@see    attrib
@see    attrib_set
*/
CREATE TABLE attrib_type (

    attrib_type_id    SMALLINT(5) UNSIGNED NOT NULL DEFAULT 0,
    code              VARCHAR(20) NOT NULL DEFAULT '',
    name              VARCHAR(255) NOT NULL DEFAULT '',
    description       TEXT,

    PRIMARY KEY (attrib_type_id),
    UNIQUE KEY code_idx (code)
);

/**
@table  attrib

@colour #FF0000
@desc   Defines various attributes used elsewhere in the database

@column attrib_id       Primary key
@column attrib_type_id  Key into the @link attrib_type table, identifies the type of this attribute
@column value           The value of this attribute

@example See below the query to display a subset of the attrib entries:
         @sql SELECT * FROM attrib WHERE attrib_type_id IN (469,470,471) ORDER BY attrib_id LIMIT 21;

@see    attrib_type
@see    attrib_set
*/
CREATE TABLE attrib (

    attrib_id           INT(11) UNSIGNED NOT NULL auto_increment,
    attrib_type_id      SMALLINT(5) UNSIGNED NOT NULL DEFAULT 0,
    value               TEXT NOT NULL,

    PRIMARY KEY (attrib_id),
    UNIQUE KEY type_val_idx (attrib_type_id, value(80))
);

/**
@table  attrib_set

@colour #FF0000
@desc   Groups related attributes together

@column attrib_set_id   Primary key
@column attrib_id       Key of an attribute in this set

@see    attrib
@see    attrib_type
*/
CREATE TABLE attrib_set (

    attrib_set_id       INT(11) UNSIGNED NOT NULL DEFAULT 0,
    attrib_id           INT(11) UNSIGNED NOT NULL DEFAULT 0,

    UNIQUE KEY set_idx (attrib_set_id, attrib_id),
    KEY attrib_idx (attrib_id)
);




/**
@header  Protein tables
@desc    These tables define the protein prediction data.
@colour  #1E90FF
*/


/**
@table  protein_function_predictions

@colour #1E90FF
@desc   Contains encoded protein function predictions for every protein-coding transcript in this species.

@column translation_md5_id  Identifies the MD5 hash corresponding to the protein sequence to which 
                            these predictions apply
@column analysis_attrib_id  Identifies the analysis (sift, polyphen etc.) that produced these predictions 
@column prediction_matrix   A compressed binary string containing the predictions for all possible 
                            amino acid substitutions in this protein. See the explanation <a href="/info/genome/variation/predicted_data.html#nsSNP">here</a>

@see    translation_md5
@see    attrib
*/

CREATE TABLE protein_function_predictions (
    translation_md5_id int(11) unsigned NOT NULL,
    analysis_attrib_id int(11) unsigned NOT NULL,
    prediction_matrix mediumblob,
    
    PRIMARY KEY (translation_md5_id, analysis_attrib_id)
);

/**
@table  translation_md5

@colour #1E90FF
@desc   Maps a hex MD5 hash of a translation sequence to an ID used for the protein function predictions

@column translation_md5_id  Primary key
@column translation_md5     Hex MD5 hash of a translation sequence

@see    protein_function_predictions
*/

CREATE TABLE translation_md5 (
    translation_md5_id int(11) NOT NULL AUTO_INCREMENT,
    translation_md5 char(32) NOT NULL,

    PRIMARY KEY (translation_md5_id),
    UNIQUE KEY md5_idx (translation_md5)
);

/**
@table  protein_function_predictions_attrib

@colour #1E90FF
@desc   Contains information on the data use in protein function predictions

@column translation_md5_id  Identifies the MD5 hash corresponding to the protein sequence to which 
                            these data use in prediction apply
@column analysis_attrib_id  Identifies the analysis (sift, polyphen etc.) that produced these values 
@column attrib_type_id      Key into the @link attrib_type table, identifies the type of this attribute
@column position_values     A compressed binary string containing data relevant to the quality of the predictions

@see    protein_function_predictions
*/
CREATE TABLE protein_function_predictions_attrib (
    translation_md5_id int(11) unsigned NOT NULL,
    analysis_attrib_id int(11) unsigned NOT NULL,
    attrib_type_id     int(11) unsigned NOT NULL,
    position_values    blob,
    
    PRIMARY KEY (translation_md5_id, analysis_attrib_id,attrib_type_id )
);

/**
@legend #FF8500 Tables containing sample, individual, population and genotype data
@legend #01D4F7	Tables containing structural variation data
@legend #FFD700	Tables containing sets of variations
@legend #7CFC00	Tables containing source and study data
@legend #BC5CEC	Tables containing metadata
@legend #3CB371	Tables containing "failed" data
@legend #FF0000	Tables containing attribute data
@legend #1E90FF	Tables concerning protein data
@legend #FF4DC8 Tables concerning the prediction of variation effect(s) in different Ensembl features
@legend #22949b Tables concerning data linked to phenotype
*/

#possible values in the failed_description table

INSERT INTO failed_description (failed_description_id,description) VALUES (1,'Variant maps to more than 3 different locations');
INSERT INTO failed_description (failed_description_id,description) VALUES (2,'None of the variant alleles match the reference allele');
INSERT INTO failed_description (failed_description_id,description) VALUES (3,'Variant has more than 3 different alleles');
INSERT INTO failed_description (failed_description_id,description) VALUES (4,'Loci with no observed variant alleles in dbSNP');
INSERT INTO failed_description (failed_description_id,description) VALUES (5,'Variant does not map to the genome');
INSERT INTO failed_description (failed_description_id,description) VALUES (6,'Variant has no genotypes');
INSERT INTO failed_description (failed_description_id,description) VALUES (7,'Genotype frequencies do not add up to 1');
INSERT INTO failed_description (failed_description_id,description) VALUES (8,'Variant has no associated sequence');
INSERT INTO failed_description (failed_description_id,description) VALUES (9,'Variant submission has been withdrawn by the 1000 genomes project due to high false positive rate');
INSERT INTO failed_description (failed_description_id,description) VALUES (11,'Additional submitted allele data from dbSNP does not agree with the dbSNP refSNP alleles'); 
INSERT INTO failed_description (failed_description_id,description) VALUES (12,'Variant has more than 3 different submitted alleles');         
INSERT INTO failed_description (failed_description_id,description) VALUES (13,'Alleles contain non-nucleotide characters');  
INSERT INTO failed_description (failed_description_id,description) VALUES (14,'Alleles contain ambiguity codes');  
INSERT INTO failed_description (failed_description_id,description) VALUES (15,'Mapped position is not compatible with reported alleles');
INSERT INTO failed_description (failed_description_id,description) VALUES (16,'Flagged as suspect by dbSNP');
INSERT INTO failed_description (failed_description_id,description) VALUES (17,'Variant can not be re-mapped to the current assembly');
INSERT INTO failed_description (failed_description_id,description) VALUES (18,'Supporting evidence can not be re-mapped to the current assembly');
INSERT INTO failed_description (failed_description_id,description) VALUES (19,'Variant maps to more than one genomic location');
INSERT INTO failed_description (failed_description_id,description) VALUES (20,'Variant at first base in sequence');

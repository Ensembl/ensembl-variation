#
# variation
#
# Central table containing actual variations (indels, SNPs etc.)  
#

# variation_id        - primary key, internal identifier
# source_id           - foreign key ref source
# name                - identifier for the variation such as the dbSNP
#                       refSNP id (rs#) or SubSNP id (ss#)
# SNPAncestralAllele  - taken from dbSNP to show ancestral allele for the variation
# failed_description_id - foreign key to failed_description table

create table variation (
	variation_id int(10) unsigned not null auto_increment, # PK
	source_id int(10) unsigned not null, 
	name varchar(255),
	validation_status SET('cluster','freq','submitter','doublehit','hapmap','1000Genome','failed','precious'),
	ancestral_allele text,

	primary key( variation_id ),
	unique ( name ),
	key source_idx (source_id)
);

#
# variation_annotation
#
# Table containing annotation associated with the variation
# such as GWAS
#

create table variation_annotation (
        variation_annotation_id int(10) unsigned not null auto_increment,
        variation_id int(10) unsigned not null,
        phenotype_id int(10) unsigned not null,
        source_id int(10) unsigned not null,
	study varchar(30) default NULL,
        study_type set('GWAS'),
        local_stable_id varchar(255),
        associated_gene varchar(255) default NULL,
        associated_variant_risk_allele varchar(255) default NULL,
        variation_names varchar(255) default NULL,
  	risk_allele_freq_in_controls varchar(30) default NULL,
  	p_value varchar(20) default NULL,
        
	primary key (variation_annotation_id),
        key variation_idx(variation_id),
        key phenotype_idx(phenotype_id),
        key source_idx(source_id)
);

#
# phenotype
# containing phenotype name and description, from GWAS so far
#
#

create table phenotype (
        phenotype_id int(10) unsigned not null auto_increment,
        name varchar(50),
        description varchar(255),

        primary key (phenotype_id),
        unique key name_idx(name)
);

#
# variation_synonym
#
# Table containing alternate identifiers for the same variation.
# For example this might be subsnp identifiers for the refsnp.
#
#

create table variation_synonym (
  variation_synonym_id int(10) unsigned not null auto_increment,
  variation_id int(10) unsigned not null,
  subsnp_id int(15) unsigned ,
  source_id int(10) unsigned not null,
  name varchar(255),
  moltype varchar(50),

  primary key(variation_synonym_id),
  key variation_idx (variation_id),
  key subsnp_idx(subsnp_id),
  unique (name, source_id),
  key source_idx (source_id)
);


#
# sample_synonym
#
# Table containing alternate identifiers for the same sample.
# For example this might be pop_id identifiers for the population in dbSNP
# or individual id identifiers for the individual in dbSNP.
#
#

create table sample_synonym (
  sample_synonym_id int(10) unsigned not null auto_increment,
  sample_id int(10) unsigned not null,
  source_id int(10) unsigned not null,
  name varchar(255),

  primary key(sample_synonym_id),
  key sample_idx (sample_id),
  key (name, source_id)
);

#
# subsnp_handle
#
# Table containing information about subsnp_id and submitter handle
#

create table subsnp_handle (
  subsnp_id int(11) unsigned not null,
  handle varchar(20),

  primary key(subsnp_id)
);

#
# allele
#
# Every allele for every variation in the database has a row in this table.
# Alleles are repeated in this table as often as necessary.  For example
# a variation may have alleles 'A' and 'T'. This would be represented by
# two rows in this table.  A different variation which also had an 'A'
# allele would require another row in this table. This way it is 
# simple to track frequency and population for each allele and hopefully not 
# too much space is wasted on the actual allele strings.
#

# allele_id     - primary key, internal identifier
# variation_id  - foreign key ref variation
# allele        - string representing an allele.  E.g. 'A', 'T'
# frequency     - the frequency of this allele in population 
# sample_id     - foreign key ref population

create table allele(
	allele_id int(10) unsigned not null auto_increment,
	variation_id int(10) unsigned not null,
        subsnp_id int(15) unsigned,
	allele varchar(255),
	frequency float,
	sample_id int(10) unsigned,

	primary key( allele_id ),
        key subsnp_idx(subsnp_id),
	key variation_idx( variation_id,allele(10) )
);

#
# sample
#
# A base class to merge the individual and population or assay in a more general]
# concept, basically to have a unique sample_id

# sample_id      - primary key, internal identifier
# name           - name or identifier of the sample
# size           - if the size is NULL its not known or not relevant for this sample
#                  eg. "european" would not have a size 
# description    - free text that describes the sample
# display		 - determines whether the sample appears by default, optionally or not at all

create table sample(
	sample_id int(10) unsigned not null auto_increment,
	name varchar(255) not null,
	size int,
	description text,
	display enum('REFERENCE','DEFAULT','DISPLAYABLE','UNDISPLAYABLE','LD') default 'UNDISPLAYABLE',

	primary key( sample_id ),
	key name_idx( name )
);

#
# population
#
# A population may be an ethnic group (e.g. caucasian, hispanic), assay group (e.g. 24 europeans),
# strain, phenotypic group (e.g. blue eyed, diabetes) etc. 
# Populations may be composed of other populations by defining relationships in the 
# population_structure table.
#

# sample_id            - primary key, internal identifier

create table population(
	sample_id int(10) unsigned not null,

	primary key( sample_id )
);


#
# population_structure
#
# Defines sub/super population relationships.  For example an assay used to determine
# allele frequency may be represented by a superpopulation of caucasions and a sub population 
# of the group of people used in the assay.
#
create table population_structure (
  super_population_sample_id int(10) unsigned not null,
  sub_population_sample_id int(10) unsigned not null,

  unique(super_population_sample_id, sub_population_sample_id),
  key sub_pop_sample_idx (sub_population_sample_id, super_population_sample_id)
);


#
# individual
#
# Table containing individuals.  An individual is a single member of a population.
#
#  sample_id             - PK, unique internal identifier
#  gender                - the sex of this individual
#  father_individual_id  - self referential id, the father of this individual if known
#  mother_individual_id  - self referential id, the mother of this individual if known
#  
#

create table individual(
  sample_id int(10) unsigned not null,
  gender enum('Male', 'Female', 'Unknown') default 'Unknown' NOT NULL,
  father_individual_sample_id int(10) unsigned,
  mother_individual_sample_id int(10) unsigned,
  individual_type_id int(10) unsigned not null,

  primary key(sample_id)
);

#
# individual_type
#
# Table containing the different types of individuals depending on the specie
#

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

#
# variation_feature
#
# This is a feature table similar to the feature tables in the core database.
# The seq_region_id references a seq_region in the core database and the
# seq_region_start, seq_region_end and seq_region_strand represent a 
# variation position on that seq_region.  This table incorporates some 
# denormalisation, taking fields from other tables so that information
# needed for feature creation can be quickly retrieved.
#
# variation_feature_id  - primary key, internal identifier
# seq_region_id         - foreign key references seq_region in core db
#                         This refers to the seq_region which this snp is
#                         on, which may be a chromosome or clone etc.
# seq_region_start      - the start position of the variation on the seq_region
# seq_region_end        - the end position of the variation on the seq_region
# seq_region_strand     - the orientation of the variation on the seq_region
# variation_id          - foreign key refs variation, the variation associated
#                         with this position
# allele_string         - this is a denormalised string taken from the 
#                         alleles in the allele table associated with this
#                         variation.  The reference allele (i.e. one on the
#                         reference genome comes first).
# variation_name        - a denormalisation taken from the variation table
#                         this is the name or identifier that is used for
#                         displaying the feature.
# map_weight            - the number of times that this variation has mapped 
#                         to the genome.  This is a denormalisation as this
#                         particular feature is one example of a mapped 
#                         location.  This can be used to limit the 
#                         the features that come back from a query.
# flags                 - possible values genotyped, to filter the selection of
#			  variations


create table variation_feature(
	variation_feature_id int(10) unsigned not null auto_increment,
	seq_region_id int(10) unsigned not null,
	seq_region_start int not null,
	seq_region_end int not null,
	seq_region_strand tinyint not null,
	variation_id int(10) unsigned not null,
	allele_string text,
        variation_name varchar(255),
	map_weight int not null,
	flags SET('genotyped'),
	source_id int(10) unsigned not null, 
	validation_status SET('cluster','freq','submitter','doublehit','hapmap','1000Genome','precious'),
	consequence_type SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','COMPLEX_INDEL',
	                      'FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','PARTIAL_CODON','SYNONYMOUS_CODING',
				    'REGULATORY_REGION','WITHIN_MATURE_miRNA','5PRIME_UTR','3PRIME_UTR','INTRONIC','NMD_TRANSCRIPT','UPSTREAM','DOWNSTREAM',
				    'WITHIN_NON_CODING_GENE','NO_CONSEQUENCE','INTERGENIC','HGMD_MUTATION')
	default "INTERGENIC" not null ,	
	primary key( variation_feature_id ),
	key pos_idx( seq_region_id, seq_region_start ),
	key variation_idx( variation_id )
);


#
# structural_variation
#
# This is a feature table similar to the variation_feature table above.
# It specifically stores structural variations, and has several
# differences with the columns seen in the variation_feature table.
# The seq_region_id references a seq_region in the core database and the
# seq_region_start, seq_region_end and seq_region_strand represent a 
# variation position on that seq_region. This table also has columns for
# bound_start and bound_end, which represent the outermost boundaries of the
# feature submitted. This table is not linked to variation as above.
#
# structural_variation_id  - primary key, internal identifier
# seq_region_id         - foreign key references seq_region in core db
#                         This refers to the seq_region which this snp is
#                         on, which may be a chromosome or clone etc.
# seq_region_start      - the start position of the variation on the seq_region
# seq_region_end        - the end position of the variation on the seq_region
# seq_region_strand     - the orientation of the variation on the seq_region
# bound_start           - the 5'-most bound defined for the feature
# bound_end             - the 3'-most bound defined for the feature
# variation_name        - the external identifier or name of the variation
# class                 - the type of structural variation feature e.g. 'CNV'

CREATE TABLE structural_variation (
  structural_variation_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  seq_region_id int(10) unsigned NOT NULL,
  seq_region_start int(11) NOT NULL,
  seq_region_end int(11) NOT NULL,
  seq_region_strand tinyint(4) NOT NULL,
  variation_name varchar(255) DEFAULT NULL,
  source_id int(10) unsigned NOT NULL,
  class varchar(255) DEFAULT NULL,
  bound_start int(11) DEFAULT NULL,
  bound_end int(11) DEFAULT NULL,
  PRIMARY KEY (structural_variation_id),
  KEY pos_idx (seq_region_id,seq_region_start),
  KEY name_idx (variation_name)
);


# A table for mapping variations to variation_sets
CREATE TABLE IF NOT EXISTS variation_set_variation (
	variation_id int(10) unsigned NOT NULL,
	variation_set_id int(10) unsigned NOT NULL,
	PRIMARY KEY (variation_id,variation_set_id),
	KEY variation_set_idx (variation_set_id,variation_id)
);

# A table containing variation_set information  
CREATE TABLE IF NOT EXISTS variation_set (
	variation_set_id int(10) unsigned NOT NULL AUTO_INCREMENT,
	name VARCHAR(255),
	description TEXT,
	PRIMARY KEY (variation_set_id),
	KEY name_idx (name)
);

 
# A table containing relashionship between variation sets
CREATE TABLE IF NOT EXISTS variation_set_structure (
	variation_set_super int(10) unsigned NOT NULL,
	variation_set_sub int(10) unsigned NOT NULL,
	PRIMARY KEY (variation_set_super,variation_set_sub),
	KEY sub_idx (variation_set_sub,variation_set_super)
);



#
# variation_group
# 
# This table defines a variation group.  Allele groups can 
# be classed into a variation_group when they are comprised
# of the same set of variations.  This is equivalent to 
# HapSets in dbSNP.
#
# variation_group_id   - primary_key, internal identifier
# name                 - the code or name of this variation_group
# source_id            - foreign key ref source
#

create table variation_group (
	variation_group_id int(10) unsigned not null auto_increment,
	name varchar(255),
	source_id int(10) unsigned not null,
  type enum('haplotype', 'tag'),

	primary key (variation_group_id),
  unique(name)
);

#
# variation_group_variation
# Keeps all the variations stored in a group (normalisation of the n..n relationship between variation and variation_group)
# variation_id - foreign key references variation
# variation_group_id - foreign key references variation_group
#

create table variation_group_variation (
	variation_id int(10) unsigned not null,
	variation_group_id int(10) unsigned not null,

	unique( variation_group_id, variation_id ),
	key variation_idx( variation_id, variation_group_id )
);

#
# variation_group_feature
#
# Places a variation_group (i.e. group of associated haplotypes) on the genome
# as a feature.
#
# variation_group_feature_id - primary key, internal identifier
# seq_region_id              - foreign key references seq_region in core db
# seq_region_start           - start position of the variation_group_feature
#                              on the referenced seq_region
# seq_region_end             - end position of the variation_group_feature
# seq_region_strand          - orientation of feature on seq_region
# variation_group_id         - foreign key references variation_group
# variation_group_name       - name of the variation_group
#

create table variation_group_feature(
  variation_group_feature_id int(10) unsigned not null auto_increment,
  seq_region_id int(10) unsigned not null,
  seq_region_start int not null,
  seq_region_end int not null,
  seq_region_strand tinyint not null,
  variation_group_id int(10) unsigned not null,
  variation_group_name varchar(255),

  primary key (variation_group_feature_id),
  key pos_idx(seq_region_id, seq_region_start),
  key variation_group_idx(variation_group_id)
);

#
# transcript_variation
# 
# This table contains a classification of variation features based on Ensembl
# predicted transcripts.  Variation features which fall into Ensembl 
# transcript regions are classified as 'ESSENTIAL_SPLICE_SITE','SPLICE_SITE',
# 'FRAMESHIFT_CODING','STOP_GAINED','STOP_LOST','NON_SYNONYMOUS_CODING','PARTIAL_CODON',
# 'SYNONYMOUS_CODING','REGULATORY_REGION','WITHIN_MATURE_miRNA','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM'
# 'WITHIN_NON_CODING_GENE', 'COMPLEX_INDEL'

#
# transcript_variation_id - primary key, internal identifier
# transcript_stable_id    - foreign key to core databases
#                           unique stable id of related transcript
# variation_feature_id    - foreign key ref variation_feature
# cdna_start              - start position of variation in cdna coordinates
# cdna_end                - end position of variation in cdna coordinates
# cds_start               - start position of variation in cds coordinates
# cds_end               - end position of variation in cds coordinates
# translation_start       - start position of variation on peptide
# translation_end         - end position of variation on peptide
# peptide_allele_string   - allele string of '/' separated amino acids
# consequence_type        - reference allele is first
# 

create table transcript_variation(
  transcript_variation_id int(10) unsigned not null auto_increment,
  transcript_stable_id varchar(128) not null,
  variation_feature_id int(10) unsigned not null,
  cdna_start int,
  cdna_end   int,
  cds_start  int,
  cds_end    int,
  translation_start int,
  translation_end int,  
  peptide_allele_string varchar(255),
  consequence_type SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','COMPLEX_INDEL',
			     'FRAMESHIFT_CODING', 'NON_SYNONYMOUS_CODING','SPLICE_SITE','PARTIAL_CODON','SYNONYMOUS_CODING',
			     'REGULATORY_REGION','WITHIN_MATURE_miRNA','5PRIME_UTR','3PRIME_UTR','INTRONIC','NMD_TRANSCRIPT','UPSTREAM',
			     'DOWNSTREAM','WITHIN_NON_CODING_GENE','HGMD_MUTATION') not null,	 
  primary key( transcript_variation_id ),
  key variation_idx( variation_feature_id ),
  key transcript_stable_idx( transcript_stable_id ),
  key consequence_type_idx(consequence_type)
	);
	

CREATE TABLE seq_region (

  seq_region_id               INT(10) UNSIGNED NOT NULL,
  name                        VARCHAR(40) NOT NULL,

  PRIMARY KEY (seq_region_id),
  UNIQUE KEY name_idx (name)

) ;


#
# allele_group
#
# This table defines haplotypes - groups of
# polymorphisms which are found together in a block.
# This is equivalent to Haps in dbSNP
#
# allele_group_id    - primary_key, internal identifier
# variation_group_id - foreign key, ref variation_group
# sample_id          - foreign key, ref population
# name               - the name of this allele group
# frequency          - the frequency of this allele_group
#                      within the referenced population
#
#

create table allele_group(
	allele_group_id int(10) unsigned not null auto_increment,
	variation_group_id int(10) unsigned not null,
	sample_id int(10) unsigned,
	name varchar(255),
	source_id int(10) unsigned,
	frequency float,

	primary key( allele_group_id ),
  unique(name)
);


#
# allele_group_allele
#
# Defines which alleles make up an allele group.  
# There is no direct link to the allele table because the allele table has 
# population and frequency data which may not correspond to this allele group
#
# allele_group_id - primary key, internal identifier
# allele - base present in the group
# variation_id - foreign key, references variation

create table allele_group_allele (
	allele_group_id int(10) unsigned not null,
	allele varchar(255) not null,
        variation_id int(10) unsigned not null,

	unique( allele_group_id, variation_id ),
	key allele_idx( variation_id, allele_group_id )
);


#
# flanking_sequence
#
# table that stores the flanking sequences from th core database. To reduce space used, takes coordinates from the sequences in the core database
# variation_id - primary key, internal identifier
# up_seq - upstream sequence, used to initially store the sequence from the core database, and in a later process get from here the position
# down_seq - similiar the one before, but for the downstream
# up_seq_region_start, down_seq_region_start - position of the starting of the sequence in the region
# up_seq_region_end, down_seq_region_end - position of the end of the sequence in the region
# seq_region_id - foreign key, references the sequence table in the core database
# seq_region_stran - strand of the seq_region in the core database
#

create table flanking_sequence (
	variation_id int(10) unsigned not null,
	up_seq text,
	down_seq text,
  	up_seq_region_start int,
  	up_seq_region_end   int,
  	down_seq_region_start int,
  	down_seq_region_end int,
  	seq_region_id int(10) unsigned,
  	seq_region_strand tinyint,

  primary key( variation_id )

) MAX_ROWS = 100000000;


#
# httag
#
# this table contains the tags of a haplotype: bases of the haplotypes that uniquely identify it
#
# httag_id - primary key, internal identifier
# variation_group_id - foreign key, references variation_group
# name - name of the tag, for web purposes
# source_id - foreign key, references source
#

create table httag(
	httag_id int(10) unsigned not null auto_increment,
	variation_group_id int(10) unsigned not null,
	name varchar(255),
	source_id int(10) unsigned not null,

	primary key( httag_id ),
	key variation_group_idx( variation_group_id )
);

#
# source
#
# this table contains sources of snps. this might be dbSNP, TSC, HGBase, etc. 
#
# source_id 	- primary key, internal identifier
# name      	- the name of the source.  e.g. 'dbSNP'
# description	- description of the source for ContigView
# url           - URL for the source
# somatic       - flag to indicate somatic mutations

create table source(
	source_id int(10) unsigned not null auto_increment,
	name varchar(255),
	version int,
	description varchar(255),
	url varchar(255),
	somatic tinyint(1) DEFAULT '0',
	
	primary key( source_id )
);



#
# population_genotype
#
# This table contains genotype frequencies estimated for populations or calculated on
# a set of individuals.
#
# population_genotype_id - primary key, internal identifier
# variation_id - foreign key, references variation table
# allele_1 - first allele in the genotype
# allele_2 - second allele in the genotype
# frequency - frequency of the genotype in the population
# sample_id - foreign key, references population table
#

create table population_genotype (
	population_genotype_id int(10) unsigned not null auto_increment,
	variation_id int(10) unsigned not null,
        subsnp_id int(15) unsigned not null,
	allele_1 varchar(255),
	allele_2 varchar(255),
	frequency float,
 	sample_id int(10) unsigned,

	primary key( population_genotype_id ),
 	key variation_idx(variation_id),
        key subsnp_idx(subsnp_id),
	key sample_idx(sample_id)
);

#
# individual_population
#
# This table contains the relations between individuals and populations (n to n relationship)
#
# individual_sample_id - FK to individual table
# population_sample_id - FK to population table

create table individual_population (
  individual_sample_id int(10) unsigned not null,
  population_sample_id int(10) unsigned not null,

  key individual_sample_idx(individual_sample_id),
  key population_sample_idx(population_sample_id)

);

#This table is only needed for create master schema when run healthcheck system
#needed for other species, but human, so keep it

CREATE TABLE tmp_individual_genotype_single_bp (
                            variation_id int(10) not null,
			    subsnp_id int(15) unsigned,   
	                    allele_1 varchar(255),allele_2 varchar(255),sample_id int,
                            key variation_idx(variation_id),
                            key subsnp_idx(subsnp_id),
                            key sample_idx(sample_id)
                            ) MAX_ROWS = 100000000
;

#
# individual_genotype_multiple_bp
#
# This table contains genotypes of individuals with more than 1 bp in the alleles.
#
# variation_id	- FK to variation table
# allele_1	- One of the alleles of the genotype
# allele_2	- The other allele of the genotype
# sample_id     - foreign key, references individual table

create table individual_genotype_multiple_bp (
  variation_id int(10) unsigned not null,
  subsnp_id int(15) unsigned,	
  allele_1 varchar(255),
  allele_2 varchar(255),
  sample_id int(10) unsigned,

  key variation_idx(variation_id),
  key subsnp_idx(subsnp_id),
  key sample_idx(sample_id)
);


#
# meta_coord
#
# Same table structure as in core database. Contains info about what coord
# systems features can be found in.
#
# table_name - name of the feature table
# coord_system_id - foreign key to core database coord_system table
#                   refers to coord system that features from this table can
#                   be found in
#

CREATE TABLE meta_coord (

  table_name                  VARCHAR(40) NOT NULL,
  coord_system_id             INT(10) UNSIGNED NOT NULL,
  max_length		      INT,

  UNIQUE(table_name, coord_system_id)

) TYPE=MyISAM;


################################################################################
#
# Table structure for table 'meta' 
#

CREATE TABLE meta (

  meta_id 		      INT(10) UNSIGNED not null auto_increment,
  species_id                  INT UNSIGNED DEFAULT 1,
  meta_key                    varchar( 40 ) not null,
  meta_value                  varchar( 255 ) not null,

  PRIMARY KEY( meta_id ),
  UNIQUE KEY species_key_value_idx (species_id, meta_key, meta_value ),
  KEY species_value_idx (species_id, meta_value )

);


# insert schema type row
INSERT INTO meta (meta_key, meta_value) VALUES ('schema_type', 'variation');


###############
#
#  Table structure for table tagged_variation_features
#
###############

CREATE TABLE tagged_variation_feature (

  variation_feature_id       INT(10) UNSIGNED not null,
  sample_id              INT(10) UNSIGNED not null,
  
  PRIMARY KEY(variation_feature_id, sample_id)
);

###############
#
# Table structure for table read_coverage
#
###############

CREATE TABLE read_coverage (
   seq_region_id int(10) unsigned not null,
   seq_region_start int not null,
   seq_region_end int not null,
   level tinyint not null,
   sample_id int(10) unsigned not null,
		  
   key seq_region_idx(seq_region_id,seq_region_start)   
);

################
#
# Table structure for table compressed_genotype_single_bp
#
################

CREATE TABLE compressed_genotype_single_bp(
  sample_id int(10) unsigned not null,
  seq_region_id int(10) unsigned not null,
  seq_region_start int not null,
  seq_region_end int not null,
  seq_region_strand tinyint not null,
  genotypes blob,

  key pos_idx(seq_region_id,seq_region_start)
) MAX_ROWS = 100000000;

#
# failed_description
#
# Contains reasons for removing some variations from the Variation database
#
# failed_description_id  - primary key, internal identifier
# description - text containing the reason why the Variation information has been removed from the 
#               Variation databse except in the Variation table
#

CREATE TABLE failed_description(

 failed_description_id int(10) unsigned not null,
 description  text not null,

 PRIMARY KEY (failed_description_id)
);

#
# failed_variation
#
# Contains all variations that did not pass the Ensembl filters
# variation_id - primary key
# failed_descriptin_id - foreign key to failed_description table
#

CREATE TABLE failed_variation(
    variation_id int(10) unsigned not null,
    failed_description_id int(10) unsigned not null,

    PRIMARY KEY(variation_id)
);

#
# strain_gtype_poly
#
# This table is populated for mouse and rat only for mart to use. Mart build need this table in both staging sever in master_schema_variation database
#

CREATE TABLE strain_gtype_poly (
  variation_id int(10) unsigned NOT NULL,
  sample_name varchar(100) DEFAULT NULL,
  
  PRIMARY KEY (variation_id)
);

#possible values in the failed_description table
INSERT INTO failed_description (failed_description_id,description) VALUES (1,'Variation maps to more than 3 different locations');
INSERT INTO failed_description (failed_description_id,description) VALUES (2,'None of the variant alleles match the reference allele');
INSERT INTO failed_description (failed_description_id,description) VALUES (3,'Variation has more than 3 different alleles');
INSERT INTO failed_description (failed_description_id,description) VALUES (4,'Loci with no observed variant alleles in dbSNP');
INSERT INTO failed_description (failed_description_id,description) VALUES (5,'Variation does not map to the genome'); 

/**
@table variation

@desc This is the schema's generic representation of a variation, defined as a genetic feature that varies between individuals of the same species.The most common type is the single nucleotide variation (SNP) though the schema also accommodates copy number variations (CNVs) and structural variations (SVs).A variation is defined by its flanking sequence rather than its mapped location on a chromosome; a variation may in fact have multiple mappings across a genome.This table stores a variation's name (commonly an ID of the form e.g. rs123456, assigned by dbSNP), along with a validation status and ancestral (or reference) allele.

@column variation_id		Primary key, internal identifier.
@column source_id			Foreign key references to the @link source table.
@column name				Name of the variation. e.g. "rs1333049".
@column validation_status	Variant discovery method and validation from dbSNP.
@column ancestral_allele	Taken from dbSNP to show ancestral allele for the variation.
@column flipped				This is set to 1 if the variant is flipped from the negative to the positive strand during import.
@column class_attrib_id		Class of the variation, key into the @link attrib table
@column somatic             flags whether this variation is known to be somatic or not

@see variation_synonym
@see flanking_sequence
@see failed_variation
@see variation_feature
@see variation_group_variation
@see allele
@see allele_group_allele
@see individual_genotype_multiple_bp
@see attrib
*/

create table variation (
	variation_id int(10) unsigned not null auto_increment, # PK
	source_id int(10) unsigned not null, 
	name varchar(255),
	validation_status SET('cluster','freq',
								 'submitter','doublehit',
								 'hapmap','1000Genome',
								 'failed','precious'),
	ancestral_allele varchar(255) DEFAULT NULL,
	flipped tinyint(1) unsigned NULL DEFAULT NULL,
	class_attrib_id int(10) unsigned not null default 0,
	somatic tinyint(1) DEFAULT 0 NOT NULL,

	primary key( variation_id ),
	unique ( name ),
	key source_idx (source_id)
);


/**
@table variation_annotation

@desc This table stores information linking genotypes and phenotypes. It stores various fields pertaining to the study conducted, along with the associated gene, risk allele frequency and a p-value.

@column variation_annotation_id					Primary key, internal identifier.
@column variation_id										Foreign key references to the @link variation table.
@column phenotype_id										Foreign key references to the @link phenotype table.
@column study_id												Foreign key references to the @link study table.
@column associated_gene									Common gene(s) name(s) associated to the variation.
@column associated_variant_risk_allele	Allele associated to the phenotype.
@column variation_names									Name of the variation. e.g. "rs1333049".
@column risk_allele_freq_in_controls		Risk allele frequency.
@column p_value													P value of the association phenotype/variation.

@see variation
@see phenotype
@see source
*/

create table variation_annotation (
	variation_annotation_id int(10) unsigned not null auto_increment,
	variation_id int(10) unsigned not null,
	phenotype_id int(10) unsigned not null,
	study_id int(10) unsigned not null,
	associated_gene varchar(255) default NULL,
	associated_variant_risk_allele varchar(255) default NULL,
	variation_names varchar(255) default NULL,
	risk_allele_freq_in_controls double default NULL,
	p_value double default NULL,
	
	primary key (variation_annotation_id),
	key variation_idx(variation_id),
	key phenotype_idx(phenotype_id),
	key study_idx(study_id)
);


/**
@table phenotype

@desc This table stores details of the phenotypes associated with variation annotations.

@column phenotype_id	Primary key, internal identifier.
@column name					Phenotype short name. e.g. "CAD".
@column description	varchar		Phenotype long name. e.g. "Coronary Artery Disease".

@see variation_annotation
*/

create table phenotype (
	phenotype_id int(10) unsigned not null auto_increment,
	name varchar(50),
	description varchar(255),

	primary key (phenotype_id),
	unique key name_idx(name)
);


/**
@table variation_synonym

@desc This table allows for a variation to have multiple IDs, generally given by multiple sources.

@column variation_synonym_id	Primary key, internal identifier.
@column variation_id					Foreign key references to the variation table.
@column subsnp_id							Foreign key references to the subsnp_handle table.
@column source_id							Foreign key references to the source table.
@column name									Name of the synonym variation. e.g. 'rs1333049'. The corresponding variation ID of this variation is different from the one stored in the column variation_id.
@column moltype								...

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
  moltype varchar(50),

  primary key(variation_synonym_id),
  key variation_idx (variation_id),
  key subsnp_idx(subsnp_id),
  unique (name, source_id),
  key source_idx (source_id)
);


/**
@table sample_synonym

@colour #FF8500
@desc Used to store alternative names for populations when data comes from multiple sources.

@column sample_synonym_id	Primary key, internal identifier.
@column sample_id					Foreign key references to the @link sample table.
@column source_id					Foreign key references to the @link source table.
@column name							Name of the synonym (a different <b>sample_id</b>).

@see sample
@see population
@see source
*/

create table sample_synonym (
  sample_synonym_id int(10) unsigned not null auto_increment,
  sample_id int(10) unsigned not null,
  source_id int(10) unsigned not null,
  name varchar(255),

  primary key(sample_synonym_id),
  key sample_idx (sample_id),
  key (name, source_id)
);


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
@table allele

@desc This table stores information about each of a variation's alleles, along with population frequencies.

@column allele_id		Primary key, internal identifier.
@column variation_id	Foreign key references to the @link variation table.
@column subsnp_id		Foreign key references to the @link subsnp_handle table.
@column allele			Allele found at the variation location, for the sample. e.g. "A".
@column frequency		Frequency of this allele in the sample.
@column sample_id		Foreign key references to the @link sample table.
@column count			Number of individuals in the sample where this allele is found.

@see variation
@see population
@see subsnp_handle
*/

create table allele(
	allele_id int(10) unsigned not null auto_increment,
	variation_id int(10) unsigned not null,
   subsnp_id int(15) unsigned,
	allele varchar(25000),
	frequency float,
	sample_id int(10) unsigned,
	count int(10) unsigned DEFAULT NULL,

	primary key( allele_id ),
   key subsnp_idx(subsnp_id),
	key variation_idx( variation_id,allele(10) )
);


/**
@table sample

@colour #FF8500
@desc Sample is used as a generic catch-all term to cover individuals, populations and strains; it contains a name and description, as well as a size if applicable to the population.

@column sample_id		Primary key, internal identifier.
@column name			Name of the sample (can be an individual or a population name).
@column size			Number of individual in the sample.
@column description	Description of the sample.
@column display		Information used by the website: samples with little information are filtered from some web displays.

@see individual
@see population
*/

create table sample(
	sample_id int(10) unsigned not null auto_increment,
	name varchar(255) not null,
	size int,
	description text,
	display enum('REFERENCE',
					 'DEFAULT',
					 'DISPLAYABLE',
					 'UNDISPLAYABLE',
					 'LD',
					 'MARTDISPLAYABLE') default 'UNDISPLAYABLE',

	primary key( sample_id ),
	key name_idx( name )
);


/**
@table population

@colour #FF8500
@desc A table consisting simply of sample_ids representing populations; all data relating to the populations are stored in separate tables (see below).<br />A population may be an ethnic group (e.g. caucasian, hispanic), assay group (e.g. 24 europeans), strain, phenotypic group (e.g. blue eyed, diabetes) etc. Populations may be composed of other populations by defining relationships in the population_structure table.

@column sample_id	int	Foreign key references to the @link sample table. Corresponds to the population ID.

@see sample
@see sample_synonym
@see individual_population
@see  population_structure
@see population_genotype
@see allele
@see allele_group
@see tagged_variation_feature
*/

create table population(
	sample_id int(10) unsigned not null,

	primary key( sample_id )
);


/**
@table population_structure

@colour #FF8500
@desc This table stores hierarchical relationships between populations by relating them as populations and sub-populations.

@column super_population_sample_id	Foreign key references to the population table.
@column sub_population_sample_id		Foreign key references to the population table.

@see population
*/

create table population_structure (
  super_population_sample_id int(10) unsigned not null,
  sub_population_sample_id int(10) unsigned not null,

  unique(super_population_sample_id, sub_population_sample_id),
  key sub_pop_sample_idx (sub_population_sample_id, super_population_sample_id)
);


/**
@table individual

@colour #FF8500
@desc Stores information about an identifiable individual, including gender and the identifiers of the individual's parents (if known).

@column sample_id							Primary key, internal identifier. See the @link sample table. Corresponds to the individual ID.
@column gender								The sex of this individual.
@column father_individual_sample_id	Self referential ID, the father of this individual if known.
@column mother_individual_sample_id	Self referential ID, the mother of this individual if known.
@column individual_type_id				Foreign key references to the @link individual_type table.


@see sample
@see individual_type
@see individual_population
@see individual_genotype_multiple_bp
@see compressed_genotype_single_bp
*/

create table individual(
  sample_id int(10) unsigned not null,
  gender enum('Male', 'Female', 'Unknown') default 'Unknown' NOT NULL,
  father_individual_sample_id int(10) unsigned,
  mother_individual_sample_id int(10) unsigned,
  individual_type_id int(10) unsigned not null,

  primary key(sample_id)
);


/**
@table individual_type

@colour #FF8500
@desc This table resolves the many-to-many relationship between the individual and population tables; i.e. samples may belong to more than one population. Hence it is composed of rows of individual and population identifiers.

@column individual_type_id	Primary key, internal identifier.
@column name					Short name of the individual type. e.g. "fully_inbred","mutant".
@column description			Long name of the individual type.


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
@table variation_feature

@desc This table represents mappings of variations to genomic locations. It stores an allele string representing the different possible alleles that are found at that locus e.g. "A/T" for a SNP, as well as a "worst case" consequence of the mutation. It also acts as part of the relationship between variations and transcripts.

@column variation_feature_id	Primary key, internal identifier.
@column seq_region_id			Foreign key references @link seq_region in core db. Refers to the seq_region which this variant is on, which may be a chromosome, a clone, etc...
@column seq_region_start		The start position of the variation on the @link seq_region.
@column seq_region_end			The end position of the variation on the @link seq_region.
@column seq_region_strand		The orientation of the variation on the @link seq_region.
@column variation_id				Foreign key references to the @link variation table.
@column allele_string			This is a denormalised string taken from the alleles in the allele table associated with this variation. The reference allele (i.e. one on the reference genome comes first).
@column variation_name			A denormalisation taken from the variation table. This is the name or identifier that is used for displaying the feature.
@column map_weight				The number of times that this variation has mapped to the genome. This is a denormalisation as this particular feature is one example of a mapped location. This can be used to limit the the features that come back from a query.
@column flags						Flag to filter the selection of variations.
@column source_id					Foreign key references to the source table.
@column validation_status		SET('cluster', 'freq', 'submitter', 'doublehit', 'hapmap', '1000Genome', 'precious')	Variant discovery method and validation from dbSNP.
@column consequence_type		The SO accession(s) representing the 'worst' consequence(s) of the variation in a transcript or regulatory region
@column variation_set_id		The variation feature can belong to a @link variation_set.
@column class_attrib_id			Class of the variation, key in the @link attrib table
@column somatic                 Flags whether this variation_feature is somatic or germline

@see variation
@see tagged_variation_feature
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
	validation_status SET(
        'cluster',
        'freq',
		'submitter',
        'doublehit',
		'hapmap',
        '1000Genome',
		'precious'
    ),
    consequence_type SET (
        'intergenic_variant',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'complex_change_in_transcript', 
        'stop_lost',
        'coding_sequence_variant',
        'non_synonymous_codon',
        'stop_gained',
        'synonymous_codon',
        'frameshift_variant',
        'nc_transcript_variant',
        'mature_miRNA_variant',
        'NMD_transcript_variant',
        '5_prime_UTR_variant',
        '3_prime_UTR_variant',
        'incomplete_terminal_codon_variant',
        'intron_variant',
        'splice_region_variant',
        '5KB_downstream_variant',
        '500B_downstream_variant',
        '5KB_upstream_variant',
        '2KB_upstream_variant',
        'initiator_codon_change',
        'stop_retained_variant',
        'inframe_codon_gain',
        'inframe_codon_loss',
        'miRNA_target_site_variant',
        'pre_miRNA_variant',
        'regulatory_region_variant',
        'increased_binding_affinity',
        'decreased_binding_affinity',
        'binding_site_variant'
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
    class_attrib_id int(10) unsigned not null default 0,
    somatic tinyint(1) DEFAULT 0 NOT NULL,

   	primary key( variation_feature_id ),
	  key pos_idx( seq_region_id, seq_region_start, seq_region_end ),
	  key variation_idx( variation_id ),
    key variation_set_idx ( variation_set_id )
);


/**
@table structural_variation

@colour #01D4F7
@desc This table stores information about structural variation.

@column structural_variation_id	Primary key, internal identifier.
@column variation_name					The external identifier or name of the variation. e.g. "esv9549".
@column source_id								Foreign key references to the @link source table.
@column study_id								Foreign key references to the @link study table.	
@column class_attrib_id					Foreign key references to the @link attrib table. Defines the type of structural variant. 
@column validation_status				Validation status of the variant.

@see source
@see study
@see attrib
*/

CREATE TABLE structural_variation (
  structural_variation_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  variation_name varchar(255) DEFAULT NULL,
  source_id int(10) unsigned NOT NULL,
  study_id int(10) unsigned DEFAULT NULL,
	class_attrib_id int(10) unsigned NOT NULL DEFAULT 0,
  validation_status ENUM('validated','not validated','high quality'),
	
  PRIMARY KEY (structural_variation_id),
  KEY name_idx (variation_name),
	KEY study_idx (study_id),
	KEY attrib_idx (class_attrib_id)
);


/**
@table supporting_structural_variation

@colour #01D4F7
@desc This table stores the name of the supporting evidence for the structural variants (e.g. DGVa structural variants).

@column supporting_structural_variation_id	Primary key, internal identifier.
@column name																The identifier or name of the supporting evidence.
@column structural_variation_id							Foreign key references to the @link structural_variation table.
@column class_attrib_id					            Foreign key references to the @link attrib table. Defines the allele type of the supporting evidence. 

@see structural_variation
@see attrib
*/

create table supporting_structural_variation (
	supporting_structural_variation_id int(10) unsigned NOT NULL AUTO_INCREMENT,

	name varchar(255) NOT NULL,
	structural_variation_id int(10) unsigned NOT NULL,
	class_attrib_id int(10) unsigned NOT NULL DEFAULT 0,
	
	PRIMARY KEY (supporting_structural_variation_id),
	KEY structural_variation_idx (structural_variation_id),
	KEY attrib_idx (class_attrib_id)
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
@column class_attrib_id					         Foreign key references to the @link attrib table. Defines the type of structural variant.
@column allele_string						         The variant allele, where known.

@see structural_variation
@see source
@see seq_region
@see attrib
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
  class_attrib_id int(10) unsigned NOT NULL DEFAULT 0,
	allele_string longtext DEFAULT NULL,
	
  PRIMARY KEY (structural_variation_feature_id),
	KEY pos_idx( seq_region_id, seq_region_start, seq_region_end ),
	KEY structural_variation_idx (structural_variation_id),
	KEY attrib_idx (class_attrib_id)
);


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
@desc This table containts the name of sets and subsets of variations stored in the database. It usually represents the name of the project or subproject where a group of variations has been identified.

@column variation_set_id			Primary key, internal identifier.
@column name									Name of the set e.g. "Phenotype-associated variations".
@column description						Description of the set.
@column short_name_attrib_id	Foreign key references to the @link attrib table. Short name used for web purpose.

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
@table variation_group

@desc This table represents the equivalent of a variation for a multi-marker variation.

@column variation_group_id	Primary key, internal identifier.
@column name					The code or name of this variation_group.
@column source_id				Foreign key references to the @link source table.
@column type					Type of the variation group.

@see variation_group_variation
@see allele_group
@see variation
@see source
@see httag
*/

create table variation_group (
	variation_group_id int(10) unsigned not null auto_increment,
	name varchar(255),
	source_id int(10) unsigned not null,
  	type enum('haplotype', 'tag'),

	primary key (variation_group_id),
  unique(name)
);


/**
@table variation_group_variation

@desc This table represents an individual variation that makes up a multi-marker variation, and resolves the many-to-many relationship between variation and variation_group.

@column variation_id				Foreign key references to the @link variation table.
@column variation_group_id		Foreign key references to the @link variation_group table.

@see variation_group
@see variation
*/

create table variation_group_variation (
	variation_id int(10) unsigned not null,
	variation_group_id int(10) unsigned not null,

	unique( variation_group_id, variation_id ),
	key variation_idx( variation_id, variation_group_id )
);


/**
@table variation_group_feature

@desc This table represents the equivalent of a variation_feature for multi-marker variations, mapping a haplotype to a chromosomal coordinate system.

@column variation_group_feature_id	Primary key, internal identifier.
@column seq_region_id					Foreign key references @link seq_region in core db. Refers to the seq_region which this variant is on, which may be a chromosome, a clone, etc...
@column seq_region_start				The start position of the variation on the @link seq_region.
@column seq_region_end					The end position of the variation on the @link seq_region.
@column seq_region_strand				The orientation of the variation on the @link seq_region.
@column variation_group_id				Foreign key references to the @link variation_group table.
@column variation_group_name			Name of the variation group.

@see variation_group
@see seq_region
*/

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


/**
@table transcript_variation

@desc This table relates a single allele of a variation_feature to a transcript (see Core documentation). It contains the consequence of the allele e.g. intron_variant, non_synonymous_codon, stop_lost etc, along with the change in amino acid in the resulting protein if applicable.

@column transcript_variation_id	    Primary key, internal identifier.
@column feature_stable_id		    Foreign key to core databases. Unique stable id of related transcript.
@column variation_feature_id		Foreign key references to the @link variation_feature table.
@column allele_string               Shows the reference sequence and variant sequence of this allele
@column somatic                     Flags if the associated variation is known to be somatic
@column consequence_types			The consequence(s) of the variant allele on this transcript.
@column cds_start					The start position of variation in cds coordinates.
@column cds_end						The end position of variation in cds coordinates.
@column cdna_start					The start position of variation in cdna coordinates.
@column cdna_end					The end position of variation in cdna coordinates.
@column translation_start			The start position of variation on peptide.
@column translation_end				The end position of variation on peptide.
@column codon_allele_string         The reference and variant codons
@column pep_allele_string           The reference and variant peptides
@column hgvs_genomic                HGVS representation of this allele with respect to the genomic sequence
@column hgvs_coding                 HGVS representation of this allele with respect to the CDS
@column hgvs_protein                HGVS representation of this allele with respect to the protein
@column polyphen_prediction         The PolyPhen prediction for the effect of this allele on the protein
@column polyphen_score              The PolyPhen score corresponding to the prediction 
@column sift_prediction             The SIFT prediction for the effect of this allele on the protein 
@column sift_score                  The SIFT score corresponsing to this prediction

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
                                            'complex_change_in_transcript', 
                                            'stop_lost',
                                            'coding_sequence_variant',
                                            'non_synonymous_codon',
                                            'stop_gained',
                                            'synonymous_codon',
                                            'frameshift_variant',
                                            'nc_transcript_variant',
                                            'mature_miRNA_variant',
                                            'NMD_transcript_variant',
                                            '5_prime_UTR_variant',
                                            '3_prime_UTR_variant',
                                            'incomplete_terminal_codon_variant',
                                            'intron_variant',
                                            'splice_region_variant',
                                            '5KB_downstream_variant',
                                            '500B_downstream_variant',
                                            '5KB_upstream_variant',
                                            '2KB_upstream_variant',
                                            'initiator_codon_change',
                                            'stop_retained_variant',
                                            'inframe_codon_gain',
                                            'inframe_codon_loss',
                                            'pre_miRNA_variant'
                                        ),
    cds_start                           int(11) unsigned,
    cds_end                             int(11) unsigned,
    cdna_start                          int(11) unsigned,
    cdna_end                            int(11) unsigned,
    translation_start                   int(11) unsigned,
    translation_end                     int(11) unsigned,
    codon_allele_string                 text,
    pep_allele_string                   text,
    hgvs_genomic                        text,
    hgvs_coding                         text,
    hgvs_protein                        text,
    polyphen_prediction                 enum('unknown', 'benign', 'possibly damaging', 'probably damaging') DEFAULT NULL,
    polyphen_score                      float DEFAULT NULL,
    sift_prediction                     enum('tolerated', 'deleterious') DEFAULT NULL,
    sift_score                          float DEFAULT NULL,

    PRIMARY KEY                         (transcript_variation_id),
    KEY variation_feature_idx           (variation_feature_id),
    KEY feature_idx                     (feature_stable_id),
    KEY consequence_type_idx            (consequence_types),
    KEY somatic_feature_idx             (feature_stable_id, somatic)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;
	
/**
@table seq_region

@desc This table stores the relationship between Ensembl's internal coordinate system identifiers and traditional chromosome names.

@column seq_region_id	Primary key. Foreign key references seq_region in core db. Refers to the seq_region which this variant is on, which may be a chromosome, a clone, etc...
@column name				The name of this sequence region.

@see variation_feature
@see variation_group_feature
@see flanking_sequence
@see compressed_genotype_single_bp
@see read_coverage
*/

CREATE TABLE seq_region (

  seq_region_id               INT(10) UNSIGNED NOT NULL,
  name                        VARCHAR(40) NOT NULL,

  PRIMARY KEY (seq_region_id),
  UNIQUE KEY name_idx (name)

) ;


/**
@table allele_group

@desc This table, along with allele_group_allele, represents a particular multi-marker allele of a given multi-marker variation, or haplotype. It stores an associated population frequency.

@column allele_group_id			Primary key, internal identifier.
@column variation_group_id		Foreign key references to the @link variation_group table.
@column sample_id					Foreign key references to the @link population table.
@column name						The name of this allele group.
@column source_id				Foreign key references to the @link source table.
@column frequency					The frequency of this allele_group within the referenced population.

@see allele_group_allele
@see variation_group
@see population
@see source
*/

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


/**
@table allele_group_allele

@desc This table represents an allele of one variation in a multi-marker variation, or haplotype. It stores a string of the allele.

@column allele_group_id	Primary key, internal identifier.
@column allele				Nucleotid presents in the group.
@column variation_id		Foreign key references to the @link variation table.

@see allele_group
@see variation
*/
create table allele_group_allele (
	allele_group_id int(10) unsigned not null,
	allele varchar(255) not null,
        variation_id int(10) unsigned not null,

	unique( allele_group_id, variation_id ),
	key allele_idx( variation_id, allele_group_id )
);


/**
@table flanking_sequence

@desc This table contains the upstream and downstream sequence surrounding a variation. Since each variation is defined by its flanking sequence, this table has a one-to-one relationship with the variation table.

@column variation_id				Primary key. Foreign key references to the variation table.
@column up_seq						Upstream sequence, used to initially store the sequence from the core database, and in a later process get from here the position.
@column down_seq					Downstream sequence, used to initially store the sequence from the core database, and in a later process get from here the position.
@column up_seq_region_start	Position of the starting of the sequence in the region.
@column up_seq_region_end		Position of the end of the sequence in the region.
@column down_seq_region_start	Position of the starting of the sequence in the region.
@column down_seq_region_end	Position of the end of the sequence in the region.
@column seq_region_id			Foreign key references @link seq_region in core db. Refers to the seq_region which this variant is on, which may be a chromosome or clone etc..
@column seq_region_strand		The orientation of the variation on the seq_region.

@see variation
@see seq_region
*/

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


/**
@table httag

@desc This table represents the equivalent of a tagged_variation_feature for multi-marker variations, representing an instance where a haplotype is tagged by or tags another marker.

@column httag_id				Primary key, internal identifier.
@column variation_group_id	Foreign key references to the @link variation_group table.
@column name					The name of the tag, for web purposes.
@column source_id				Foreign key references to the @link source table.

@see variation_group
@see source
*/

create table httag(
	httag_id int(10) unsigned not null auto_increment,
	variation_group_id int(10) unsigned not null,
	name varchar(255),
	source_id int(10) unsigned not null,

	primary key( httag_id ),
	key variation_group_idx( variation_group_id )
);


/**
@table source

@colour #7CFC00
@desc This table contains details of the source from which a variation is derived. Most commonly this is NCBI's dbSNP; other sources include SNPs called by Ensembl.

@column source_id		Primary key, internal identifier.
@column name			Name of the source. e.g. "dbSNP"
@column version		Version number of the source (if available). e.g. "132"
@column description	Description of the source.
@column url				URL of the source.
@column type			Define the type of the source, e.g. 'chip'
@column somatic_status  Indicates if this source includes somatic or germline mutations, or a mixture

@see variation
@see variation_synonym
@see variation_feature
@see variation_annotation
@see variation_group
@see allele_group
@see sample_synonym
@see httag
@see structural_variation
@see study
*/

create table source(
	source_id int(10) unsigned not null auto_increment,
	name varchar(255),
	version int,
	description varchar(255),
	url varchar(255),
	type ENUM('chip','lsdb') DEFAULT NULL,
    somatic_status ENUM ('germline','somatic','mixed') DEFAULT 'germline',
	
	primary key( source_id )
);


/**
@table study

@colour #7CFC00
@desc This table contains details of the studies.
			The studies information can come from internal studies (DGVa, EGA) or from external studies (Uniprot, NHGRI, ...).

@column study_id						Primary key, internal identifier.
@column source_id						Foreign key references to the @link source table.
@column name								Name of the study. e.g. "EGAS00000000001"
@column description					Description of the study.
@column url									URL to find the study data (http or ftp).
@column external_reference	The pubmed/id or project name associated with this study.
@column study_type					Displays if a study comes from a genome-wide association study or not.

@see source
@see variation_annotation
@see structural_variation
*/

create table study (
	study_id int(10) unsigned not null auto_increment,
	source_id int(10) unsigned not null,
	name varchar(255) DEFAULT null,
	description varchar(255) DEFAULT NULL,
	url varchar(255) DEFAULT NULL,
	external_reference varchar(255) DEFAULT NULL,
	study_type set('GWAS'),
	
	primary key( study_id ),
	key source_idx (source_id)
);


/**
@table associate_study

@colour #7CFC00
@desc This table contains identifiers of associated studies (e.g. NHGRI and EGA studies with the same pubmed identifier).

@column study1_id		Primary key. Foreign key references to the @link study table.
@column study2_id		Primary key. Foreign key references to the @link study table.

@see study
*/
create table associate_study (
	study1_id int(10) unsigned not null,
	study2_id int(10) unsigned not null,
	
	primary key( study1_id,study2_id )
);


/**
@table population_genotype

@colour #FF8500
@desc This table stores alleles and frequencies for variations in given populations.

@column population_genotype_id	Primary key, internal identifier.
@column variation_id					Foreign key references to the @link variation table.
@column subsnp_id						Foreign key references to the subsnp_handle table.
@column allele_1						First allele in the genotype.
@column allele_2						Second allele in the genotype.
@column frequency						Frequency of the genotype in the population.
@column sample_id						Foreign key references to the @link population table.
@column count							Number of individuals who have this genotype, in this population.

@see population
@see variation
@see subsnp_handle
*/

create table population_genotype (
	population_genotype_id int(10) unsigned not null auto_increment,
	variation_id int(10) unsigned not null,
   subsnp_id int(15) unsigned DEFAULT NULL,
	allele_1 varchar(25000),
	allele_2 varchar(25000),
	frequency float,
 	sample_id int(10) unsigned,
	count int(10) unsigned DEFAULT NULL,

	primary key( population_genotype_id ),
 	key variation_idx(variation_id),
    key subsnp_idx(subsnp_id),
	key sample_idx(sample_id)
);


/**
@table individual_population

@colour #FF8500
@desc This table resolves the many-to-many relationship between the individual and population tables; i.e. samples may belong to more than one population. Hence it is composed of rows of individual and population identifiers.

@column individual_sample_id	Foreign key references to the @link individual table.
@column population_sample_id	Foreign key references to the @link population table.

@see individual
@see population
*/

create table individual_population (
  individual_sample_id int(10) unsigned not null,
  population_sample_id int(10) unsigned not null,

  key individual_sample_idx(individual_sample_id),
  key population_sample_idx(population_sample_id)

);


/**
@table tmp_individual_genotype_single_bp

@colour #FF8500
@desc his table is only needed for create master schema when run healthcheck system. Needed for other species, but human, so keep it.

@column variation_id	Primary key. Foreign key references to the @link variation table.
@column subsnp_id		Foreign key references to the @link subsnp_handle table.
@column allele_1		One of the alleles of the genotype, e.g. "TAG".
@column allele_2		The other allele of the genotype.
@column sample_id		Foreign key references to the @link individual table.

@see individual
@see variation
@see subsnp_handle
*/

CREATE TABLE tmp_individual_genotype_single_bp (
	variation_id int(10) not null,
	subsnp_id int(15) unsigned,   
	allele_1 char(1),
	allele_2 char(1),
	sample_id int,

	key variation_idx(variation_id),
   key subsnp_idx(subsnp_id),
   key sample_idx(sample_id)
) MAX_ROWS = 100000000;


/**
@table individual_genotype_multiple_bp

@colour #FF8500
@desc This table holds uncompressed genotypes for given variations.

@column variation_id	Primary key. Foreign key references to the @link variation table.
@column subsnp_id		Foreign key references to the @link subsnp_handle table.
@column allele_1		One of the alleles of the genotype, e.g. "TAG".
@column allele_2		The other allele of the genotype.
@column sample_id		Foreign key references to the @link individual table.

@see individual
@see variation
@see subsnp_handle
*/

create table individual_genotype_multiple_bp (
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
@table meta_coord

@colour #DA70D6
@desc This table gives the coordinate system used by various tables in the database.

@column table_name			Name of the feature table, e.g. "variation_feature".
@column coord_system_id		Foreign key to core database coord_system table refers to coordinate system that features from this table can be found in.
@column max_length			Maximun length of the feature. 
*/

CREATE TABLE meta_coord (

  table_name                  VARCHAR(40) NOT NULL,
  coord_system_id             INT(10) UNSIGNED NOT NULL,
  max_length		      INT,

  UNIQUE(table_name, coord_system_id)

) TYPE=MyISAM;


/**
@table meta

@colour #DA70D6
@desc This table stores various metadata relating to the database, generally used by the Ensembl web code.

@column meta_id		Primary key, internal identifier.
@column species_id	...
@column meta_key		Name of the meta entry, e.g. "schema_version".
@column meta_value	Corresponding value of the key, e.g. "61".
*/

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


/**
@table tagged_variation_feature

@desc This table lists variation features that are tagged by another variation feature. Tag pairs are defined as having an r<sup>2</sup> &gt; 0.99.

@column variation_feature_id	Primary key. Foreign key references to the @link variation_feature table.
@column sample_id					Primary key. Foreign key references to the @link sample table.

@see variation_feature
@see population
*/

CREATE TABLE tagged_variation_feature (

  variation_feature_id       INT(10) UNSIGNED not null,
  sample_id              INT(10) UNSIGNED not null,
  
  PRIMARY KEY(variation_feature_id, sample_id)
);


/**
@table read_coverage

@colour #FF8500
@desc This table stores the read coverage in the resequencing of individuals. Each row contains an individual ID, chromosomal coordinates and a read coverage level.

@column seq_region_id		Foreign key references @link seq_region in core db. ers to the seq_region which this variant is on, which may be a chromosome, a clone, etc...
@column seq_region_start	The start position of the variation on the @link seq_region.
@column seq_region_end		The end position of the variation on the @link seq_region.
@column level					Minimum number of reads.
@column sample_id				Foreign key references to the @link individual table.

@see individual
@see seq_region
*/

CREATE TABLE read_coverage (
   seq_region_id int(10) unsigned not null,
   seq_region_start int not null,
   seq_region_end int not null,
   level tinyint not null,
   sample_id int(10) unsigned not null,
		  
   key seq_region_idx(seq_region_id,seq_region_start)   
);

/**
@table compressed_genotype_single_bp

@colour #FF8500
@desc This table holds genotypes compressed using the pack() method in Perl. These genotypes are mapped to particular genomic locations rather than variation objects. The data have been compressed to reduce table size and increase the speed of the web code.

@column sample_id				Primary key. Foreign key references to the sample table.
@column seq_region_id		Foreign key references @link seq_region in core db. ers to the seq_region which this variant is on, which may be a chromosome, a clone, etc...
@column seq_region_start	The start position of the variation on the @link seq_region.
@column seq_region_end		The end position of the variation on the @link seq_region.
@column seq_region_strand	The orientation of the variation on the @link seq_region.
@column genotypes				Encoded representation of the genotype data:<br />Each row in the compressed table stores genotypes from one individual in one fixed-size region of the genome (arbitrarily defined as 100 Kb). The compressed string (using Perl's pack method) consisting of a repeating triplet of elements: a distance in base pairs from the previous genotype followed by a pair of alleles.<br />For example, a given row may have a start position of 1000, indicating the chromosomal position of the first genotype in this row. The unpacked genotypes field then may contain the following elements:<br />0, A, G, 20, C, C, 35, G, T, 320, A, A, ...<br />The first genotype has a position of 1000 + 0 = 1000 and alleles A and G.<br />The second genotype has a position of 1000 + 20 = 1020 and alleles C and C.<br />The third genotype similarly has a position of 1055 and alleles G and T, and so on.

@see individual
@see seq_region
*/

CREATE TABLE compressed_genotype_single_bp(
  sample_id int(10) unsigned not null,
  seq_region_id int(10) unsigned not null,
  seq_region_start int not null,
  seq_region_end int not null,
  seq_region_strand tinyint not null,
  genotypes blob,

  key pos_idx(seq_region_id,seq_region_start)
) MAX_ROWS = 100000000;


/**
@table failed_description

@colour #3CB371
@desc This table contains descriptions of reasons for a variation being flagged as failed.

@column failed_description_id	Primary key, internal identifier.
@column description				Text containing the reason why the Variation has been flagged as failed. e.g. "Variation does not map to the genome".

@see failed_variation
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

/**
@table  attrib_type

@colour #FF0000
@desc   Defines the set of possible attribute types used in the attrib table

@column attrib_type_id  Primary key
@column code            A short codename for this type (indexed, so should be used for lookups)
@column name            The name of this type
@column description     Longer description of this type

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

@see    attrib_type
@see    attrib_set
*/
CREATE TABLE attrib (

    attrib_id           INT(11) UNSIGNED NOT NULL DEFAULT 0,
    attrib_type_id      SMALLINT(5) UNSIGNED NOT NULL DEFAULT 0,
    value               TEXT NOT NULL,

    PRIMARY KEY (attrib_id),
    UNIQUE KEY type_val_idx (attrib_type_id, value(40))
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
@table  protein_function_predictions

@colour #1E90FF
@desc   Contains encoded sift and polyphen predictions for every protein-coding transcript in this species

@column translation_stable_id   The stable ID of the translation
@column transcript_stable_id    The stable ID of the transcript from which this protein is translated
@column translation_md5         A hexidecimal string representing the MD5 hash of the protein sequence
@column polyphen_predictions    A compressed binary string containing the PolyPhen-2 predictions for all
                                possible amino acid substitutions in this translation
@column sift_predictions        A similarly encoded string for SIFT predictions
*/

CREATE TABLE protein_function_predictions (

    translation_stable_id   VARCHAR(128) NOT NULL,
    transcript_stable_id    VARCHAR(128) NOT NULL,
    translation_md5         CHAR(32) NOT NULL,
    polyphen_predictions    MEDIUMBLOB,
    sift_predictions        MEDIUMBLOB,
    
    PRIMARY KEY (translation_stable_id),
    KEY transcript_idx (transcript_stable_id)
);

#possible values in the failed_description table

INSERT INTO failed_description (failed_description_id,description) VALUES (1,'Variation maps to more than 3 different locations');
INSERT INTO failed_description (failed_description_id,description) VALUES (2,'None of the variant alleles match the reference allele');
INSERT INTO failed_description (failed_description_id,description) VALUES (3,'Variation has more than 3 different alleles');
INSERT INTO failed_description (failed_description_id,description) VALUES (4,'Loci with no observed variant alleles in dbSNP');
INSERT INTO failed_description (failed_description_id,description) VALUES (5,'Variation does not map to the genome');
INSERT INTO failed_description (failed_description_id,description) VALUES (6,'Variation has no genotypes');
INSERT INTO failed_description (failed_description_id,description) VALUES (7,'Genotype frequencies do not add up to 1');
INSERT INTO failed_description (failed_description_id,description) VALUES (8,'Variation has no associated sequence');
INSERT INTO failed_description (failed_description_id,description) VALUES (9,'Variation submission has been withdrawn by the 1000 genomes project due to high false positive rate');
INSERT INTO failed_description (failed_description_id,description) VALUES (11,'Allele string obtained from dbSNP does not agree with the submitted alleles'); 
INSERT INTO failed_description (failed_description_id,description) VALUES (12,'Variation has more than 3 different submitted alleles');         
INSERT INTO failed_description (failed_description_id,description) VALUES (13,'Alleles contain non-nucleotide characters');  
INSERT INTO failed_description (failed_description_id,description) VALUES (14,'Alleles contain ambiguity codes');  
INSERT INTO failed_description (failed_description_id,description) VALUES (15,'Mapped position is not compatible with reported alleles');


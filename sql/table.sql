create table variation(
	variation_id int not null auto_increment,
	source_id int not null,
	name varchar(255),
	parent_variation_id int,

	primary key( variation_id ),
	key parent_var_idx( parent_variation_id ),
	key name_idx( name, source_id )
);

create table allele(
	allele_id int not null auto_increment,
	variation_id int not null,
	allele text,
	frequency float,
	population_id int,

	primary key( allele_id ),
	key variation_idx( variation_id )
);

create table population(
	population_id int not null auto_increment,
	name varchar(255) not null,
	method varchar(255) not null,
	parent_population_id int,
	size int,

	primary key( population_id ),
	key parent_pop_idx( parent_population_id ) 
);

create table variation_feature(
	variation_feature_id int not null auto_increment,
	seq_region_id int not null,
	seq_region_start int not null,
	seq_region_end int not null,
	seq_region_strand enum( "-1", "0", "1" ) default "0"  not null,
	variation_id int not null,
	allele_string text,
	method varchar(255),
	variation_name varchar(255),
	map_weight int not null,

	primary key( variation_feature_id ),
	key pos_idx( seq_region_id, seq_region_start, method ),
	key variation_idx( variation_id )
);


create table transcript_variation(
	transcript_variation_id int not null auto_increment,
	variation_feature_id int not null,
	amino_acid_change varchar(255),
	amino_acid_position int,
	cdna_position int,
	type enum( "INTRONIC", "5PRIME", "3PRIME", "SYNONYMOUS_CODING",
	"NON_SYNONYMOUS_CODING", "5PRIME_UTR", "3PRIME_UTR" ),
	
	primary key( transcript_variation_id ),
	key variation_idx( variation_feature_id )
	);

create table allele_group(
	allele_group_id int not null auto_increment,
	population_id int,
	name varchar(255),
	frequency float,
	type enum( "HAPLOBLOCK", "GENOTYPE" ),

	primary key( allele_group_id )

);

create table allele_group_allele (
	allele_id int not null,
	allele_group_id int not null,

	unique( allele_group_id, allele_id ),
	key allele_idx( allele_id, allele_group_id )
);

create table flanking_sequence (
	variation_id int not null,
	upstream varchar(255),
	downstream varchar(255),

	primary key( variation_id )
);


create table httag(
	httag_id int not null auto_increment,
	variation_id int not null,
	name varchar(255),
	source_id int not null,

	primary key( httag_id ),
	key variation_idx( variation_id )
);

create table source(
	source_id int not null auto_increment,
	name varchar(255),
	
	primary key( source_id )
);

## Replace sample table with more comprehensive versions of individual and population table

## Rename individual/population table to old_individual/old_population
## Create new individual/population table with new columns from sample table  
## Populate individual/population table from old_individual/old_population and sample tables
## Drop old_individual and old_population


RENAME TABLE individual to old_individual;

CREATE TABLE individual(
    individual_id int(10) unsigned not null auto_increment,
    name varchar(255),
    description text,
    gender enum('Male', 'Female', 'Unknown') default 'Unknown' NOT NULL,
    father_individual_id int(10) unsigned,
    mother_individual_id int(10) unsigned,
    individual_type_id int(10) unsigned NOT NULL DEFAULT 0,
    display enum('REFERENCE', 'DEFAULT', 'DISPLAYABLE', 'UNDISPLAYABLE', 'LD', 'MARTDISPLAYABLE') default 'UNDISPLAYABLE',

    primary key(individual_id)
);

INSERT INTO individual(individual_id, name, description, gender, father_individual_id, mother_individual_id, individual_type_id, display)
SELECT oi.sample_id, s.name, s.description, oi.gender, oi.father_individual_sample_id, oi.mother_individual_sample_id, oi.individual_type_id, s.display
FROM old_individual oi, sample s
WHERE oi.sample_id = s.sample_id;

## Now for population
RENAME TABLE population to old_population;

CREATE TABLE population(
    population_id int(10) unsigned not null auto_increment,
    name varchar(255),
    size int(10),
    description text,
    collection tinyint(1) default 0,
    freqs_from_gts tinyint(1),
    display enum('LD', 'MARTDISPLAYABLE', 'UNDISPLAYABLE') default 'UNDISPLAYABLE',

    primary key(population_id)
);

INSERT INTO population(population_id, name, size, description, collection, freqs_from_gts, display)
SELECT op.sample_id, s.name, s.size, s.description, 0, s.freqs_from_gts, s.display
FROM old_population op, sample s
WHERE op.sample_id = s.sample_id;

DROP TABLE old_individual;
DROP TABLE old_population;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_71_72_c.sql|Move data from sample table to new individual and population tables.');

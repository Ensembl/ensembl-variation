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


##Rename column sample_id which was relating to old sample table to either individual_id or population_id
##Rename indexes 

DELIMITER $$

DROP PROCEDURE IF EXISTS drop_index $$
CREATE PROCEDURE drop_index(in TableName varchar(128), in IndexName varchar(128))
BEGIN
    IF((SELECT COUNT(*) FROM information_schema.statistics WHERE TABLE_SCHEMA = DATABASE() and table_name = TableName AND index_name = IndexName)  > 0) THEN
        SET @s = CONCAT('DROP INDEX ', IndexName, ' ON ', TableName);
        PREPARE stmt FROM @s;
        EXECUTE stmt;
    END IF;
END $$

DROP PROCEDURE IF EXISTS change_column_name $$
CREATE PROCEDURE change_column_name(in TableName varchar(128), in OldColumnName varchar(128), in NewColumnName varchar(128), in DataType varchar(128))
BEGIN
    IF ((SELECT COUNT(*) FROM information_schema.tables WHERE TABLE_SCHEMA = DATABASE() and table_name = TableName) = 1) THEN
        IF((SELECT COUNT(*) FROM information_schema.columns WHERE TABLE_SCHEMA = DATABASE() and table_name = TableName AND column_name = NewColumnName) = 0) THEN
            SET @s = CONCAT('ALTER TABLE ', TableName, ' CHANGE ', OldColumnName, ' ', NewColumnName, ' ', DataType);
            PREPARE stmt FROM @s;
            EXECUTE stmt;
        END IF;
    END IF;
END $$

DROP PROCEDURE IF EXISTS create_index $$
CREATE PROCEDURE create_index(in TableName varchar(128), in IndexName varchar(128), in ColumnName varchar(128), in IndexType varchar(128))
BEGIN
    IF ((SELECT COUNT(*) FROM information_schema.tables WHERE TABLE_SCHEMA = DATABASE() and table_name = TableName) = 1) THEN   
        IF ((SELECT COUNT(*) FROM information_schema.statistics WHERE TABLE_SCHEMA = DATABASE() and table_name = TableName AND index_name = IndexName)  = 0) THEN
            SET @s = CONCAT('CREATE ', IndexType, ' INDEX ', IndexName, ' ON ', TableName, '(', ColumnName, ')');
            PREPARE stmt FROM @s;
            EXECUTE stmt;
        END IF;
    END IF;
END $$

CALL drop_index('allele', 'sample_idx');
CALL drop_index('individual_genotype_multiple_bp', 'sample_idx');
CALL drop_index('individual_population', 'individual_sample_idx');
CALL drop_index('individual_population', 'population_sample_idx');
CALL drop_index('population_genotype', 'sample_idx');
CALL drop_index('population_structure', 'super_population_sample_id');
CALL drop_index('population_structure', 'sub_pop_sample_idx');
CALL drop_index('tagged_variation_feature', 'sample_idx');
CALL drop_index('tmp_individuale_genotypes_single_bp', 'sample_idx');

CALL change_column_name('allele', 'sample_id', 'population_id', 'INT');
CALL change_column_name('individual_genotype_multiple_bp', 'sample_id', 'individual_id', 'INT UNSIGNED');
CALL change_column_name('individual_population', 'individual_sample_id', 'individual_id', 'INT UNSIGNED NOT NULL');
CALL change_column_name('individual_population', 'population_sample_id', 'population_id', 'INT UNSIGNED NOT NULL');
CALL change_column_name('population_genotype', 'sample_id', 'population_id', 'INT UNSIGNED');
CALL change_column_name('population_structure', 'super_population_sample_id', 'super_population_id', 'INT UNSIGNED NOT NULL');
CALL change_column_name('population_structure', 'sub_population_sample_id', 'sub_population_id', 'INT UNSIGNED NOT NULL');
CALL change_column_name('read_coverage', 'sample_id', 'individual_id', 'INT UNSIGNED NOT NULL');
CALL change_column_name('tagged_variation_feature', 'sample_id', 'population_id', 'INT UNSIGNED NOT NULL');
CALL change_column_name('tmp_individual_genotype_single_bp', 'sample_id', 'individual_id', 'INT UNSIGNED');

CALL create_index('allele', 'population_idx', 'population_id', '');
CALL create_index('individual_genotype_multiple_bp', 'individual_idx', 'individual_id', '');
CALL create_index('individual_population', 'individual_idx', 'individual_id', '');
CALL create_index('individual_population', 'population_idx', 'population_id', '');
CALL create_index('population_genotype', 'population_idx', 'population_id', '');
CALL create_index('population_structure', 'super_population_idx', 'super_population_id, sub_population_id', 'UNIQUE');
CALL create_index('population_structure', 'sub_population_idx', 'sub_population_id', '');
CALL create_index('tagged_variation_feature', 'population_idx', 'population_id', '');
CALL create_index('tmp_individual_genotype_single_bp', 'individual_idx', 'individual_id', '');

DROP PROCEDURE IF EXISTS drop_index;
DROP PROCEDURE IF EXISTS change_column_name;
DROP PROCEDURE IF EXISTS create_index;

##compressed_genotype_region: copy and rename is much faster than alter table..
CREATE TABLE compressed_genotype_region_tmp LIKE compressed_genotype_region;
DROP INDEX sample_idx ON compressed_genotype_region_tmp;
ALTER TABLE compressed_genotype_region_tmp CHANGE sample_id individual_id INT UNSIGNED NOT NULL;
CREATE INDEX individual_idx ON compressed_genotype_region_tmp(individual_id);
ALTER TABLE compressed_genotype_region_tmp DISABLE KEYS;
INSERT INTO compressed_genotype_region_tmp(individual_id,seq_region_id,seq_region_start,seq_region_end,seq_region_strand,genotypes)
SELECT sample_id,seq_region_id,seq_region_start,seq_region_end,seq_region_strand,genotypes FROM compressed_genotype_region;
ALTER TABLE compressed_genotype_region_tmp ENABLE KEYS;
RENAME TABLE compressed_genotype_region TO compressed_genotype_region_old;
RENAME TABLE compressed_genotype_region_tmp TO compressed_genotype_region;
DROP TABLE compressed_genotype_region_old;

##patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_71_72_b.sql|Changes for sample table redesign: Rename columns and indexes');

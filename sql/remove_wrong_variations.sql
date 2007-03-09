#set of SQL statements to remove from the database variations that have been flagged as wrong
#by the Ensembl team, and to update the Variation table to reflect the reason
DELETE a FROM allele a, failed_variation v WHERE a.variation_id = v.variation_id;
DELETE f FROM flanking_sequence f, failed_variation v WHERE f.variation_id = v.variation_id; 
DELETE i FROM individual_genotype_multiple_bp i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE tv FROM transcript_variation tv, variation_feature vf, failed_variation v where tv.variation_feature_id = vf.variation_feature_id and vf.variation_id = v.variation_id;
DELETE vs FROM variation_synonym vs, failed_variation v WHERE vs.variation_id = v.variation_id;
DELETE i FROM tmp_individual_genotype_single_bp i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE i FROM individual_genotype_multiple_bp i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE i FROM population_genotype i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE tv FROM tagged_variation_feature tv, variation_feature vf, failed_variation v where tv.variation_feature_id = vf.variation_feature_id and vf.variation_id = v.variation_id;
DELETE vf FROM variation_feature vf, failed_variation v WHERE vf.variation_id = v.variation_id;
SELECT 'Remember to run again compressed_genotype script to remove variations from the compressed table' as '';
TRUNCATE TABLE compressed_genotype_single_bp; #emtpy the compressed table for future recreation
DELETE FROM meta_coord where table_name = 'compressed_genotype_single_bp';

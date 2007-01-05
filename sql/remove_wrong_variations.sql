#set of SQL statements to remove from the database variations that have been flagged as wrong
#by the Ensembl team, and to update the Variation table to reflect the reason
DELETE a FROM allele a, failed_variation v WHERE a.variation_id = v.variation_id;
DELETE f FROM flanking_sequence f, failed_variation v WHERE f.variation_id = v.variation_id; 
DELETE i FROM individual_genotype_multiple_bp i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE vf FROM variation_feature vf, failed_variation v WHERE vf.variation_id = v.variation_id;
DELETE vs FROM variation_synonym vs, failed_variation v WHERE vs.variation_id = v.variation_id;
DELETE i FROM tmp_individual_genotype_single_bp i, failed_variation v WHERE i.variation_id = v.variation_id;
SELECT 'Remember to run again compressed_genotype script to remove variations from the compressed table' as '';
#after removing suspicious variations, update the Variation table to reflect changes
UPDATE variation v, failed_variation w set v.validation_status = 'failed', v.ancestral_allele = '', v.failed_description_id = w.failed_description_id WHERE v.variation_id = w.variation_id;
#and finally, remove table with the variations
DROP TABLE failed_variation;

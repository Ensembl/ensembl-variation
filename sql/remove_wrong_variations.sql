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


#set of SQL statements to remove from the database variations that have been flagged as wrong
#by the Ensembl team, and to update the Variation table to reflect the reason
#start delete here normally
#
DELETE fv FROM failed_variation fv, variation v WHERE v.variation_id = fv.variation_id AND v.validation_status LIKE '%precious%';
#DELETE fv FROM failed_variation fv, variation_feature vf WHERE fv.variation_id = vf.variation_id AND vf.map_weight != 1;
DELETE fv FROM failed_variation fv, variation_set vst, variation_set_variation vsv WHERE fv.variation_id = vsv.variation_id AND vst.variation_set_id = vsv.variation_set_id AND vst.name = 'Clinical/LSDB variations from dbSNP' ;
UPDATE variation v, failed_variation w set v.validation_status = 'failed', v.ancestral_allele = NULL WHERE v.variation_id = w.variation_id;
DELETE a FROM allele a, failed_variation v WHERE a.variation_id = v.variation_id;
DELETE f FROM flanking_sequence f, failed_variation v WHERE f.variation_id = v.variation_id; 
DELETE i FROM individual_genotype_multiple_bp i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE tv FROM transcript_variation tv, variation_feature vf, failed_variation v where tv.variation_feature_id = vf.variation_feature_id and vf.variation_id = v.variation_id;
DELETE vs FROM variation_synonym vs, failed_variation v WHERE vs.variation_id = v.variation_id;
DELETE i FROM tmp_individual_genotype_single_bp i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE i FROM population_genotype i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE tv FROM tagged_variation_feature tv, variation_feature vf, failed_variation v where tv.variation_feature_id = vf.variation_feature_id and vf.variation_id = v.variation_id;
DELETE vf FROM variation_feature vf, failed_variation v WHERE vf.variation_id = v.variation_id;
DELETE vf FROM variation_set_variation vf, failed_variation v WHERE vf.variation_id = v.variation_id;
SELECT 'Remember to run again compressed_genotype script to remove variations from the compressed table' as '';
#TRUNCATE TABLE compressed_genotype_single_bp; #emtpy the compressed table for future recreation
#DELETE FROM meta_coord where table_name = 'compressed_genotype_single_bp';

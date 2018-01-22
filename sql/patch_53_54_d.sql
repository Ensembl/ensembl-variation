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


# patch_53_54_d.sql
#
# title: add no_mappings in failed_description table
#
# description:
# If a variation don't have a entry in variation_feature table, this means no mappings for the variation
# need to change validation_status='failed' and remove it from rest of tables

INSERT INTO failed_description (failed_description_id,description) VALUES (5,'Variation do not have genome mappings');
#then, copy the information from the variation, variation_feature table
INSERT IGNORE INTO failed_variation (variation_id,failed_description_id) 
   SELECT v.variation_id, 5
   FROM  variation v LEFT JOIN variation_feature vf
   ON v.variation_id=vf.variation_id
   WHERE vf.variation_id is null;
   #AND v.validation_status != 'failed';

#delete newly failed_variation from other tables
UPDATE variation v, failed_variation w set v.validation_status = 'failed', v.ancestral_allele = NULL WHERE v.variation_id = w.variation_id;
DELETE a FROM allele a, failed_variation v WHERE a.variation_id = v.variation_id;
DELETE f FROM flanking_sequence f, failed_variation v WHERE f.variation_id = v.variation_id; 
DELETE i FROM individual_genotype_multiple_bp i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE tv FROM transcript_variation tv, variation_feature vf, failed_variation v where tv.variation_feature_id = vf.variation_feature_id and vf.variation_id = v.variation_id;
DELETE vs FROM variation_synonym vs, failed_variation v WHERE vs.variation_id = v.variation_id;
DELETE i FROM tmp_individual_genotype_single_bp i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE i FROM individual_genotype_multiple_bp i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE i FROM population_genotype i, failed_variation v WHERE i.variation_id = v.variation_id;
DELETE tv FROM tagged_variation_feature tv, variation_feature vf, failed_variation v where tv.variation_feature_id = vf.variation_feature_id and vf.variation_id = v.variation_id;

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_53_54_d.sql|add no_mappings in failed_description table');

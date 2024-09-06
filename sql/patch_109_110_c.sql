-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

# Add new clinical_significance
ALTER TABLE variation MODIFY clinical_significance 
  SET('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective','affects','likely pathogenic low penetrance','pathogenic low penetrance','uncertain risk allele','likely risk allele','established risk allele')
  DEFAULT NULL;

# Add new clinical_significance
ALTER TABLE variation_feature MODIFY clinical_significance 
  SET('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective','affects','likely pathogenic low penetrance','pathogenic low penetrance','uncertain risk allele','likely risk allele','established risk allele')
  DEFAULT NULL;

# Add new clinical_significance
ALTER TABLE structural_variation MODIFY clinical_significance
  SET('uncertain significance','not provided','benign','likely benign','likely pathogenic','pathogenic','drug response','histocompatibility','other','confers sensitivity','risk factor','association','protective','affects','likely pathogenic low penetrance','pathogenic low penetrance','uncertain risk allele','likely risk allele','established risk allele')
  DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_109_110_c.sql|Add new clinical_significance values to variation, variation_feature and structural_variation');

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


# add HGMD_MUTATION consequence to variation_feature
ALTER TABLE `variation_feature` CHANGE `consequence_type` `consequence_type` set('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','COMPLEX_INDEL','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','PARTIAL_CODON','SYNONYMOUS_CODING','REGULATORY_REGION','WITHIN_MATURE_miRNA','5PRIME_UTR','3PRIME_UTR','INTRONIC','NMD_TRANSCRIPT','UPSTREAM','DOWNSTREAM','WITHIN_NON_CODING_GENE','NO_CONSEQUENCE','INTERGENIC','HGMD_MUTATION') NOT NULL DEFAULT 'INTERGENIC';

# add HGMD_MUTATION consequence to transcript_variation
ALTER TABLE `transcript_variation` CHANGE `consequence_type` `consequence_type` set('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','COMPLEX_INDEL','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','PARTIAL_CODON','SYNONYMOUS_CODING','REGULATORY_REGION','WITHIN_MATURE_miRNA','5PRIME_UTR','3PRIME_UTR','INTRONIC','NMD_TRANSCRIPT','UPSTREAM','DOWNSTREAM','WITHIN_NON_CODING_GENE','NO_CONSEQUENCE','INTERGENIC','HGMD_MUTATION') NOT NULL DEFAULT 'INTERGENIC';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_57_58_e.sql|new consequence type HGMD_MUTATION');

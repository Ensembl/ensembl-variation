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


#this patch is suposed to change the order in the set column

ALTER TABLE variation_feature ADD column new_consequence_type SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM','INTERGENIC') DEFAULT 'INTERGENIC' not null;
UPDATE variation_feature SET new_consequence_type = consequence_type;
ALTER TABLE variation_feature DROP COLUMN consequence_type, CHANGE new_consequence_type consequence_type  SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM','INTERGENIC') DEFAULT 'INTERGENIC' not null;

#do the same in the trancript_variation table

ALTER TABLE transcript_variation ADD column new_consequence_type SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM') not null;
UPDATE transcript_variation SET new_consequence_type = consequence_type;
ALTER TABLE transcript_variation DROP COLUMN consequence_type, CHANGE new_consequence_type consequence_type  SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM') not null;

INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_39_40_b.sql|modify_set_column');

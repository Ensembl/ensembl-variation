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


# patch_51_52_b.sql
#
# title: change consequence_type column in variation_feature 
# and transcript_variation tables to add two columns :
# 'WITHIN_MATURE_miRNA' and 'WITHIN_NON_CODING_GENE'
#


ALTER table variation_feature CHANGE consequence_type consequence_type SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','COMPLEX_INDEL','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION','WITHIN_MATURE_miRNA','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM','WITHIN_NON_CODING_GENE','INTERGENIC')
        default "INTERGENIC" not null ; 

ALTER table transcript_variation CHANGE consequence_type consequence_type SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','COMPLEX_INDEL','FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION','WITHIN_MATURE_miRNA','5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM','WITHIN_NON_CODING_GENE') not null;

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_51_52_d.sql|change consequence_type');

alter table seq_region drop KEY name_cs_idx ;
alter table seq_region add UNIQUE KEY name_idx(name);

delete from meta where meta_value = "patch_51_52
_b.sql|add seq_region table";

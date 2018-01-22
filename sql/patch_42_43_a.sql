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


#this patch changes database schema version from 42 -> 43


#first, add a new column in the Variation table

ALTER TABLE variation ADD failed_description_id int(10) unsigned not null, MODIFY validation_status SET('cluster','freq','submitter','doublehit','hapmap','failed');

#
# create the failed_description table
#
# Contains reasons for removing some variations from the Variation database
#
# failed_description_id  - primary key, internal identifier
# description - text containing the reason why the Variation information has been removed from the 
#               Variation databse except in the Variation table
#

CREATE TABLE failed_description(

 failed_description_id int(10) unsigned not null,
 description  text not null,

 PRIMARY KEY (failed_description_id)
);

#possible values in the failed_description table
INSERT INTO failed_description (failed_description_id,description) VALUES (1,'Variation has more than 3 different locations');
INSERT INTO failed_description (failed_description_id,description) VALUES (2,'Reference allele not present in the alleles of the variation');
INSERT INTO failed_description (failed_description_id,description) VALUES (3,'Variation containing more than 3 alleles');
INSERT INTO failed_description (failed_description_id,description) VALUES (4,'Variation with \'NoVariation\' alleles');


#update the VariationFeature and the TranscriptVariation tables with the new consequenceType
ALTER TABLE variation_feature MODIFY consequence_type SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','COMPLEX_INDEL',
			'FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION',
			'5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM','INTERGENIC')
	default "INTERGENIC" not null;

ALTER TABLE transcript_variation MODIFY consequence_type SET ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','COMPLEX_INDEL',
			'FRAMESHIFT_CODING','NON_SYNONYMOUS_CODING','SPLICE_SITE','SYNONYMOUS_CODING','REGULATORY_REGION',
			'5PRIME_UTR','3PRIME_UTR','INTRONIC','UPSTREAM','DOWNSTREAM') not null;

#and update the meta table with the patch applied
INSERT INTO meta (meta_key,meta_value) VALUES ('patch','patch_42_43_a.sql|update database schema');	

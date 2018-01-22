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


# structural variation schema changes

DROP TABLE supporting_structural_variation;


ALTER TABLE structural_variation ADD COLUMN is_evidence TINYINT(4) DEFAULT 0;

ALTER TABLE structural_variation ADD INDEX source_idx (source_id);
ALTER TABLE structural_variation_feature ADD INDEX source_idx (source_id);

CREATE TABLE structural_variation_association (
  structural_variation_id int(10) unsigned NOT NULL,
  supporting_structural_variation_id int(10) unsigned NOT NULL,
	
  PRIMARY KEY (structural_variation_id, supporting_structural_variation_id)
);


CREATE TABLE structural_variation_annotation (
	structural_variation_annotation_id int(10) unsigned NOT NULL auto_increment,
	structural_variation_id int(10) unsigned NOT NULL,
	clinical_attrib_id int(10) unsigned DEFAULT NULL,
	phenotype_id int(10) unsigned DEFAULT NULL,
	sample_id int(10) unsigned DEFAULT NULL,
	strain_id int(10) unsigned DEFAULT NULL,
	
	primary key (structural_variation_annotation_id),
	key structural_variation_idx(structural_variation_id),
	key clinical_attrib_idx(clinical_attrib_id),
	key phenotype_idx(phenotype_id),
	key sample_idx(sample_id),
	key strain_idx(strain_id)
);


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_64_65_b.sql|structural variation schema changes');

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


# structural variation changes

create table structural_variation_feature (
	structural_variation_feature_id int(10) unsigned NOT NULL AUTO_INCREMENT,
	seq_region_id int(10) unsigned NOT NULL,
	outer_start int,	
	seq_region_start int NOT NULL,
	inner_start int,
	inner_end int,
	seq_region_end int NOT NULL,
	outer_end int,
	seq_region_strand tinyint NOT NULL,
	structural_variation_id int(10) unsigned NOT NULL,
  variation_name varchar(255),
	source_id int(10) unsigned NOT NULL, 
  class_attrib_id int(10) unsigned NOT NULL DEFAULT 0,
	allele_string longtext DEFAULT NULL,
	
  PRIMARY KEY (structural_variation_feature_id),
	KEY pos_idx( seq_region_id, seq_region_start, seq_region_end ),
	KEY structural_variation_idx (structural_variation_id),
	KEY attrib_idx (class_attrib_id)
);

ALTER TABLE supporting_structural_variation ADD COLUMN class_attrib_id int(10) unsigned NOT NULL DEFAULT 0;
ALTER TABLE supporting_structural_variation ADD INDEX attrib_idx (class_attrib_id);

ALTER TABLE structural_variation ADD COLUMN somatic tinyint(1) NOT NULL DEFAULT 0;

ALTER TABLE structural_variation DROP COLUMN seq_region_id;
ALTER TABLE structural_variation DROP COLUMN seq_region_start;
ALTER TABLE structural_variation DROP COLUMN seq_region_end;
ALTER TABLE structural_variation DROP COLUMN seq_region_strand;
ALTER TABLE structural_variation DROP COLUMN inner_start;
ALTER TABLE structural_variation DROP COLUMN inner_end;
ALTER TABLE structural_variation DROP COLUMN allele_string;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_63_64_b.sql|structural variation changes');

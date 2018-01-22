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


# add MAF columns to variation_feature

# create tmp table and add columns
CREATE TABLE tmp_vf LIKE variation_feature;
ALTER TABLE tmp_vf ADD COLUMN `minor_allele` char(1) DEFAULT NULL;
ALTER TABLE tmp_vf ADD COLUMN `minor_allele_freq` float DEFAULT NULL;
ALTER TABLE tmp_vf ADD COLUMN `minor_allele_count` int(10) unsigned DEFAULT NULL;

# copy data with join to variation
INSERT INTO tmp_vf
SELECT vf.*, v.minor_allele, v.minor_allele_freq, v.minor_allele_count
FROM variation v, variation_feature vf
WHERE v.variation_id = vf.variation_id;

# drop and rename
DROP TABLE variation_feature;
RENAME TABLE tmp_vf TO variation_feature;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_67_68_e.sql|add MAF columns to variation_feature');
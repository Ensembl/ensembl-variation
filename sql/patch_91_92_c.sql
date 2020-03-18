-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

# Change the column zygosity in the table structural_variation_sample 

ALTER TABLE structural_variation_sample CHANGE zygosity zygosity_bak ENUM ('homozygous', 'heterozygous') DEFAULT NULL;

ALTER TABLE structural_variation_sample ADD COLUMN zygosity TINYINT(1) DEFAULT NULL;

UPDATE structural_variation_sample SET zygosity=1 WHERE zygosity_bak='heterozygous';
UPDATE structural_variation_sample SET zygosity=2 WHERE zygosity_bak='homozygous';

ALTER TABLE structural_variation_sample DROP COLUMN zygosity_bak;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_91_92_c.sql|Change the column zygosity in the table structural_variation_sample');


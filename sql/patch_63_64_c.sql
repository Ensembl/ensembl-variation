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


# update to the new protein function prediction schema

DROP TABLE protein_info;
DROP TABLE protein_position;
DROP TABLE polyphen_prediction;
DROP TABLE sift_prediction;

CREATE TABLE protein_function_predictions (

    translation_stable_id   VARCHAR(128) NOT NULL,
    transcript_stable_id    VARCHAR(128) NOT NULL,
    translation_md5         CHAR(32) NOT NULL,
    polyphen_predictions    MEDIUMBLOB,
    sift_predictions        MEDIUMBLOB,
    
    PRIMARY KEY (translation_stable_id),
    KEY transcript_idx (transcript_stable_id)
);

# add sift and polyphen scores to the transcript_variation table for biomart

ALTER TABLE transcript_variation ADD COLUMN polyphen_score FLOAT DEFAULT NULL; 
ALTER TABLE transcript_variation ADD COLUMN sift_score FLOAT DEFAULT NULL; 

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_63_64_c.sql|update to new protein function prediction schema');

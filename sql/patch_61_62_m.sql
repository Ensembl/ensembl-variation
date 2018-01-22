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


# add new schema for nsSNP predictions 

CREATE TABLE protein_info (
    protein_info_id         int(10) unsigned NOT NULL AUTO_INCREMENT,
    transcript_stable_id    varchar(128) NOT NULL,
    transcript_version      smallint unsigned NOT NULL,
    translation_md5         char(32) NOT NULL,

    PRIMARY KEY         (protein_info_id),
    KEY transcript_idx  (transcript_stable_id, transcript_version)
);

#
# protein_position
#
# Table with a row for each position in every ensembl translation, used by the
# nsSNP prediction tables
#

CREATE TABLE protein_position (
    protein_position_id             int(10) unsigned NOT NULL AUTO_INCREMENT,
    protein_info_id                 int(10) unsigned NOT NULL,
    position                        mediumint unsigned NOT NULL,
    amino_acid                      char(1) NOT NULL,
    sift_median_conservation        float,
    sift_num_sequences_represented  smallint,

    PRIMARY KEY (protein_position_id),
    KEY pos_idx (protein_info_id, position)
);

#
# polyphen_prediction
#
# Table storing the polyphen prediction for every possible amino acid
# substitution in the ensembl proteome
#

CREATE TABLE polyphen_prediction (
    polyphen_prediction_id  int(10) unsigned NOT NULL AUTO_INCREMENT,
    protein_position_id     int(10) unsigned NOT NULL,
    amino_acid              char(1) NOT NULL,
    prediction              enum('unknown', 'benign', 'possibly damaging', 'probably damaging') NOT NULL,
    probability             float NOT NULL,
    compressed_result_hash  blob,
    
    PRIMARY KEY     (polyphen_prediction_id),
    KEY pos_aa_idx  (protein_position_id, amino_acid)
);

#
# sift_prediction
#
# Table storing the sift prediction for every possible amino acid
# substitution in the ensembl proteome
#

CREATE table sift_prediction (
    sift_prediction_id      int(10) unsigned NOT NULL AUTO_INCREMENT,
    protein_position_id     int(10) NOT NULL,
    amino_acid              char(1) NOT NULL,
    prediction              enum('tolerated', 'deleterious') NOT NULL,
    score                   float NOT NULL,

    PRIMARY KEY     (sift_prediction_id),
    KEY pos_aa_idx  (protein_position_id, amino_acid)
);


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_m.sql|add new schema for nsSNP predictions');


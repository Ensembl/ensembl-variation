-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2025] EMBL-European Bioinformatics Institute
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

/**
Use MyISAM storage engine
*/
SET default_storage_engine=MYISAM;

# Links variation to PMID (PubMed identifier) from dbSNP
# Used to populate tables publication and variation_citation
# on completion of load
CREATE TABLE tmp_variation_citation (
  variation_id  INT UNSIGNED NOT NULL,
  pmid          INT UNSIGNED NOT NULL,
  PRIMARY KEY (variation_id,pmid)
);

# Stores SPDI information used to determine variation_feature.allele_string
CREATE TABLE placement_allele (
  placement_allele_id INT UNSIGNED NOT NULL AUTO_INCREMENT, # PK
  variation_id INT UNSIGNED,
  seq_id VARCHAR(30),
  position INT,
  deleted_sequence VARCHAR(255),
  inserted_sequence VARCHAR(255),
  hgvs VARCHAR(255),
  PRIMARY KEY (placement_allele_id),
  KEY (variation_id)
);

# Stores information on batch load
CREATE TABLE batch (
  batch_id INT UNSIGNED NOT NULL AUTO_INCREMENT,
  filename VARCHAR(255) NOT NULL,
  parent_filename VARCHAR(255),
  current VARCHAR(3) DEFAULT 'yes',
  PRIMARY KEY (batch_id),
  UNIQUE KEY file_ndx(filename, parent_filename)
);

# Links variation to batch
CREATE TABLE batch_variation (
  batch_id INT UNSIGNED,
  variation_id INTEGER UNSIGNED NOT NULL,
  variant_type VARCHAR(30),
  PRIMARY KEY (batch_id, variation_id),
  KEY (variation_id)
);

# Stores descriptions of reasons that SPDI cannot be converted
# to an allele_string
CREATE TABLE spdi_failed_description (
  spdi_failed_description_id INT UNSIGNED NOT NULL AUTO_INCREMENT,
  description text NOT NULL,
  PRIMARY KEY (spdi_failed_description_id)
);

INSERT INTO spdi_failed_description (spdi_failed_description_id, description) VALUES (1,'No SPDI for reference allele');
INSERT INTO spdi_failed_description (spdi_failed_description_id, description) VALUES (2,'SPDIs deleted_sequence differ');
INSERT INTO spdi_failed_description (spdi_failed_description_id, description) VALUES (3,'SPDIs position differ');
INSERT INTO spdi_failed_description (spdi_failed_description_id, description) VALUES (4,'No variation in SPDI');

# Stores information about a variation_feature that had failures for converting 
# SPDI to allele_string 
CREATE TABLE failed_variation_feature_spdi (
  variation_feature_id int(10) unsigned NOT NULL,
  spdi_failed_description_id int(10) unsigned NOT NULL,
  PRIMARY KEY (variation_feature_id, spdi_failed_description_id)
);


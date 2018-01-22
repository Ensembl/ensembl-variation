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


## copy data to phenotype_feature and phenotype_feature_attrib

# ensure we have the necessary attrib_type entries
INSERT IGNORE INTO `attrib_type` (`attrib_type_id`, `code`, `name`, `description`)
VALUES
	(13, 'associated_gene', 'Associated gene', 'ID of gene(s) linked by a phenotype association'),
	(14, 'risk_allele', 'Risk allele', 'Risk allele in phenotype association'),
	(15, 'p_value', 'P-value', 'P-value denoting significance of an observed phenotype annotation'),
	(16, 'variation_names', 'Variation names', 'ID of variant(s) linked with a phenotype association'),
	(17, 'sample_id', 'Sample ID', 'Sample ID for source of phenotype association'),
    (18, 'strain_id', 'Strain ID', 'Strain ID for source of phenotype association');
    

# temporarily alter phenotype_feature
ALTER TABLE `phenotype_feature` ADD `variation_annotation_id` INT  NULL  DEFAULT NULL  AFTER `seq_region_strand`;
ALTER TABLE `phenotype_feature` ADD INDEX (`variation_annotation_id`);

# populate phenotype_feature from variation_annotation
INSERT INTO phenotype_feature(phenotype_id, source_id, study_id, object_id, type, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, variation_annotation_id)
SELECT phenotype_id, s.source_id, va.study_id, vf.variation_name, 'Variation', vf.seq_region_id, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, va.variation_annotation_id
FROM variation_annotation va, variation_feature vf, study s
WHERE va.variation_id = vf.variation_id
AND va.study_id = s.study_id;

# populate phenotype_feature_attrib
INSERT INTO phenotype_feature_attrib
SELECT p.phenotype_feature_id, t.attrib_type_id, va.associated_gene
FROM variation_annotation va, phenotype_feature p, attrib_type t
WHERE va.variation_annotation_id = p.variation_annotation_id
AND t.code = "associated_gene"
AND va.associated_gene IS NOT NULL;

INSERT INTO phenotype_feature_attrib
SELECT p.phenotype_feature_id, t.attrib_type_id, va.associated_variant_risk_allele
FROM variation_annotation va, phenotype_feature p, attrib_type t
WHERE va.variation_annotation_id = p.variation_annotation_id
AND t.code = "risk_allele"
AND va.associated_variant_risk_allele IS NOT NULL;

INSERT INTO phenotype_feature_attrib
SELECT p.phenotype_feature_id, t.attrib_type_id, cast(va.p_value as char)
FROM variation_annotation va, phenotype_feature p, attrib_type t
WHERE va.variation_annotation_id = p.variation_annotation_id
AND t.code = "p_value"
AND va.p_value IS NOT NULL;

INSERT INTO phenotype_feature_attrib
SELECT p.phenotype_feature_id, t.attrib_type_id, va.variation_names
FROM variation_annotation va, phenotype_feature p, attrib_type t
WHERE va.variation_annotation_id = p.variation_annotation_id
AND t.code = "variation_names"
AND va.variation_names IS NOT NULL;


## structural_variation_annotation
INSERT INTO phenotype_feature(phenotype_id, object_id, type, is_significant, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, variation_annotation_id)
SELECT phenotype_id, vf.variation_name, if(vf.is_evidence = 0, 'StructuralVariation', 'SupportingStructuralVariation'), if(phenotype_id IS NULL, 0, 1), vf.seq_region_id, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, va.structural_variation_annotation_id
FROM structural_variation_annotation va, structural_variation_feature vf
WHERE va.structural_variation_id = vf.structural_variation_id;

# populate phenotype_feature_attrib
INSERT INTO phenotype_feature_attrib
SELECT p.phenotype_feature_id, t.attrib_type_id, a.value
FROM structural_variation_annotation va, phenotype_feature p, attrib a, attrib_type t
WHERE va.structural_variation_annotation_id = p.variation_annotation_id
AND va.clinical_attrib_id = a.attrib_id
AND t.code = "dgva_clin_sig";

INSERT INTO phenotype_feature_attrib
SELECT p.phenotype_feature_id, t.attrib_type_id, va.sample_id
FROM structural_variation_annotation va, phenotype_feature p, attrib_type t
WHERE va.structural_variation_annotation_id = p.variation_annotation_id
AND va.sample_id IS NOT NULL
AND t.code = "sample_id";

INSERT INTO phenotype_feature_attrib
SELECT p.phenotype_feature_id, t.attrib_type_id, va.strain_id
FROM structural_variation_annotation va, phenotype_feature p, attrib_type t
WHERE va.structural_variation_annotation_id = p.variation_annotation_id
AND va.strain_id IS NOT NULL
AND t.code = "strain_id";

# revert phenotype_feature
ALTER TABLE `phenotype_feature` DROP INDEX `variation_annotation_id`;
ALTER TABLE `phenotype_feature` DROP COLUMN `variation_annotation_id`;


# update meta_coord
INSERT IGNORE INTO meta_coord
SELECT 'phenotype_feature', a.*, b.*
FROM (
    SELECT DISTINCT(coord_system_id)
    FROM meta_coord
) AS a,
(
    SELECT max(seq_region_end - seq_region_start)
    FROM phenotype_feature
    WHERE seq_region_start <= seq_region_end
) AS b;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_h.sql|copy data to phenotype_feature and phenotype_feature_attrib');

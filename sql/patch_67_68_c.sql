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


## change consequence_type columns in transcript_variation and variation_feature

## TRANSCRIPT VARIATION
# create tmp table
CREATE TABLE tmp_tv LIKE transcript_variation;

# alter its consequence type field
ALTER TABLE tmp_tv CHANGE consequence_types consequence_types SET('splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','nc_transcript_variant','non_coding_exon_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','initiator_codon_variant','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation') DEFAULT NULL;

# insert into tmp table, replacing old terms with new as we go
INSERT INTO tmp_tv SELECT transcript_variation_id, variation_feature_id, feature_stable_id, allele_string, somatic,
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(consequence_types,
        'non_synonymous_codon', 'missense_variant'),
        'synonymous_codon', 'synonymous_variant'),
        'initiator_codon_change', 'initiator_codon_variant'),
        'inframe_codon_loss', 'inframe_deletion'),
        'inframe_codon_gain', 'inframe_insertion'),
        '5KB_downstream_variant', 'downstream_gene_variant'),
        '500B_downstream_variant', 'downstream_gene_variant'),
        '5KB_upstream_variant', 'upstream_gene_variant'),
        '2KB_upstream_variant', 'upstream_gene_variant'),
cds_start, cds_end, cdna_start, cdna_end, translation_start, translation_end, codon_allele_string, pep_allele_string, hgvs_genomic, hgvs_transcript, hgvs_protein, polyphen_prediction, polyphen_score, sift_prediction, sift_score
FROM transcript_variation;

DROP TABLE transcript_variation;
RENAME TABLE tmp_tv TO transcript_variation;

## VARIATION_FEATURE
# create tmp table
CREATE TABLE tmp_vf LIKE variation_feature;

# alter its consequence type field
ALTER TABLE tmp_vf CHANGE consequence_type consequence_type SET('intergenic_variant','splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','nc_transcript_variant','non_coding_exon_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','initiator_codon_variant','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation') NOT NULL DEFAULT 'intergenic_variant';

# insert into tmp table, replacing old terms with new as we go
INSERT INTO tmp_vf SELECT variation_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, variation_id, allele_string, variation_name, map_weight, flags, source_id, validation_status,
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(
    REPLACE(consequence_type,
        'non_synonymous_codon', 'missense_variant'),
        'synonymous_codon', 'synonymous_variant'),
        'initiator_codon_change', 'initiator_codon_variant'),
        'inframe_codon_loss', 'inframe_deletion'),
        'inframe_codon_gain', 'inframe_insertion'),
        '5KB_downstream_variant', 'downstream_gene_variant'),
        '500B_downstream_variant', 'downstream_gene_variant'),
        '5KB_upstream_variant', 'upstream_gene_variant'),
        '2KB_upstream_variant', 'upstream_gene_variant'),
variation_set_id, class_attrib_id, somatic
FROM variation_feature;

DROP TABLE variation_feature;
RENAME TABLE tmp_vf TO variation_feature;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_67_68_c.sql|change consequence_type sets in transcript_variation and variation_feature');

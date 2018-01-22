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


# change consequence_types column in motif_feature_variation, regulatory_feature_variation, transcript_variation, variation_feature
#    - nc_transcript_variant will be updated to non_coding_transcript_variant
#    - non_coding_exon_variant will be updated to non_coding_transcript_exon_variant


## motif_feature_variation
# create tmp table
CREATE TABLE tmp_patch_mfv LIKE motif_feature_variation;

# alter its consequence type field
ALTER TABLE tmp_patch_mfv CHANGE consequence_types consequence_types SET('splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','non_coding_transcript_variant','non_coding_transcript_exon_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','initiator_codon_variant','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TF_binding_site_variant','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_variant','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation') DEFAULT NULL;

# insert into tmp table, replacing old terms with new as we go
INSERT INTO tmp_patch_mfv SELECT motif_feature_variation_id, variation_feature_id, feature_stable_id, motif_feature_id, allele_string, somatic,
    REPLACE(
    REPLACE(consequence_types,
        'nc_transcript_variant', 'non_coding_transcript_variant'),
        'non_coding_exon_variant', 'non_coding_transcript_exon_variant'),
motif_name, motif_start, motif_end, motif_score_delta, in_informative_position
FROM motif_feature_variation;

DROP TABLE motif_feature_variation;
RENAME TABLE tmp_patch_mfv TO motif_feature_variation;
## -------------------------------------------------------------------------------

## regulatory_feature_variation
# create tmp table
CREATE TABLE tmp_patch_rfv LIKE regulatory_feature_variation;

# alter its consequence type field
ALTER TABLE tmp_patch_rfv CHANGE consequence_types consequence_types SET('splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','non_coding_transcript_variant','non_coding_transcript_exon_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','initiator_codon_variant','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TF_binding_site_variant','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_variant','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation') DEFAULT NULL;

# insert into tmp table, replacing old terms with new as we go
INSERT INTO tmp_patch_rfv SELECT regulatory_feature_variation_id, variation_feature_id, feature_stable_id, feature_type, allele_string, somatic,
    REPLACE(
    REPLACE(consequence_types,
        'nc_transcript_variant', 'non_coding_transcript_variant'),
        'non_coding_exon_variant', 'non_coding_transcript_exon_variant')
FROM regulatory_feature_variation;

DROP TABLE regulatory_feature_variation;
RENAME TABLE tmp_patch_rfv TO regulatory_feature_variation;
## -------------------------------------------------------------------------------

## transcript_variation
# create tmp table
CREATE TABLE tmp_patch_tv LIKE transcript_variation;

# alter its consequence type field
ALTER TABLE tmp_patch_tv CHANGE consequence_types consequence_types SET('splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','non_coding_transcript_variant','non_coding_transcript_exon_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','initiator_codon_variant','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation') DEFAULT NULL;

# insert into tmp table, replacing old terms with new as we go
INSERT INTO tmp_patch_tv SELECT transcript_variation_id, variation_feature_id, feature_stable_id, allele_string, somatic,
    REPLACE(
    REPLACE(consequence_types,
        'nc_transcript_variant', 'non_coding_transcript_variant'),
        'non_coding_exon_variant', 'non_coding_transcript_exon_variant'),
cds_start, cds_end, cdna_start, cdna_end, translation_start, translation_end, distance_to_transcript, codon_allele_string, pep_allele_string, hgvs_genomic, hgvs_transcript, hgvs_protein, polyphen_prediction, polyphen_score, sift_prediction, sift_score
FROM transcript_variation;

DROP TABLE transcript_variation;
RENAME TABLE tmp_patch_tv TO transcript_variation;
## -------------------------------------------------------------------------------

## variation_feature
# create tmp table
CREATE TABLE tmp_patch_vf LIKE variation_feature;

# alter its consequence type field
ALTER TABLE tmp_patch_vf CHANGE consequence_types consequence_types SET('intergenic_variant','splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','non_coding_transcript_variant','non_coding_transcript_exon_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','initiator_codon_variant','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation','regulatory_region_variant','TF_binding_site_variant') NOT NULL DEFAULT 'intergenic_variant';

# insert into tmp table, replacing old terms with new as we go
INSERT INTO tmp_patch_vf SELECT variation_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, variation_id, allele_string, variation_name, map_weight, flags, source_id, validation_status,
    REPLACE(
    REPLACE(consequence_types,
        'nc_transcript_variant', 'non_coding_transcript_variant'),
        'non_coding_exon_variant', 'non_coding_transcript_exon_variant'),
variation_set_id, class_attrib_id, somatic, minor_allele, minor_allele_freq, minor_allele_count, alignment_quality, clinical_significance, evidence_attribs
FROM variation_feature;

DROP TABLE variation_feature;
RENAME TABLE tmp_patch_vf TO variation_feature;
## -------------------------------------------------------------------------------


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_76_77_c.sql|update SO consequence terms');


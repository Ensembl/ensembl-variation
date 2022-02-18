-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

## Remove TFBS and regulatory consequences from transcript_variation table

## TRANSCRIPT VARIATION
ALTER TABLE transcript_variation CHANGE consequence_types consequence_types SET (
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_lost',
    'coding_sequence_variant',
    'missense_variant',
    'stop_gained',
    'synonymous_variant',
    'frameshift_variant',
    'non_coding_transcript_variant',
    'non_coding_transcript_exon_variant',
    'mature_miRNA_variant',
    'NMD_transcript_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'incomplete_terminal_codon_variant',
    'intron_variant',
    'splice_region_variant',
    'downstream_gene_variant',
    'upstream_gene_variant',
    'start_lost',
    'stop_retained_variant',
    'inframe_insertion',
    'inframe_deletion',
    'transcript_ablation',
    'transcript_fusion',
    'transcript_amplification',
    'transcript_translocation',
    'feature_elongation',
    'feature_truncation',
    'protein_altering_variant',
    'start_retained_variant',
    'splice_donor_5th_base_variant',
    'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant'
);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_106_107_c.sql|remove regulatory and TFBS consequences from consequence_types in transcript_variation');

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


# change the consequence_type column in variation_feature to use SO accessions

ALTER TABLE variation_feature MODIFY consequence_type SET (
    'intergenic_variant',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'complex_change_in_transcript', 
    'stop_lost',
    'coding_sequence_variant',
    'non_synonymous_codon',
    'stop_gained',
    'synonymous_codon',
    'frameshift_variant',
    'nc_transcript_variant',
    'mature_miRNA_variant',
    'NMD_transcript_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'incomplete_terminal_codon_variant',
    'intron_variant',
    'splice_region_variant',
    '5KB_downstream_variant',
    '500B_downstream_variant',
    '5KB_upstream_variant',
    '2KB_upstream_variant',
    'initiator_codon_change',
    'stop_retained_variant',
    'inframe_codon_gain',
    'inframe_codon_loss',
    'miRNA_target_site_variant',
    'pre_miRNA_variant',
    'regulatory_region_variant',
    'increased_binding_affinity',
    'decreased_binding_affinity',
    'binding_site_variant'
) DEFAULT 'intergenic_variant' NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_i.sql|change the consequence_type column in variation_feature to use SO terms');

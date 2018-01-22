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


# create new regulatory region variation tables
CREATE TABLE IF NOT EXISTS `motif_feature_variation` (
    `motif_feature_variation_id`          int(11) unsigned NOT NULL AUTO_INCREMENT,
    `variation_feature_id`                int(11) unsigned NOT NULL,
    `feature_stable_id`                   varchar(128) DEFAULT NULL,
    `motif_feature_id`                    int(11) unsigned NOT NULL,
    `allele_string`                       text,
    `somatic`                             tinyint(1) NOT NULL DEFAULT 0,
    `consequence_types`                   set (
                                            'splice_acceptor_variant',
                                            'splice_donor_variant',
                                            'stop_lost',
                                            'coding_sequence_variant',
                                            'missense_variant',
                                            'stop_gained',
                                            'synonymous_variant',
                                            'frameshift_variant',
                                            'nc_transcript_variant',
                                            'non_coding_exon_variant',
                                            'mature_miRNA_variant',
                                            'NMD_transcript_variant',
                                            '5_prime_UTR_variant',
                                            '3_prime_UTR_variant',
                                            'incomplete_terminal_codon_variant',
                                            'intron_variant',
                                            'splice_region_variant',
                                            'downstream_gene_variant',
                                            'upstream_gene_variant',
                                            'initiator_codon_variant',
                                            'stop_retained_variant',
                                            'inframe_insertion',
                                            'inframe_deletion',
                                            'transcript_ablation',
                                            'transcript_fusion',
                                            'transcript_amplification',
                                            'transcript_translocation',
                                            'TF_binding_site_variant',
                                            'TFBS_ablation',
                                            'TFBS_fusion',
                                            'TFBS_amplification',
                                            'TFBS_translocation',
                                            'regulatory_region_variant',
                                            'regulatory_region_ablation',
                                            'regulatory_region_fusion',
                                            'regulatory_region_amplification',
                                            'regulatory_region_translocation',
                                            'feature_elongation',
                                            'feature_truncation'
                                        ),
    `motif_name`                          text,
    `motif_start`                         int(11) unsigned,
    `motif_end`                           int(11) unsigned,
    `motif_score_delta`                   float DEFAULT NULL,
    `in_informative_position`             tinyint(1) NOT NULL DEFAULT 0,

    PRIMARY KEY                           (`motif_feature_variation_id`),
    KEY `variation_feature_idx`           (`variation_feature_id`),
    KEY `feature_idx`                     (`feature_stable_id`),
    KEY `consequence_type_idx`            (`consequence_types`),
    KEY `somatic_feature_idx`             (`feature_stable_id`, `somatic`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS `regulatory_feature_variation` (
    `regulatory_feature_variation_id`     int(11) unsigned NOT NULL AUTO_INCREMENT,
    `variation_feature_id`                int(11) unsigned NOT NULL,
    `feature_stable_id`                   varchar(128) DEFAULT NULL,
    `feature_type`                        text, 
    `allele_string`                       text,
    `somatic`                             tinyint(1) NOT NULL DEFAULT 0,
    `consequence_types`                   set (
                                            'splice_acceptor_variant',
                                            'splice_donor_variant',
                                            'stop_lost',
                                            'coding_sequence_variant',
                                            'missense_variant',
                                            'stop_gained',
                                            'synonymous_variant',
                                            'frameshift_variant',
                                            'nc_transcript_variant',
                                            'non_coding_exon_variant',
                                            'mature_miRNA_variant',
                                            'NMD_transcript_variant',
                                            '5_prime_UTR_variant',
                                            '3_prime_UTR_variant',
                                            'incomplete_terminal_codon_variant',
                                            'intron_variant',
                                            'splice_region_variant',
                                            'downstream_gene_variant',
                                            'upstream_gene_variant',
                                            'initiator_codon_variant',
                                            'stop_retained_variant',
                                            'inframe_insertion',
                                            'inframe_deletion',
                                            'transcript_ablation',
                                            'transcript_fusion',
                                            'transcript_amplification',
                                            'transcript_translocation',
                                            'TF_binding_site_variant',
                                            'TFBS_ablation',
                                            'TFBS_fusion',
                                            'TFBS_amplification',
                                            'TFBS_translocation',
                                            'regulatory_region_variant',
                                            'regulatory_region_ablation',
                                            'regulatory_region_fusion',
                                            'regulatory_region_amplification',
                                            'regulatory_region_translocation',
                                            'feature_elongation',
                                            'feature_truncation'
                                        ),

    PRIMARY KEY                           (`regulatory_feature_variation_id`),
    KEY `variation_feature_idx`           (`variation_feature_id`),
    KEY `feature_idx`                     (`feature_stable_id`),
    KEY `consequence_type_idx`            (`consequence_types`),
    KEY `somatic_feature_idx`             (`feature_stable_id`, `somatic`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_69_70_d.sql|add regulatory region variation tables');

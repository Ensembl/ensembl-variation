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


# update to new transcript_variation schema

DROP TABLE IF EXISTS transcript_variation;

CREATE TABLE attrib_type (

    attrib_type_id    SMALLINT(5) UNSIGNED NOT NULL DEFAULT 0,
    code              VARCHAR(20) NOT NULL DEFAULT '',
    name              VARCHAR(255) NOT NULL DEFAULT '',
    description       TEXT,

    PRIMARY KEY (attrib_type_id),
    UNIQUE KEY code_idx (code)
);

CREATE TABLE attrib (

    attrib_id           INT(11) UNSIGNED NOT NULL DEFAULT 0,
    attrib_type_id      SMALLINT(5) UNSIGNED NOT NULL DEFAULT 0,
    value               TEXT NOT NULL,

    PRIMARY KEY (attrib_id),
    KEY type_val_idx (attrib_type_id, value(40))
);

CREATE TABLE attrib_set (

    attrib_set_id       INT(11) UNSIGNED NOT NULL DEFAULT 0,
    attrib_id           INT(11) UNSIGNED NOT NULL DEFAULT 0,

    UNIQUE KEY set_idx (attrib_set_id, attrib_id),
    KEY attrib_idx (attrib_id)
);


CREATE TABLE transcript_variation (
    transcript_variation_id             int(11) unsigned NOT NULL AUTO_INCREMENT,
    variation_feature_id                int(11) unsigned NOT NULL,
    feature_stable_id                   varchar(128) DEFAULT NULL,
    allele_string                       text,
    somatic                             tinyint(1) NOT NULL DEFAULT 0,
    consequence_types                   set (
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
                                            'pre_miRNA_variant'
                                        ),
    cds_start                           int(11) unsigned,
    cds_end                             int(11) unsigned,
    cdna_start                          int(11) unsigned,
    cdna_end                            int(11) unsigned,
    translation_start                   int(11) unsigned,
    translation_end                     int(11) unsigned,
    codon_allele_string                 text,
    pep_allele_string                   text,
    hgvs_genomic                        text,
    hgvs_coding                         text,
    hgvs_protein                        text,
    polyphen_prediction                 enum('unknown', 'benign', 'possibly damaging', 'probably damaging') DEFAULT NULL,
    sift_prediction                     enum('tolerated', 'deleterious') DEFAULT NULL,

    PRIMARY KEY                         (transcript_variation_id),
    KEY variation_feature_idx           (variation_feature_id),
    KEY feature_idx                     (feature_stable_id),
    KEY consequence_type_idx            (consequence_types),
    KEY somatic_feature_idx             (feature_stable_id, somatic)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_l.sql|update to new transcript_variation schema');

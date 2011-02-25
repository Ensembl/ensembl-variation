# update to new transcript_variation schema

DROP TABLE IF EXISTS transcript_variation;

CREATE TABLE attrib_type (

    attrib_type_id    SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
    code              VARCHAR(20) NOT NULL DEFAULT '',
    name              VARCHAR(255) NOT NULL DEFAULT '',
    description       TEXT,

    PRIMARY KEY (attrib_type_id),
    UNIQUE KEY code_idx (code)
);

CREATE TABLE attrib (

    attrib_id           INT(11) UNSIGNED NOT NULL AUTO_INCREMENT,
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
                                            'miRNA_target_site_variant',
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

SET 
    @SO_accession        = 1, 
    @SO_term             = 2, 
    @display_term        = 3, 
    @ncbi_term           = 4, 
    @feature_SO_term     = 5, 
    @conseq_rank         = 6, 
    @conseq_predicate    = 7, 
    @ens_feature_class   = 8, 
    @ens_feature_subtype = 9, 
    @ens_variant_class   = 10
;

INSERT INTO attrib_type 
    (attrib_type_id, code, description)
VALUES 
    (@SO_accession,         'SO_accession',         'Sequence Ontology accession'),
    (@SO_term,              'SO_term',              'Sequence Ontology term'),
    (@display_term,         'display_term',         'Ensembl display term'),
    (@ncbi_term,            'ncbi_term',            'NCBI term'),
    (@feature_SO_term,      'feature_SO_term',      'Sequence Ontology term for the associated feature'),
    (@conseq_rank,          'conseq_rank',          'Relative severity of this variation consequence'),
    (@conseq_predicate,     'conseq_predicate',     'Predicate function which tests an associated variation consequence'),
    (@ens_feature_class,    'ens_feature_class',    'Associated Ensembl feature class'),
    (@ens_feature_subtype,  'ens_feature_subtype',  'Subtype of the associated Ensembl feature'),
    (@ens_variant_class,    'ens_variant_class',    'Associated Ensembl Variation feature class')
;

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (1, @SO_accession, 'SO:0001628');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (2, @SO_term, 'intergenic_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (3, @display_term, 'INTERGENIC');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (4, @conseq_rank, '100');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (1, 1);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (1, 2);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (1, 3);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (1, 4);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (5, @SO_accession, 'SO:0001635');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (6, @SO_term, '5KB_upstream_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (7, @display_term, 'UPSTREAM');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (8, @feature_SO_term, 'transcript');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (9, @conseq_rank, '20');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (10, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream_5KB');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (2, 5);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (2, 6);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (2, 7);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (2, 8);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (2, 9);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (2, 10);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (11, @SO_accession, 'SO:0001633');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (12, @SO_term, '5KB_downstream_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (13, @display_term, 'DOWNSTREAM');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (14, @conseq_rank, '21');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (15, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream_5KB');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (3, 11);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (3, 12);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (3, 13);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (3, 8);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (3, 14);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (3, 15);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (16, @SO_accession, 'SO:0001636');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (17, @SO_term, '2KB_upstream_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (18, @ncbi_term, 'near-gene-5');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (19, @conseq_rank, '18');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (20, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream_2KB');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (4, 16);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (4, 17);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (4, 7);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (4, 18);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (4, 8);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (4, 19);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (4, 20);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (21, @SO_accession, 'SO:0001634');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (22, @SO_term, '500B_downstream_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (23, @ncbi_term, 'near-gene-3');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (24, @conseq_rank, '19');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (25, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream_500B');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (5, 21);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (5, 22);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (5, 13);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (5, 23);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (5, 8);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (5, 24);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (5, 25);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (26, @SO_accession, 'SO:0001575');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (27, @SO_term, 'splice_donor_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (28, @display_term, 'ESSENTIAL_SPLICE_SITE');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (29, @ncbi_term, 'splice-5');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (30, @feature_SO_term, 'primary_transcript');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (31, @conseq_rank, '1');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (32, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::donor_splice_site');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (6, 26);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (6, 27);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (6, 28);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (6, 29);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (6, 30);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (6, 31);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (6, 32);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (33, @SO_accession, 'splice_acceptor_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (34, @SO_term, 'splice_acceptor_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (35, @ncbi_term, 'splice-3');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (36, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::acceptor_splice_site');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (7, 33);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (7, 34);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (7, 28);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (7, 35);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (7, 30);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (7, 31);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (7, 36);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (37, @SO_accession, 'SO:0001630');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (38, @SO_term, 'splice_region_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (39, @display_term, 'SPLICE_SITE');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (40, @conseq_rank, '8');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (41, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::splice_region');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (8, 37);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (8, 38);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (8, 39);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (8, 30);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (8, 40);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (8, 41);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (42, @SO_accession, 'SO:0001627');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (43, @SO_term, 'intron_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (44, @display_term, 'INTRONIC');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (45, @ncbi_term, 'intron');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (46, @conseq_rank, '15');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (47, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_intron');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (9, 42);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (9, 43);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (9, 44);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (9, 45);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (9, 30);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (9, 46);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (9, 47);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (48, @SO_accession, 'SO:0001623');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (49, @SO_term, '5_prime_UTR_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (50, @display_term, '5PRIME_UTR');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (51, @ncbi_term, 'untranslated_5');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (52, @feature_SO_term, 'mRNA');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (53, @conseq_rank, '13');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (54, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_5_prime_utr');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (10, 48);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (10, 49);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (10, 50);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (10, 51);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (10, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (10, 53);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (10, 54);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (55, @SO_accession, 'SO:0001624');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (56, @SO_term, '3_prime_UTR_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (57, @display_term, '3PRIME_UTR');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (58, @ncbi_term, 'untranslated_3');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (59, @conseq_rank, '14');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (60, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_3_prime_utr');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (11, 55);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (11, 56);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (11, 57);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (11, 58);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (11, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (11, 59);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (11, 60);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (61, @SO_accession, 'SO:0001577');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (62, @SO_term, 'complex_change_in_transcript');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (63, @display_term, 'COMPLEX_INDEL');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (64, @conseq_rank, '5');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (65, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::complex_indel');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (12, 61);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (12, 62);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (12, 63);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (12, 30);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (12, 64);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (12, 65);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (66, @SO_accession, 'SO:0001588');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (67, @SO_term, 'synonymous_codon');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (68, @display_term, 'SYNONYMOUS_CODING');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (69, @ncbi_term, 'cds-synon');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (70, @conseq_rank, '10');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (71, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::synonymous_codon');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (13, 66);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (13, 67);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (13, 68);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (13, 69);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (13, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (13, 70);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (13, 71);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (72, @SO_term, 'non_synonymous_codon');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (73, @display_term, 'NON_SYNONYMOUS_CODING');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (74, @ncbi_term, 'missense');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (75, @conseq_rank, '7');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (76, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::non_synonymous_codon');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (14, 66);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (14, 72);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (14, 73);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (14, 74);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (14, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (14, 75);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (14, 76);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (77, @SO_accession, 'SO:0001651');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (78, @SO_term, 'inframe_codon_gain');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (79, @conseq_rank, '6');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (80, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_codon_gain');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (15, 77);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (15, 78);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (15, 73);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (15, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (15, 79);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (15, 80);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (81, @SO_accession, 'SO:0001652');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (82, @SO_term, 'inframe_codon_loss');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (83, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_codon_loss');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (16, 81);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (16, 82);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (16, 73);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (16, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (16, 79);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (16, 83);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (84, @SO_accession, 'SO:0001587');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (85, @SO_term, 'stop_gained');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (86, @display_term, 'STOP_GAINED');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (87, @ncbi_term, 'nonsense');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (88, @conseq_rank, '3');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (89, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_gained');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (17, 84);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (17, 85);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (17, 86);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (17, 87);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (17, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (17, 88);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (17, 89);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (90, @SO_accession, 'SO:0001578');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (91, @SO_term, 'stop_lost');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (92, @display_term, 'STOP_LOST');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (93, @conseq_rank, '4');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (94, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_lost');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (18, 90);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (18, 91);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (18, 92);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (18, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (18, 93);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (18, 94);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (95, @SO_accession, 'SO:0001567');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (96, @SO_term, 'stop_retained_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (97, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_retained');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (19, 95);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (19, 96);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (19, 68);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (19, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (19, 70);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (19, 97);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (98, @SO_accession, 'SO:0001582');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (99, @SO_term, 'initiator_codon_change');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (100, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::affects_start_codon');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (20, 98);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (20, 99);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (20, 73);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (20, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (20, 75);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (20, 100);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (101, @SO_accession, 'SO:0001589');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (102, @SO_term, 'frameshift_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (103, @display_term, 'FRAMESHIFT_CODING');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (104, @ncbi_term, 'frameshift');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (105, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::frameshift');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (21, 101);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (21, 102);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (21, 103);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (21, 104);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (21, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (21, 79);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (21, 105);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (106, @SO_accession, 'SO:0001626');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (107, @SO_term, 'incomplete_terminal_codon_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (108, @display_term, 'PARTIAL_CODON');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (109, @conseq_rank, '9');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (110, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::partial_codon');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (22, 106);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (22, 107);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (22, 108);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (22, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (22, 109);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (22, 110);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (111, @SO_accession, 'SO:0001621');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (112, @SO_term, 'NMD_transcript_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (113, @display_term, 'NMD_TRANSCRIPT');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (114, @conseq_rank, '16');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (115, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_nmd_transcript');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (23, 111);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (23, 112);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (23, 113);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (23, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (23, 114);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (23, 115);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (116, @SO_accession, 'SO:0001619');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (117, @SO_term, 'nc_transcript_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (118, @display_term, 'WITHIN_NON_CODING_GENE');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (119, @feature_SO_term, 'ncRNA');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (120, @conseq_rank, '17');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (121, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_non_coding_gene');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (24, 116);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (24, 117);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (24, 118);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (24, 119);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (24, 120);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (24, 121);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (122, @SO_accession, 'SO:0001620');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (123, @SO_term, 'mature_miRNA_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (124, @display_term, 'WITHIN_MATURE_miRNA');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (125, @feature_SO_term, 'miRNA');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (126, @conseq_rank, '12');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (127, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_mature_miRNA');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (25, 122);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (25, 123);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (25, 124);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (25, 125);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (25, 126);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (25, 127);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (128, @SO_accession, 'SO:0001580');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (129, @SO_term, 'coding_sequence_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (130, @display_term, 'CODING_UNKNOWN');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (131, @conseq_rank, '11');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (132, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::coding_unknown');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (26, 128);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (26, 129);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (26, 130);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (26, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (26, 131);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (26, 132);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (133, @SO_accession, 'SO:0001566');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (134, @SO_term, 'regulatory_region_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (135, @display_term, 'REGULATORY_REGION');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (136, @feature_SO_term, 'regulatory_region');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (137, @conseq_rank, '50');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (138, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_regulatory_feature');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (27, 133);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (27, 134);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (27, 135);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (27, 136);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (27, 137);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (27, 138);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (139, @SO_accession, 'SO:X000005');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (140, @SO_term, 'pre_miRNA_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (141, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_miRNA');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (28, 139);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (28, 140);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (28, 118);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (28, 125);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (28, 53);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (28, 141);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (142, @SO_accession, 'SO:X000004');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (143, @SO_term, 'miRNA_target_site_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (144, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_miRNA_target_site');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (29, 142);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (29, 143);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (29, 135);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (29, 8);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (29, 53);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (29, 144);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (145, @SO_accession, 'SO:X000003');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (146, @SO_term, 'binding_site_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (147, @feature_SO_term, 'binding_site');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (148, @conseq_rank, '49');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (149, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_motif_feature');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (30, 145);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (30, 146);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (30, 135);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (30, 147);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (30, 148);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (30, 149);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (150, @SO_accession, 'SO:X000002');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (151, @SO_term, 'decreased_binding_affinity');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (152, @conseq_rank, '47');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (153, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::decreased_binding_affinity');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (31, 150);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (31, 151);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (31, 135);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (31, 147);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (31, 152);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (31, 153);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (154, @SO_accession, 'SO:X000001');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (155, @SO_term, 'increased_binding_affinity');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (156, @conseq_rank, '48');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (157, @conseq_predicate, 'Bio::EnsEMBL::Variation::Utils::VariationEffect::increased_binding_affinity');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (32, 154);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (32, 155);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (32, 135);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (32, 147);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (32, 156);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (32, 157);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (158, @SO_accession, 'SO:0001483');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (159, @SO_term, 'SNV');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (160, @display_term, 'SNP');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (33, 158);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (33, 159);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (33, 160);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (161, @SO_accession, 'SO:1000002');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (162, @SO_term, 'substitution');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (34, 161);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (34, 162);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (163, @SO_accession, 'SO:0001019');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (164, @SO_term, 'copy_number_variation');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (165, @display_term, 'CNV');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (35, 163);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (35, 164);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (35, 165);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (166, @SO_accession, 'SO:0000667');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (167, @SO_term, 'insertion');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (36, 166);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (36, 167);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (168, @SO_accession, 'SO:0000159');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (169, @SO_term, 'deletion');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (37, 168);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (37, 169);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (170, @SO_accession, 'SO:1000032');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (171, @SO_term, 'indel');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (38, 170);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (38, 171);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (172, @SO_accession, 'SO:0000705');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (173, @SO_term, 'tandem_repeat');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (39, 172);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (39, 173);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (174, @SO_accession, 'SO:0001059');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (175, @SO_term, 'sequence_alteration');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (40, 174);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (40, 175);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (176, @SO_accession, 'SO:0000234');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (177, @SO_term, 'mRNA');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (178, @ens_feature_class, 'Bio::EnsEMBL::Transcript');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (179, @ens_feature_subtype, 'protein_coding');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (180, @ens_variant_class, 'Bio::EnsEMBL::Variation::TranscriptVariationNew');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (41, 176);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (41, 177);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (41, 178);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (41, 179);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (41, 180);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (181, @SO_accession, 'SO:0000673');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (182, @SO_term, 'transcript');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (42, 181);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (42, 182);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (42, 178);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (42, 180);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (183, @SO_accession, 'SO:0000185');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (184, @SO_term, 'primary_transcript');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (43, 183);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (43, 184);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (43, 178);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (43, 180);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (185, @SO_accession, 'SO:0000655');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (186, @SO_term, 'ncRNA');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (44, 185);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (44, 186);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (44, 178);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (44, 180);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (187, @SO_accession, 'SO:0000276');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (188, @SO_term, 'miRNA');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (45, 187);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (45, 188);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (45, 178);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (45, 180);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (189, @SO_accession, 'SO:0005836');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (190, @SO_term, 'regulatory_region');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (191, @ens_feature_class, 'Bio::EnsEMBL::Funcgen::RegulatoryFeature');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (192, @ens_variant_class, 'Bio::EnsEMBL::Variation::RegulatoryFeatureVariation');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (46, 189);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (46, 190);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (46, 191);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (46, 192);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (193, @SO_accession, 'SO:0000409');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (194, @SO_term, 'binding_site');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (195, @ens_feature_class, 'Bio::EnsEMBL::Funcgen::MotifFeature');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (196, @ens_variant_class, 'Bio::EnsEMBL::Variation::MotifFeatureVariation');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (47, 193);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (47, 194);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (47, 195);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (47, 196);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (197, @ens_feature_class, 'Bio::EnsEMBL::ExternalFeature');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (198, @ens_feature_subtype, 'VISTA enhancer set');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (199, @ens_variant_class, 'Bio::EnsEMBL::Variation::ExternalFeatureVariation');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (48, 189);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (48, 190);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (48, 197);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (48, 198);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (48, 199);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (200, @ens_feature_subtype, 'cisRED motif');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (49, 189);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (49, 190);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (49, 197);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (49, 200);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (49, 199);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_61_62_l.sql|update to new transcript_variation schema');

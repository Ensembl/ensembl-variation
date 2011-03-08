INSERT INTO attrib_type (attrib_type_id, code, name, description) VALUES (1, 'SO_accession', '', 'Sequence Ontology accession');
INSERT INTO attrib_type (attrib_type_id, code, name, description) VALUES (2, 'SO_term', '', 'Sequence Ontology term');
INSERT INTO attrib_type (attrib_type_id, code, name, description) VALUES (3, 'display_term', '', 'Ensembl display term');
INSERT INTO attrib_type (attrib_type_id, code, name, description) VALUES (4, 'NCBI_term', '', 'NCBI term');
INSERT INTO attrib_type (attrib_type_id, code, name, description) VALUES (5, 'feature_SO_term', '', 'Sequence Ontology term for the associated feature');
INSERT INTO attrib_type (attrib_type_id, code, name, description) VALUES (6, 'rank', '', 'Relative severity of this variation consequence');

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (1, 1, 'SO:0001483');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (2, 2, 'SNV');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (3, 3, 'SNP');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (1, 1);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (1, 2);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (1, 3);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (4, 1, 'SO:1000002');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (5, 2, 'substitution');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (2, 4);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (2, 5);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (6, 1, 'SO:0001019');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (7, 2, 'copy_number_variation');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (8, 3, 'CNV');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (3, 6);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (3, 7);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (3, 8);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (9, 1, 'SO:0000667');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (10, 2, 'insertion');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (4, 9);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (4, 10);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (11, 1, 'SO:0000159');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (12, 2, 'deletion');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (5, 11);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (5, 12);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (13, 1, 'SO:1000032');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (14, 2, 'indel');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (6, 13);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (6, 14);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (15, 1, 'SO:0000705');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (16, 2, 'tandem_repeat');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (7, 15);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (7, 16);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (17, 1, 'SO:0001059');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (18, 2, 'sequence_alteration');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (8, 17);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (8, 18);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (19, 1, 'SO:0001628');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (20, 2, 'intergenic_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (21, 3, 'INTERGENIC');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (9, 19);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (9, 20);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (9, 21);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (22, 1, 'SO:0001635');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (23, 2, '5KB_upstream_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (24, 3, 'UPSTREAM');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (10, 22);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (10, 23);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (10, 24);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (25, 1, 'SO:0001633');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (26, 2, '5KB_downstream_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (27, 3, 'DOWNSTREAM');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (11, 25);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (11, 26);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (11, 27);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (28, 1, 'SO:0001636');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (29, 2, '2KB_upstream_variant');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (12, 28);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (12, 29);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (12, 24);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (30, 1, 'SO:0001634');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (31, 2, '500B_downstream_variant');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (13, 30);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (13, 31);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (13, 27);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (32, 1, 'SO:0001575');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (33, 2, 'splice_donor_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (34, 3, 'ESSENTIAL_SPLICE_SITE');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (14, 32);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (14, 33);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (14, 34);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (35, 1, 'splice_acceptor_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (36, 2, 'splice_acceptor_variant');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (15, 35);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (15, 36);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (15, 34);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (37, 1, 'SO:0001630');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (38, 2, 'splice_region_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (39, 3, 'SPLICE_SITE');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (16, 37);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (16, 38);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (16, 39);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (40, 1, 'SO:0001627');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (41, 2, 'intron_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (42, 3, 'INTRONIC');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (17, 40);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (17, 41);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (17, 42);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (43, 1, 'SO:0001623');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (44, 2, '5_prime_UTR_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (45, 3, '5PRIME_UTR');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (18, 43);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (18, 44);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (18, 45);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (46, 1, 'SO:0001624');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (47, 2, '3_prime_UTR_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (48, 3, '3PRIME_UTR');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (19, 46);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (19, 47);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (19, 48);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (49, 1, 'SO:0001577');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (50, 2, 'complex_change_in_transcript');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (51, 3, 'COMPLEX_INDEL');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (20, 49);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (20, 50);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (20, 51);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (52, 1, 'SO:0001588');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (53, 2, 'synonymous_codon');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (54, 3, 'SYNONYMOUS_CODING');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (21, 52);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (21, 53);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (21, 54);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (55, 1, 'SO:0001583');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (56, 2, 'non_synonymous_codon');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (57, 3, 'NON_SYNONYMOUS_CODING');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (22, 55);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (22, 56);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (22, 57);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (58, 1, 'SO:0001651');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (59, 2, 'inframe_codon_gain');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (23, 58);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (23, 59);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (23, 57);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (60, 1, 'SO:0001652');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (61, 2, 'inframe_codon_loss');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (24, 60);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (24, 61);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (24, 57);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (62, 1, 'SO:0001587');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (63, 2, 'stop_gained');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (64, 3, 'STOP_GAINED');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (25, 62);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (25, 63);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (25, 64);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (65, 1, 'SO:0001578');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (66, 2, 'stop_lost');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (67, 3, 'STOP_LOST');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (26, 65);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (26, 66);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (26, 67);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (68, 1, 'SO:0001567');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (69, 2, 'stop_retained_variant');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (27, 68);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (27, 69);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (27, 54);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (70, 1, 'SO:0001582');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (71, 2, 'initiator_codon_change');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (28, 70);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (28, 71);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (28, 57);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (72, 1, 'SO:0001589');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (73, 2, 'frameshift_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (74, 3, 'FRAMESHIFT_CODING');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (29, 72);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (29, 73);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (29, 74);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (75, 1, 'SO:0001626');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (76, 2, 'incomplete_terminal_codon_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (77, 3, 'PARTIAL_CODON');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (30, 75);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (30, 76);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (30, 77);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (78, 1, 'SO:0001621');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (79, 2, 'NMD_transcript_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (80, 3, 'NMD_TRANSCRIPT');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (31, 78);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (31, 79);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (31, 80);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (81, 1, 'SO:0001619');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (82, 2, 'nc_transcript_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (83, 3, 'WITHIN_NON_CODING_GENE');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (32, 81);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (32, 82);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (32, 83);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (84, 1, 'SO:0001620');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (85, 2, 'mature_miRNA_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (86, 3, 'WITHIN_MATURE_miRNA');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (33, 84);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (33, 85);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (33, 86);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (87, 1, 'SO:0001580');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (88, 2, 'coding_sequence_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (89, 3, 'CODING_UNKNOWN');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (34, 87);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (34, 88);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (34, 89);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (90, 1, 'SO:0001566');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (91, 2, 'regulatory_region_variant');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (92, 3, 'REGULATORY_REGION');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (35, 90);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (35, 91);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (35, 92);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (93, 1, 'SO:X000004');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (94, 2, 'miRNA_target_site_variant');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (36, 93);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (36, 94);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (36, 92);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (95, 1, 'SO:X000003');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (96, 2, 'binding_site_variant');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (37, 95);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (37, 96);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (37, 92);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (97, 1, 'SO:0000234');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (98, 2, 'mRNA');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (38, 97);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (38, 98);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (99, 1, 'SO:0000673');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (100, 2, 'transcript');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (39, 99);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (39, 100);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (101, 1, 'SO:0000185');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (102, 2, 'primary_transcript');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (40, 101);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (40, 102);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (103, 1, 'SO:0000655');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (104, 2, 'ncRNA');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (41, 103);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (41, 104);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (105, 1, 'SO:0000276');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (106, 2, 'miRNA');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (42, 105);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (42, 106);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (107, 1, 'SO:0005836');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (108, 2, 'regulatory_region');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (43, 107);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (43, 108);

INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (109, 1, 'SO:0000409');
INSERT INTO attrib (attrib_id, attrib_type_id, value) VALUES (110, 2, 'binding_site');
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (44, 109);
INSERT INTO attrib_set (attrib_set_id, attrib_id) VALUES (44, 110);


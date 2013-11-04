=head1 LICENSE

Copyright 2013 Ensembl

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::Variation::Utils::Config;

use base qw(Exporter); 

our @EXPORT_OK = qw(
    @ATTRIB_TYPES 
    %ATTRIBS
    @ATTRIB_SETS
    @VARIATION_CLASSES 
    @OVERLAP_CONSEQUENCES 
    @FEATURE_TYPES 
    $OVERLAP_CONSEQUENCE_CLASS
    $MAX_ATTRIB_CODE_LENGTH
);

our $OVERLAP_CONSEQUENCE_CLASS = 'Bio::EnsEMBL::Variation::OverlapConsequence';

our $MAX_ATTRIB_CODE_LENGTH = 20;

our @short_names = qw(1kg_hct 1kg_hct_ceu 1kg_hct_yri 1kg_hce 1kg_hce_ceu 1kg_hce_chb
                      1kg_hce_chd 1kg_hce_jpt 1kg_hce_lwk 1kg_hce_tsi 1kg_hce_yri 1kg_lc
	                      1kg_lc_ceu 1kg_lc_chb_jpt 1kg_lc_yri 1kg_hq hapmap
                      1kg 1kg_afr 1kg_amr 1kg_asn 1kg_eur 1kg_com 1kg_afr_com 1kg_amr_com 1kg_asn_com 1kg_eur_com
                      ind_venter ind_watson ind_gill ind_ak1 ind_irish ind_angrist
                      ind_gates_jr ind_gates_sr ind_kriek ind_quake ind_saqqaq ind_saqqaq_hc ind_sjk ind_yh
                      fail_all fail_nonref fail_ambig fail_gt_fq fail_incons_map fail_mult_map
                      fail_no_alleles fail_no_gt fail_no_map fail_no_seq fail_non_nt fail_mult_alleles fail_dbsnp_suspect
                      ph_hgmd_pub ph_nhgri ph_omim ph_variants ph_uniprot
                      ph_cosmic ph_ega  hapmap_ceu hapmap_hcb hapmap_jpt hapmap_yri
                      Affy_500K Affy_SNP6 Cardio-Metabo_Chip HumanOmni1-Quad Illumina_1M-duo Illumina_660Q
                      Human610_Quad HumanHap550 HumanHap650Y HumanOmni2.5 PorcineSNP60
                      esp_6500 clin_assoc 
                      Chicken600K EquineSNP50 BovineHD BovineLD BovineSNP50  
                      phencode HumanOmni5 OvineSNP50 OvineHDSNP
                     );

our @dbsnp_clinical_significance_types = qw(
    unknown
    untested
    non-pathogenic
    probable-non-pathogenic
    probable-pathogenic
    pathogenic
    drug-response
    histocompatibility
    other
);

our @dgva_clinical_significance_types = (
    'Not tested',
    'Benign',
    'Pathogenic',
    'Uncertain Significance',
    'Uncertain Significance: likely benign',
    'Uncertain Significance: likely pathogenic'
);

our @VARIATION_CLASSES = (
    {
        SO_accession => 'SO:0001483',
        SO_term => 'SNV',
        display_term => 'SNP',
        somatic_display_term => 'somatic_SNV',
    },
    {
        SO_accession => 'SO:1000002',
        SO_term => 'substitution',
    },
    {
        SO_accession => 'SO:0001019',
        SO_term => 'copy_number_variation',
        display_term => 'CNV',
    },
    {
        SO_accession => 'SO:0000667',
        SO_term => 'insertion',
    },
    {
        SO_accession => 'SO:0000159',
        SO_term => 'deletion',
    },
    {
        SO_accession => 'SO:1000032',
        SO_term => 'indel',
    },
    {
        SO_accession => 'SO:0000705',
        SO_term => 'tandem_repeat',
    },
    {
        SO_accession => 'SO:0001059',
        SO_term => 'sequence_alteration',
    },
    # Structural variation classes
    {
        SO_accession => 'SO:0001537',
        SO_term => 'structural_variant',
        display_term => 'SV',
    },
    {
        SO_accession => 'SO:0000051',
        SO_term => 'probe',
        display_term => 'CNV_PROBE',
    },
    {
        SO_accession => 'SO:0001742',
        SO_term => 'copy_number_gain',
        display_term => 'Gain',
    },
    {
        SO_accession => 'SO:0001743',
        SO_term => 'copy_number_loss',
        display_term => 'Loss',
    },
    {
        SO_accession => 'SO:1000036',
        SO_term => 'inversion',
    },
    {
        SO_accession => 'SO:0001784',
        SO_term => 'complex_structural_alteration',
        display_term => 'Complex',
    },
    {
        SO_accession => 'SO:1000173',
        SO_term => 'tandem_duplication',
        display_term => 'Tandem duplication',
    },
    {
        SO_accession => 'SO:0001837',
        SO_term => 'mobile_element_insertion',
        display_term => 'Mobile element insertion',
    },
    {
        SO_accession => 'SO:0001873',
        SO_term => 'interchromosomal_breakpoint',
        display_term => 'Interchromosomal breakpoint',
    },   
    {
        SO_accession => 'SO:0001874',
        SO_term => 'intrachromosomal_breakpoint',
        display_term => 'Intrachromosomal breakpoint',
    },
    {
        SO_accession => 'SO:0000199',
        SO_term => 'translocation',
    },
    {
        SO_accession => 'SO:1000035',
        SO_term => 'duplication',
        display_term => 'Duplication',
    },
);

our @OVERLAP_CONSEQUENCES = (
    {
        SO_accession => 'SO:0001628',
        SO_term => 'intergenic_variant',
        display_term => 'INTERGENIC',
        rank => '38',
        tier => '4',
        description => 'A sequence variant located in the intergenic region, between genes',
        label => 'Intergenic variant',
        is_default => 1,
    },
    {
        SO_accession => 'SO:0001631',
        SO_term => 'upstream_gene_variant',
        display_term => 'UPSTREAM',
        feature_SO_term => 'transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '24',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream',
        description => 'A sequence variant located 5\' of a gene',
        label => 'Upstream gene variant',
    },
    {
        SO_accession => 'SO:0001632',
        SO_term => 'downstream_gene_variant',
        display_term => 'DOWNSTREAM',
        feature_SO_term => 'transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '25',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream',
        description => 'A sequence variant located 3\' of a gene',
        label => 'Downstream gene variant',
    },
    {
        SO_accession => 'SO:0001575',
        SO_term => 'splice_donor_variant',
        display_term => 'ESSENTIAL_SPLICE_SITE',
        NCBI_term => 'splice-5',
        feature_SO_term => 'primary_transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '3',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::donor_splice_site',
        description => 'A splice variant that changes the 2 base region at the 5\' end of an intron',
        label => 'Splice donor variant',
    },
    {
        SO_accession => 'SO:0001574',
        SO_term => 'splice_acceptor_variant',
        display_term => 'ESSENTIAL_SPLICE_SITE',
        NCBI_term => 'splice-3',
        feature_SO_term => 'primary_transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '3',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::acceptor_splice_site',
        description => 'A splice variant that changes the 2 base region at the 3\' end of an intron',
        label => 'Splice acceptor variant',
    },
    {
        SO_accession => 'SO:0001630',
        SO_term => 'splice_region_variant',
        display_term => 'SPLICE_SITE',
        feature_SO_term => 'primary_transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '13',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::splice_region',
        description => 'A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron',
        label => 'Splice region variant',
    },
    {
        SO_accession => 'SO:0001627',
        SO_term => 'intron_variant',
        display_term => 'INTRONIC',
        NCBI_term => 'intron',
        feature_SO_term => 'primary_transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '21',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_intron',
        description => 'A transcript variant occurring within an intron',
        label => 'Intron variant',
    },
    {
        SO_accession => 'SO:0001623',
        SO_term => '5_prime_UTR_variant',
        display_term => '5PRIME_UTR',
        NCBI_term => 'untranslated_5',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '18',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_5_prime_utr',
        description => 'A UTR variant of the 5\' UTR',
        label => '5 prime UTR variant',
    },
    {
        SO_accession => 'SO:0001624',
        SO_term => '3_prime_UTR_variant',
        display_term => '3PRIME_UTR',
        NCBI_term => 'untranslated_3',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '19',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_3_prime_utr',
        description => 'A UTR variant of the 3\' UTR',
        label => '3 prime UTR variant',
    },
#    {
#        SO_accession => 'SO:0001577',
#        SO_term => 'complex_change_in_transcript',
#        display_term => 'COMPLEX_INDEL',
#        feature_SO_term => 'primary_transcript',
#        feature_class => 'Bio::EnsEMBL::Transcript',
#        variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
#        rank => '4',
#        tier => '3',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::complex_indel',
#        description => 'Insertion or deletion that spans an exon/intron or coding sequence/UTR border',
#        label => 'Complex change in transcript',
#    },
    {
        SO_accession => 'SO:0001819',
        SO_term => 'synonymous_variant',
        display_term => 'SYNONYMOUS_CODING',
        NCBI_term => 'cds-synon',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '15',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::synonymous_variant',
        description => 'A sequence variant where there is no resulting change to the encoded amino acid',
        label => 'Synonymous variant',
    },
    {
        SO_accession => 'SO:0001583',
        SO_term => 'missense_variant',
        display_term => 'NON_SYNONYMOUS_CODING',
        NCBI_term => 'missense',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '12',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::missense_variant',
        description => 'A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved',
        label => 'Missense variant',
    },
    {
        SO_accession => 'SO:0001821',
        SO_term => 'inframe_insertion',
        display_term => 'NON_SYNONYMOUS_CODING',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '10',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_insertion',
        description => 'An inframe non synonymous variant that inserts bases into in the coding sequence',
        label => 'Inframe insertion',
    },
    {
        SO_accession => 'SO:0001822',
        SO_term => 'inframe_deletion',
        display_term => 'NON_SYNONYMOUS_CODING',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '11',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_deletion',
        description => 'An inframe non synonymous variant that deletes bases from the coding sequence',
        label => 'Inframe deletion',
    },
    {
        SO_accession => 'SO:0001587',
        SO_term => 'stop_gained',
        display_term => 'STOP_GAINED',
        NCBI_term => 'nonsense',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '4',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_gained',
        description => 'A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript',
        label => 'Stop gained',
    },
    {
        SO_accession => 'SO:0001578',
        SO_term => 'stop_lost',
        display_term => 'STOP_LOST',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '6',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_lost',
        description => 'A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript',
        label => 'Stop lost',
    },
    {
        SO_accession => 'SO:0001567',
        SO_term => 'stop_retained_variant',
        display_term => 'SYNONYMOUS_CODING',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '15',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_retained',
        description => 'A sequence variant where at least one base in the terminator codon is changed, but the terminator remains',
        label => 'Stop retained variant',
    },
    {
        SO_accession => 'SO:0001582',
        SO_term => 'initiator_codon_variant',
        display_term => 'NON_SYNONYMOUS_CODING',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '7',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::affects_start_codon',
        description => 'A codon variant that changes at least one base of the first codon of a transcript',
        label => 'Initiator codon variant',
    },
    {
        SO_accession => 'SO:0001589',
        SO_term => 'frameshift_variant',
        display_term => 'FRAMESHIFT_CODING',
        NCBI_term => 'frameshift',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '5',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::frameshift',
        description => 'A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three',
        label => 'Frameshift variant',
    },
    {
        SO_accession => 'SO:0001626',
        SO_term => 'incomplete_terminal_codon_variant',
        display_term => 'PARTIAL_CODON',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '14',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::partial_codon',
        description => 'A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed',
        label => 'Incomplete terminal codon variant',
    },
    {
        SO_accession => 'SO:0001621',
        SO_term => 'NMD_transcript_variant',
        display_term => 'NMD_TRANSCRIPT',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '22',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_nmd_transcript',
        description => 'A variant in a transcript that is the target of NMD',
        label => 'NMD transcript variant',
    },
    {
        SO_accession => 'SO:0001619',
        SO_term => 'nc_transcript_variant',
        display_term => 'WITHIN_NON_CODING_GENE',
        feature_SO_term => 'ncRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '23',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_non_coding_gene',
        description => 'A transcript variant of a non coding RNA',
        label => 'NC transcript variant',
    },
    {
        SO_accession => 'SO:0001792',
        SO_term => 'non_coding_exon_variant',
        display_term => 'WITHIN_NON_CODING_GENE',
        feature_SO_term => 'ncRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '20',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::non_coding_exon_variant',
        description => 'A sequence variant that changes non-coding exon sequence',
        label => 'Non coding exon variant',
    },
    {
        SO_accession => 'SO:0001620',
        SO_term => 'mature_miRNA_variant',
        display_term => 'WITHIN_MATURE_miRNA',
        feature_SO_term => 'miRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '17',
        tier => '2',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_mature_miRNA',
        description => 'A transcript variant located with the sequence of the mature miRNA',
        label => 'Mature miRNA variant',
    },
    {
        SO_accession => 'SO:0001580',
        SO_term => 'coding_sequence_variant',
        display_term => 'CODING_UNKNOWN',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '16',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::coding_unknown',
        description => 'A sequence variant that changes the coding sequence',
        label => 'Coding sequence variant',
    },
    {
        SO_accession => 'SO:0001566',
        SO_term => 'regulatory_region_variant',
        display_term => 'REGULATORY_REGION',
        feature_SO_term => 'regulatory_region',
        feature_class => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '36',
        tier => '2',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_regulatory_feature',
        description => 'A sequence variant located within a regulatory region',
        label => 'Regulatory region variant',
    },
#    {
#        SO_accession => 'SO:X000005',
#        SO_term => 'pre_miRNA_variant',
#        display_term => 'WITHIN_NON_CODING_GENE',
#        feature_SO_term => 'miRNA',
#        feature_class => 'Bio::EnsEMBL::Transcript',
#        rank => '13',
#        tier => '2',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_miRNA',
#    },
#    {
#        SO_accession => 'SO:X000004',
#        SO_term => 'miRNA_target_site_variant',
#        display_term => 'REGULATORY_REGION',
#        feature_SO_term => 'binding_site',
#        feature_class => 'Bio::EnsEMBL::Funcgen::ExternalFeature',
#        rank => '13',
#        tier => '2',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_miRNA_target_site',
#        description => 'In regulatory region annotated by Ensembl',
#        label => 'Regulatory region',
#    },
    {
        SO_accession => 'SO:0001782',
        SO_term => 'TF_binding_site_variant',
        display_term => 'REGULATORY_REGION',
        feature_SO_term => 'TF_binding_site',
        feature_class => 'Bio::EnsEMBL::Funcgen::MotifFeature',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '30',
        tier => '2',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_motif_feature',
        description => 'A sequence variant located within a transcription factor binding site',
        label => 'TF binding site',
    },

#    {
#        SO_accession => 'SO:X000002',
#        SO_term => 'decreased_binding_affinity',
#        display_term => 'REGULATORY_REGION',
#        feature_SO_term => 'binding_site',
#        feature_class => 'Bio::EnsEMBL::Funcgen::MotifFeature',
#        rank => '47',
#        tier => '2',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::decreased_binding_affinity',
#    },
#    {
#        SO_accession => 'SO:X000001',
#        SO_term => 'increased_binding_affinity',
#        display_term => 'REGULATORY_REGION',
#        feature_SO_term => 'binding_site',
#        feature_class => 'Bio::EnsEMBL::Funcgen::MotifFeature',
#        rank => '48',
#        tier => '2',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::increased_binding_affinity',
#    },


    ## NEW FOR 68
    #############
    
    {
        SO_accession => 'SO:0001893',
        SO_term => 'transcript_ablation',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '1',
        tier => '1',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
        description => 'A feature ablation whereby the deleted region includes a transcript feature',
        label => 'Transcript ablation',
    },
#    {
#        SO_accession => 'SO:0001886',
#        SO_term => 'transcript_fusion',
#        feature_SO_term => 'mRNA',
#        feature_class => 'Bio::EnsEMBL::Transcript',
#        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
#        rank => '2',
#        tier => '2',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::transcript_fusion',
#        description => 'A feature fusion where the deletion brings together transcript regions',
#        label => 'Transcript fusion',
#    },
    {
        SO_accession => 'SO:0001889',
        SO_term => 'transcript_amplification',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '8',
        tier => '1',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
        description => 'A feature amplification of a region containing a transcript',
        label => 'Transcript amplification',
    },
#    {
#        SO_accession => 'SO:0001883',
#        SO_term => 'transcript_translocation',
#        feature_SO_term => 'mRNA',
#        feature_class => 'Bio::EnsEMBL::Transcript',
#        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
#        rank => '9',
        #tier => '2',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::transcript_translocation',
#        description => 'A feature translocation where the region contains a transcript',
#        label => 'Transcript translocation',
#    },
    {
        SO_accession => 'SO:0001895',
        SO_term => 'TFBS_ablation',
        feature_SO_term => 'TF_binding_site',
        feature_class => 'Bio::EnsEMBL::Funcgen::MotifFeature',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '26',
        tier => '2',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
        description => 'A feature ablation whereby the deleted region includes a transcription factor binding site',
        label => 'TFBS ablation',
    },
#    {
#        SO_accession => 'SO:0001888',
#        SO_term => 'TFBS_fusion',
#        feature_SO_term => 'TF_binding_site',
#        feature_class => 'Bio::EnsEMBL::Funcgen::MotifFeature',
#        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
#        rank => '27',
        #tier => '2',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::motif_feature_fusion',
#        description => 'A fusion where the deletion brings together transcription factor binding sites',
#        label => 'TFBS fusion',
#    },
    {
        SO_accession => 'SO:0001892',
        SO_term => 'TFBS_amplification',
        feature_SO_term => 'TF_binding_site',
        feature_class => 'Bio::EnsEMBL::Funcgen::MotifFeature',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '28',
        tier => '2',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
        description => 'A feature amplification of a region containing a transcription factor binding site',
        label => 'TFBS amplification',
    },
#    {
#        SO_accession => 'SO:0001885',
#        SO_term => 'TFBS_translocation',
#        feature_SO_term => 'TF_binding_site',
#        feature_class => 'Bio::EnsEMBL::Funcgen::MotifFeature',
#        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
#        rank => '29',
        #tier => '2',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::motif_feature_translocation',
#        description => 'A feature translocation where the region contains a transcription factor binding site',
#        label => 'TFBS translocation',
#    },
    {
        SO_accession => 'SO:0001894',
        SO_term => 'regulatory_region_ablation',
        feature_SO_term => 'TF_binding_site',
        feature_class => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '31',
        tier => '2',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
        description => 'A feature ablation whereby the deleted region includes a regulatory region',
        label => 'Regulatory region ablation',
    },
#    {
#        SO_accession => 'SO:0001887',
#        SO_term => 'regulatory_region_fusion',
#        feature_SO_term => 'TF_binding_site',
#        feature_class => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
#        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
#        rank => '32',
        #tier => '2',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::regulatory_feature_fusion',
#        description => 'A fusion where the deletion brings together regulatory regions',
#        label => 'Regulatory region fusion',
#    },
    {
        SO_accession => 'SO:0001891',
        SO_term => 'regulatory_region_amplification',
        feature_SO_term => 'TF_binding_site',
        feature_class => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '33',
        tier => '2',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
        description => 'A feature amplification of a region containing a regulatory region',
        label => 'Regulatory region amplification',
    },
#    {
#        SO_accession => 'SO:0001884',
#        SO_term => 'regulatory_region_translocation',
#        feature_SO_term => 'TF_binding_site',
#        feature_class => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
#        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
#        rank => '34',
        #tier => '2',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::regulatory_feature_translocation',
#        description => 'A feature translocation where the region contains a regulatory region',
#        label => 'Regulatory region translocation',
#    },
    {
        SO_accession => 'SO:0001907',
        SO_term => 'feature_elongation',
        feature_SO_term => 'sequence_feature',
        feature_class => 'Bio::EnsEMBL::Feature',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '36',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_elongation',
        description => 'A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence',
        label => 'Feature elongation',
    },
    {
        SO_accession => 'SO:0001906',
        SO_term => 'feature_truncation',
        feature_SO_term => 'sequence_feature',
        feature_class => 'Bio::EnsEMBL::Feature',
        variant_feature_class => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
        rank => '37',
        tier => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_truncation',
        description => 'A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence',
        label => 'Feature truncation',
    },
);

our @FEATURE_TYPES = (
    {
        SO_accession => 'SO:0000234',
        SO_term => 'mRNA',
        ens_feature_class => 'Bio::EnsEMBL::Transcript',
        ens_feature_subtype => 'protein_coding',
        ens_variant_class => 'Bio::EnsEMBL::Variation::TranscriptVariation',
    },
    {
        SO_accession => 'SO:0000673',
        SO_term => 'transcript',
        ens_feature_class => 'Bio::EnsEMBL::Transcript',
        ens_variant_class => 'Bio::EnsEMBL::Variation::TranscriptVariation',
    },
    {
        SO_accession => 'SO:0000185',
        SO_term => 'primary_transcript',
        ens_feature_class => 'Bio::EnsEMBL::Transcript',
        ens_variant_class => 'Bio::EnsEMBL::Variation::TranscriptVariation',
    },
    {
        SO_accession => 'SO:0000655',
        SO_term => 'ncRNA',
        ens_feature_class => 'Bio::EnsEMBL::Transcript',
        ens_variant_class => 'Bio::EnsEMBL::Variation::TranscriptVariation',
    },
    {
        SO_accession => 'SO:0000276',
        SO_term => 'miRNA',
        ens_feature_class => 'Bio::EnsEMBL::Transcript',
        ens_variant_class => 'Bio::EnsEMBL::Variation::TranscriptVariation',
    },
    {
        SO_accession => 'SO:0005836',
        SO_term => 'regulatory_region',
        ens_feature_class => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
        ens_variant_class => 'Bio::EnsEMBL::Variation::RegulatoryFeatureVariation',
    },
    {
        SO_accession => 'SO:0000409',
        SO_term => 'binding_site',
        ens_feature_class => 'Bio::EnsEMBL::Funcgen::MotifFeature',
        ens_variant_class => 'Bio::EnsEMBL::Variation::MotifFeatureVariation',
    },
    {
        SO_accession => 'SO:0005836',
        SO_term => 'regulatory_region',
        ens_feature_class => 'Bio::EnsEMBL::Funcgen::ExternalFeature',
        ens_variant_class => 'Bio::EnsEMBL::Variation::ExternalFeatureVariation',
        ens_feature_subtype => 'VISTA enhancer set',
    },
    {
        SO_accession => 'SO:0000409',
        SO_term => 'binding_site',
        ens_feature_class => 'Bio::EnsEMBL::Funcgen::ExternalFeature',
        ens_variant_class => 'Bio::EnsEMBL::Variation::ExternalFeatureVariation',
        ens_feature_subtype => 'cisRED motif',
    },
    {
        SO_accession => 'SO:0005836',
        SO_term => 'regulatory_region',
        ens_feature_class => 'Bio::EnsEMBL::Funcgen::ExternalFeature',
        ens_variant_class => 'Bio::EnsEMBL::Variation::ExternalFeatureVariation',
        ens_feature_subtype => 'miRanda miRNA target',
    },
    {
        SO_accession => 'SO:0000110',
        SO_term => 'sequence_feature',
        ens_feature_class => 'Bio::EnsEMBL::Feature',
        ens_variant_class => 'Bio::EnsEMBL::Variation::StructuralVariationFeatureOverlap',
    },
);

# attrib_types are specified as hashrefs in the @ATTRIB_TYPES array. Each hashref should have a value for the key 'code' and optionally values for the keys 'name' and 'description'
our @ATTRIB_TYPES = (
    {
        code => 'SO_accession',
        description => 'Sequence Ontology accession',
    }, 
    {
        code => 'SO_term',
        description => 'Sequence Ontology term',
    },
    {
        code => 'display_term',
        description => 'Ensembl display term',
    },
    {
        code => 'NCBI_term',
        description => 'NCBI term',
    },
    {
        code => 'feature_SO_term',
        description => 'Sequence Ontology term for the associated feature',
    },
    {
        code => 'rank',
        description => 'Relative severity of this variation consequence',
    },
    {
        code => 'polyphen_prediction',
        description => 'PolyPhen-2 prediction',
    },
    {
        code => 'sift_prediction',
        description => 'SIFT prediction',
    },
    {
        code => 'short_name',
        name => 'Short name',
        description => 'A shorter name for an instance, e.g. a VariationSet',
    },
    {
        code => 'dbsnp_clin_sig',
        name => 'dbSNP/ClinVar clinical significance',
        description => 'The clinical significance of a variant as reported by ClinVar and dbSNP',
    },
    {
        code => 'dgva_clin_sig',
        name => 'DGVa clinical significance',
        description => 'The clinical significance of a structural variant as reported by DGVa',
    },
    {
        code => 'prot_func_analysis',
        name => 'Protein function analysis ',
        description => 'The program used to make protein function predictions',
    },
    {
        code => 'associated_gene',
        name => 'Associated gene',
        description => 'ID of gene(s) linked by phenotype association',
    },
    {
        code => 'risk_allele',
        name => 'Risk allele',
        description => 'Risk allele in phenotype association',
    },
    {
        code => 'p_value',
        name => 'P-value',
        description => 'P-value denoting significance of an observed phenotype annotation',
    },
    {
        code => 'variation_names',
        name => 'Variation names',
        description => 'Variant ID(s) linked with a phenotype association',
    },
	{
        code => 'sample_id',
        name => 'Sample ID',
        description => 'Sample ID for source of phenotype association',
    },
	{
        code => 'strain_id',
        name => 'Strain ID',
        description => 'Strain ID for source of phenotype association',
    },
    {
        code => 'lod_score',
        name => 'LOD score',
        description => 'Log Of Odds score',
    },
    {
        code => 'variance',
        name => 'Variance',
        description => 'Variance statistic',
    },
    {
        code => 'inheritance_type',
        name => 'Inheritance type',
        description => 'Inheritance type of a trait',
    },
    {
        code => 'external_id',
        name => 'External ID',
        description => 'External identifier for an entity',
    },
    {
        code => 'odds_ratio',
        name => 'Odds ratio',
        description => 'Odds ratio used to denote significance of an observed phenotype annotation',
    },
    {
        code => 'beta_coef',
        name => 'Beta coefficient',
        description => 'Beta coefficient (or standardized coefficient) used to denote significance of an observed phenotype annotation',
    },
	{
		code => 'allele_symbol',
		name => 'Allele symbol',
		description => 'Allele symbol linked with phenotype association',
	},
	{
		code => 'allele_accession_id',
		name => 'Allele accession ID',
		description => 'Allele accession ID linked with phenotype association',
	},
	{
		code => 'marker_accession_id',
		name => 'Marker accession ID',
		description => 'Marker ID linked with phenotype association',
	},

);

# attribs are specified in the %ATTRIBS hash, having the attrib_type code as hash key and a listref containing the attribs that will be loaded as value
our %ATTRIBS = (
   'short_name'          => \@short_names,
   'dbsnp_clin_sig'      => \@dbsnp_clinical_significance_types,
   'dgva_clin_sig'       => \@dgva_clinical_significance_types,
   'polyphen_prediction' => ['probably damaging', 'possibly damaging', 'benign', 'unknown'],
   'sift_prediction'     => [qw(tolerated deleterious)],
   'prot_func_analysis'  => [qw(sift polyphen_humvar polyphen_humdiv)],
);

# attrib sets are specified by putting a hashref in the @ATTRIB_SETS array having the attrib_type code as key and the attrib as value. new attrib entries will be inserted as necessary
our @ATTRIB_SETS = (
    @VARIATION_CLASSES,
    @OVERLAP_CONSEQUENCES,
    @FEATURE_TYPES
);

1;

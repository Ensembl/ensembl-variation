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
                      1kg_lc_ceu 1kg_lc_chb_jpt 1kg_lc_yri hapmap
                      ind_venter ind_watson ind_gill ind_ak1 ind_irish ind_angrist
                      ind_gates_jr ind_gates_sr ind_kriek ind_quake ind_saqqaq ind_saqqaq_hc ind_sjk ind_yh
                      fail_all fail_nonref fail_ambig fail_gt_fq fail_incons_map fail_mult_map
                      fail_no_alleles fail_no_gt fail_no_map fail_no_seq fail_non_nt fail_mult_alleles fail_dbsnp_suspect
                      ph_hgmd_pub ph_johnson_et_al ph_nhgri ph_omim ph_variants ph_uniprot
                      ph_cosmic ph_ega precious hapmap_ceu hapmap_hcb hapmap_jpt hapmap_yri
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
    
);

our @OVERLAP_CONSEQUENCES = (
    {
        SO_accession => 'SO:0001628',
        SO_term => 'intergenic_variant',
        display_term => 'INTERGENIC',
        rank => '26',
        description => 'More than 5 kb either upstream or downstream of a transcript',
        label => 'Intergenic',
        is_default => 1,
    },
    {
        SO_accession => 'SO:0001635',
        SO_term => '5KB_upstream_variant',
        display_term => 'UPSTREAM',
        feature_SO_term => 'transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '22',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream_5KB',
        description => 'Within 5 kb upstream of the 5 prime end of a transcript',
        label => 'Upstream',
    },
    {
        SO_accession => 'SO:0001633',
        SO_term => '5KB_downstream_variant',
        display_term => 'DOWNSTREAM',
        feature_SO_term => 'transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '23',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream_5KB',
        description => 'Within 5 kb downstream of the 3 prime end of a transcript',
        label => 'Downstream',
    },
    {
        SO_accession => 'SO:0001636',
        SO_term => '2KB_upstream_variant',
        display_term => 'UPSTREAM',
        NCBI_term => 'near-gene-5',
        feature_SO_term => 'transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '20',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream_2KB',
        description => 'Within 5 kb upstream of the 5 prime end of a transcript',
        label => 'Upstream',
    },
    {
        SO_accession => 'SO:0001634',
        SO_term => '500B_downstream_variant',
        display_term => 'DOWNSTREAM',
        NCBI_term => 'near-gene-3',
        feature_SO_term => 'transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '21',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream_500B',
        description => 'Within 5 kb downstream of the 3 prime end of a transcript',
        label => 'Downstream',
    },
    {
        SO_accession => 'SO:0001575',
        SO_term => 'splice_donor_variant',
        display_term => 'ESSENTIAL_SPLICE_SITE',
        NCBI_term => 'splice-5',
        feature_SO_term => 'primary_transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '1',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::donor_splice_site',
        description => 'In the first 2 or the last 2 basepairs of an intron',
        label => 'Essential splice site',
    },
    {
        SO_accession => 'SO:0001574',
        SO_term => 'splice_acceptor_variant',
        display_term => 'ESSENTIAL_SPLICE_SITE',
        NCBI_term => 'splice-3',
        feature_SO_term => 'primary_transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '1',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::acceptor_splice_site',
        description => 'In the first 2 or the last 2 basepairs of an intron',
        label => 'Essential splice site',
    },
    {
        SO_accession => 'SO:0001630',
        SO_term => 'splice_region_variant',
        display_term => 'SPLICE_SITE',
        feature_SO_term => 'primary_transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '10',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::splice_region',
        description => '1-3 bps into an exon or 3-8 bps into an intron',
        label => 'Splice site',
    },
    {
        SO_accession => 'SO:0001627',
        SO_term => 'intron_variant',
        display_term => 'INTRONIC',
        NCBI_term => 'intron',
        feature_SO_term => 'primary_transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '17',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_intron',
        description => 'In intron',
        label => 'Intronic',
    },
    {
        SO_accession => 'SO:0001623',
        SO_term => '5_prime_UTR_variant',
        display_term => '5PRIME_UTR',
        NCBI_term => 'untranslated_5',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '15',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_5_prime_utr',
        description => 'In 5 prime untranslated region',
        label => '5 prime UTR',
    },
    {
        SO_accession => 'SO:0001624',
        SO_term => '3_prime_UTR_variant',
        display_term => '3PRIME_UTR',
        NCBI_term => 'untranslated_3',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '16',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_3_prime_utr',
        description => 'In 3 prime untranslated region',
        label => '3 prime UTR',
    },
    {
        SO_accession => 'SO:0001577',
        SO_term => 'complex_change_in_transcript',
        display_term => 'COMPLEX_INDEL',
        feature_SO_term => 'primary_transcript',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '4',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::complex_indel',
        description => 'Insertion or deletion that spans an exon/intron or coding sequence/UTR border',
        label => 'Complex in/del',
    },
    {
        SO_accession => 'SO:0001588',
        SO_term => 'synonymous_codon',
        display_term => 'SYNONYMOUS_CODING',
        NCBI_term => 'cds-synon',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '12',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::synonymous_codon',
        description => 'In coding sequence, not resulting in an amino acid change (silent mutation)',
        label => 'Synonymous coding',
    },
    {
        SO_accession => 'SO:0001583',
        SO_term => 'non_synonymous_codon',
        display_term => 'NON_SYNONYMOUS_CODING',
        NCBI_term => 'missense',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '9',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::non_synonymous_codon',
        description => 'In coding sequence and results in an amino acid change in the encoded peptide sequence',
        label => 'Non-synonymous coding',
    },
    {
        SO_accession => 'SO:0001651',
        SO_term => 'inframe_codon_gain',
        display_term => 'NON_SYNONYMOUS_CODING',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '8',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_codon_gain',
        description => 'In coding sequence and results in an amino acid change in the encoded peptide sequence',
        label => 'Non-synonymous coding',
    },
    {
        SO_accession => 'SO:0001652',
        SO_term => 'inframe_codon_loss',
        display_term => 'NON_SYNONYMOUS_CODING',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '7',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_codon_loss',
        description => 'In coding sequence and results in an amino acid change in the encoded peptide sequence',
        label => 'Non-synonymous coding',
    },
    {
        SO_accession => 'SO:0001587',
        SO_term => 'stop_gained',
        display_term => 'STOP_GAINED',
        NCBI_term => 'nonsense',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '2',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_gained',
        description => 'In coding sequence, resulting in the gain of a stop codon',
        label => 'Stop gained',
    },
    {
        SO_accession => 'SO:0001578',
        SO_term => 'stop_lost',
        display_term => 'STOP_LOST',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '3',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_lost',
        description => 'In coding sequence, resulting in the loss of a stop codon',
        label => 'Stop lost',
    },
    {
        SO_accession => 'SO:0001567',
        SO_term => 'stop_retained_variant',
        display_term => 'SYNONYMOUS_CODING',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '12',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_retained',
        description => 'In coding sequence, not resulting in an amino acid change (silent mutation)',
        label => 'Synonymous coding',
    },
    {
        SO_accession => 'SO:0001582',
        SO_term => 'initiator_codon_change',
        display_term => 'NON_SYNONYMOUS_CODING',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '6',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::affects_start_codon',
        description => 'In coding sequence and results in an amino acid change in the encoded peptide sequence',
        label => 'Non-synonymous coding',
    },
    {
        SO_accession => 'SO:0001589',
        SO_term => 'frameshift_variant',
        display_term => 'FRAMESHIFT_CODING',
        NCBI_term => 'frameshift',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '5',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::frameshift',
        description => 'In coding sequence, resulting in a frameshift',
        label => 'Frameshift coding',
    },
    {
        SO_accession => 'SO:0001626',
        SO_term => 'incomplete_terminal_codon_variant',
        display_term => 'PARTIAL_CODON',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '11',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::partial_codon',
        description => 'Located within the final, incomplete codon of a transcript whose end coordinate is unknown',
        label => 'Partial codon',
    },
    {
        SO_accession => 'SO:0001621',
        SO_term => 'NMD_transcript_variant',
        display_term => 'NMD_TRANSCRIPT',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '18',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_nmd_transcript',
        description => 'Located within a transcript predicted to undergo nonsense-mediated decay',
        label => 'NMD transcript',
    },
    {
        SO_accession => 'SO:0001619',
        SO_term => 'nc_transcript_variant',
        display_term => 'WITHIN_NON_CODING_GENE',
        feature_SO_term => 'ncRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '19',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_non_coding_gene',
        description => 'Located within a gene that does not code for a protein',
        label => 'Within non-coding gene',
    },
    {
        SO_accession => 'SO:0001620',
        SO_term => 'mature_miRNA_variant',
        display_term => 'WITHIN_MATURE_miRNA',
        feature_SO_term => 'miRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '14',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_mature_miRNA',
        description => 'Located within a microRNA',
        label => 'Within mature miRNA',
    },
    {
        SO_accession => 'SO:0001580',
        SO_term => 'coding_sequence_variant',
        display_term => 'CODING_UNKNOWN',
        feature_SO_term => 'mRNA',
        feature_class => 'Bio::EnsEMBL::Transcript',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '13',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::coding_unknown',
        description => 'In coding sequence with indeterminate effect',
        label => 'Coding unknown',
    },
    {
        SO_accession => 'SO:0001566',
        SO_term => 'regulatory_region_variant',
        display_term => 'REGULATORY_REGION',
        feature_SO_term => 'regulatory_region',
        feature_class => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '25',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_regulatory_feature',
        description => 'In regulatory region annotated by Ensembl',
        label => 'Regulatory region',
    },
#    {
#        SO_accession => 'SO:X000005',
#        SO_term => 'pre_miRNA_variant',
#        display_term => 'WITHIN_NON_CODING_GENE',
#        feature_SO_term => 'miRNA',
#        feature_class => 'Bio::EnsEMBL::Transcript',
#        rank => '13',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_miRNA',
#    },
#    {
#        SO_accession => 'SO:X000004',
#        SO_term => 'miRNA_target_site_variant',
#        display_term => 'REGULATORY_REGION',
#        feature_SO_term => 'binding_site',
#        feature_class => 'Bio::EnsEMBL::Funcgen::ExternalFeature',
#        rank => '13',
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
		variant_feature_class => 'Bio::EnsEMBL::Variation::VariationFeature',
        rank => '24',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_motif_feature',
        description => 'In regulatory region annotated by Ensembl',
        label => 'Regulatory region',
    },
#    {
#        SO_accession => 'SO:0001566',
#        SO_term => 'regulatory_region_variant',
#        display_term => 'REGULATORY_REGION',
#        feature_SO_term => 'regulatory_region',
#        feature_class => 'Bio::EnsEMBL::Funcgen::ExternalFeature',
#        rank => '50',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_external_feature',
#        description => 'In regulatory region annotated by Ensembl',
#        label => 'Regulatory region',
#    },

#    {
#        SO_accession => 'SO:X000002',
#        SO_term => 'decreased_binding_affinity',
#        display_term => 'REGULATORY_REGION',
#        feature_SO_term => 'binding_site',
#        feature_class => 'Bio::EnsEMBL::Funcgen::MotifFeature',
#        rank => '47',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::decreased_binding_affinity',
#    },
#    {
#        SO_accession => 'SO:X000001',
#        SO_term => 'increased_binding_affinity',
#        display_term => 'REGULATORY_REGION',
#        feature_SO_term => 'binding_site',
#        feature_class => 'Bio::EnsEMBL::Funcgen::MotifFeature',
#        rank => '48',
#        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::increased_binding_affinity',
#    },


	## STRUCTURAL VARIATION
	#######################
    {
        SO_accession => 'SO:X000101',
        SO_term => 'partial_overlap',
        feature_SO_term => 'sequence_feature',
        feature_class => 'Bio::EnsEMBL::Feature',
		variant_feature_class => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
        rank => '999',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::partial_overlap_feature',
    },
    {
        SO_accession => 'SO:X000102',
        SO_term => 'entirely_within_feature',
        feature_SO_term => 'sequence_feature',
        feature_class => 'Bio::EnsEMBL::Feature',
		variant_feature_class => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
        rank => '999',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::complete_within_feature',
    },
    {
        SO_accession => 'SO:X000103',
        SO_term => 'complete_overlap_feature',
        feature_SO_term => 'sequence_feature',
        feature_class => 'Bio::EnsEMBL::Feature',
		variant_feature_class => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
        rank => '999',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::complete_overlap_feature',
    },
    {
        SO_accession => 'SO:0000159',
        SO_term => 'deletion',
        feature_SO_term => 'sequence_feature',
        feature_class => 'Bio::EnsEMBL::Feature',
		variant_feature_class => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
        rank => '999',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::deletion',
    },
    {
        SO_accession => 'SO:1000035',
        SO_term => 'duplication',
        feature_SO_term => 'sequence_feature',
        feature_class => 'Bio::EnsEMBL::Feature',
		variant_feature_class => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
        rank => '999',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::duplication',
    },
    {
        SO_accession => 'SO:1000173',
        SO_term => 'tandem_duplication',
        feature_SO_term => 'sequence_feature',
        feature_class => 'Bio::EnsEMBL::Feature',
		variant_feature_class => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
        rank => '999',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::tandem_duplication',
    },
    {
        SO_accession => 'SO:0000667',
        SO_term => 'insertion',
        feature_SO_term => 'sequence_feature',
        feature_class => 'Bio::EnsEMBL::Feature',
		variant_feature_class => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
        rank => '999',
        predicate => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::insertion',
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
        name => 'dbSNP clinical significance',
        description => 'The clinical significance of a variant as reported by dbSNP',
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

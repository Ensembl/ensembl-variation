package Bio::EnsEMBL::Variation::Utils::Constants;

#####################################################################
# NB: THIS FILE HAS BEEN AUTOMATICALLY GENERATED, EDIT WITH CAUTION #
#####################################################################

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(%OVERLAP_CONSEQUENCES %VARIATION_CLASSES $DEFAULT_OVERLAP_CONSEQUENCE $SO_ACC_MAPPER SO_TERM_TF_BINDING_SITE_VARIANT ATTRIB_TYPE_PROT_FUNC_ANALYSIS SO_TERM_START_LOST ATTRIB_TYPE_RANK ATTRIB_TYPE_DISPLAY_TERM SO_TERM_INTERGENIC_VARIANT SO_TERM_COPY_NUMBER_LOSS SO_TERM_ALU_INSERTION SO_TERM_ALU_DELETION SO_TERM_HERV_DELETION SO_TERM_HERV_INSERTION SO_TERM_LINE1_DELETION SO_TERM_LINE1_INSERTION SO_TERM_SVA_DELETION SO_TERM_SVA_INSERTION SO_TERM_COPY_NUMBER_GAIN SO_TERM_PROTEIN_ALTERING_VARIANT SO_TERM_SPLICE_REGION_VARIANT ATTRIB_TYPE_ALLELE_ACCESSION_ID ATTRIB_TYPE_EVIDENCE ATTRIB_TYPE_SEQUENCE_NUMBER SO_TERM_TANDEM_REPEAT ATTRIB_TYPE_MARKER_ACCESSION_ID ATTRIB_TYPE_VARIATION_NAMES SO_TERM_INTRACHROMOSOMAL_TRANSLOCATION ATTRIB_TYPE_SIFT_PREDICTION SO_TERM_INSERTION SO_TERM_REGULATORY_REGION_AMPLIFICATION ATTRIB_TYPE_LOD_SCORE SO_TERM_TRANSLOCATION SO_TERM_COMPLEX_SUBSTITUTION SO_TERM_SNV SO_TERM_SPLICE_DONOR_VARIANT SO_TERM_REGULATORY_REGION_VARIANT ATTRIB_TYPE_DBSNP_CLIN_SIG ATTRIB_TYPE_STRAIN_ID SO_TERM_INDEL SO_TERM_DELETION SO_TERM_NMD_TRANSCRIPT_VARIANT ATTRIB_TYPE_SO_TERM SO_TERM_INFRAME_DELETION SO_TERM_TFBS_ABLATION SO_TERM_MOBILE_ELEMENT_INSERTION SO_TERM_SUBSTITUTION SO_TERM_INTRON_VARIANT ATTRIB_TYPE_NCBI_TERM SO_TERM_MOBILE_ELEMENT_DELETION ATTRIB_TYPE_ODDS_RATIO ATTRIB_TYPE_SO_ACCESSION SO_TERM_INVERSION SO_TERM_COPY_NUMBER_VARIATION ATTRIB_TYPE_EXTERNAL_ID SO_TERM_FEATURE_TRUNCATION ATTRIB_TYPE_CONSERVATION_SCORE SO_TERM_COMPLEX_STRUCTURAL_ALTERATION SO_TERM_INTRACHROMOSOMAL_BREAKPOINT ATTRIB_TYPE_VARIANCE SO_TERM_NON_CODING_TRANSCRIPT_EXON_VARIANT SO_TERM_TFBS_AMPLIFICATION ATTRIB_TYPE_REVIEW_STATUS ATTRIB_TRAIT SO_TERM_STOP_RETAINED_VARIANT SO_TERM_TANDEM_DUPLICATION SO_TERM_STRUCTURAL_VARIANT SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT ATTRIB_TYPE_FEATURE_SO_TERM ATTRIB_TYPE_BASED_ON SO_TERM_DOWNSTREAM_GENE_VARIANT SO_TERM_SEQUENCE_VARIANT ATTRIB_TYPE_DGVA_CLIN_SIG SO_TERM_NON_CODING_TRANSCRIPT_VARIANT SO_TERM_5_PRIME_UTR_VARIANT SO_TERM_FRAMESHIFT_VARIANT SO_TERM_CODING_SEQUENCE_VARIANT ATTRIB_TYPE_SHORT_NAME SO_TERM_UPSTREAM_GENE_VARIANT SO_TERM_INFRAME_INSERTION SO_TERM_PROBE ATTRIB_TYPE_BETA_COEF ATTRIB_TYPE_CLINVAR_CLIN_SIG SO_TERM_SHORT_TANDEM_REPEAT_VARIATION SO_TERM_LOSS_OF_HETEROZYGOSITY SO_TERM_TRANSCRIPT_ABLATION SO_TERM_SYNONYMOUS_VARIANT SO_TERM_DUPLICATION SO_TERM_INTERCHROMOSOMAL_TRANSLOCATION SO_TERM_GENETIC_MARKER ATTRIB_TYPE_INHERITANCE_TYPE SO_TERM_SPLICE_ACCEPTOR_VARIANT SO_TERM_SEQUENCE_ALTERATION SO_TERM_MATURE_MIRNA_VARIANT ATTRIB_TYPE_ALLELE_SYMBOL SO_TERM_MISSENSE_VARIANT ATTRIB_TYPE_ASSOCIATED_GENE SO_TERM_TRANSCRIPT_AMPLIFICATION ATTRIB_TYPE_POLYPHEN_PREDICTION SO_TERM_NOVEL_SEQUENCE_INSERTION SO_TERM_3_PRIME_UTR_VARIANT SO_TERM_STOP_LOST SO_TERM_INTERCHROMOSOMAL_BREAKPOINT SO_TERM_START_RETAINED_VARIANT ATTRIB_TYPE_RISK_ALLELE SO_TERM_STOP_GAINED SO_TERM_REGULATORY_REGION_ABLATION ATTRIB_TYPE_P_VALUE SO_TERM_FEATURE_ELONGATION ATTRIB_TYPE_SAMPLE_ID ATTRIB_TYPE_PHENOTYPE_TYPE);

our %EXPORT_TAGS = ( attrib_types => [qw(ATTRIB_TYPE_VARIATION_NAMES ATTRIB_TYPE_LOD_SCORE ATTRIB_TYPE_SIFT_PREDICTION ATTRIB_TYPE_SHORT_NAME ATTRIB_TYPE_EVIDENCE ATTRIB_TYPE_SEQUENCE_NUMBER ATTRIB_TYPE_MARKER_ACCESSION_ID ATTRIB_TYPE_CLINVAR_CLIN_SIG ATTRIB_TYPE_DBSNP_CLIN_SIG ATTRIB_TYPE_STRAIN_ID ATTRIB_TYPE_BETA_COEF ATTRIB_TYPE_BASED_ON ATTRIB_TYPE_PROT_FUNC_ANALYSIS ATTRIB_TYPE_DISPLAY_TERM ATTRIB_TYPE_RANK ATTRIB_TYPE_FEATURE_SO_TERM ATTRIB_TYPE_ALLELE_ACCESSION_ID ATTRIB_TYPE_DGVA_CLIN_SIG ATTRIB_TYPE_POLYPHEN_PREDICTION ATTRIB_TYPE_CONSERVATION_SCORE ATTRIB_TYPE_EXTERNAL_ID ATTRIB_TYPE_P_VALUE ATTRIB_TYPE_SAMPLE_ID ATTRIB_TYPE_REVIEW_STATUS ATTRIB_TYPE_VARIANCE ATTRIB_TYPE_RISK_ALLELE ATTRIB_TYPE_INHERITANCE_TYPE ATTRIB_TYPE_NCBI_TERM ATTRIB_TYPE_SO_TERM ATTRIB_TYPE_SO_ACCESSION ATTRIB_TYPE_ASSOCIATED_GENE ATTRIB_TYPE_ALLELE_SYMBOL ATTRIB_TYPE_ODDS_RATIO ATTRIB_TYPE_PHENOTYPE_TYPE)], SO_consequence_terms => [qw(SO_TERM_REGULATORY_REGION_ABLATION SO_TERM_FEATURE_ELONGATION SO_TERM_TFBS_AMPLIFICATION SO_TERM_NON_CODING_TRANSCRIPT_EXON_VARIANT SO_TERM_START_RETAINED_VARIANT SO_TERM_STOP_RETAINED_VARIANT SO_TERM_STOP_GAINED SO_TERM_3_PRIME_UTR_VARIANT SO_TERM_STOP_LOST SO_TERM_FEATURE_TRUNCATION SO_TERM_TRANSCRIPT_AMPLIFICATION SO_TERM_MATURE_MIRNA_VARIANT SO_TERM_MISSENSE_VARIANT SO_TERM_SPLICE_ACCEPTOR_VARIANT SO_TERM_INTRON_VARIANT SO_TERM_NMD_TRANSCRIPT_VARIANT SO_TERM_TRANSCRIPT_ABLATION SO_TERM_INFRAME_DELETION SO_TERM_SYNONYMOUS_VARIANT SO_TERM_TFBS_ABLATION SO_TERM_REGULATORY_REGION_VARIANT SO_TERM_SPLICE_DONOR_VARIANT SO_TERM_UPSTREAM_GENE_VARIANT SO_TERM_REGULATORY_REGION_AMPLIFICATION SO_TERM_INFRAME_INSERTION SO_TERM_FRAMESHIFT_VARIANT SO_TERM_SPLICE_REGION_VARIANT SO_TERM_CODING_SEQUENCE_VARIANT SO_TERM_SEQUENCE_VARIANT SO_TERM_NON_CODING_TRANSCRIPT_VARIANT SO_TERM_5_PRIME_UTR_VARIANT SO_TERM_PROTEIN_ALTERING_VARIANT SO_TERM_INTERGENIC_VARIANT SO_TERM_DOWNSTREAM_GENE_VARIANT SO_TERM_TF_BINDING_SITE_VARIANT SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT SO_TERM_START_LOST)], SO_class_terms => [qw(SO_TERM_GENETIC_MARKER SO_TERM_INTERCHROMOSOMAL_TRANSLOCATION SO_TERM_DUPLICATION SO_TERM_MOBILE_ELEMENT_INSERTION SO_TERM_ALU_INSERTION SO_TERM_ALU_DELETION SO_TERM_HERV_DELETION SO_TERM_HERV_INSERTION SO_TERM_LINE1_DELETION SO_TERM_LINE1_INSERTION SO_TERM_SVA_DELETION SO_TERM_SVA_INSERTION SO_TERM_COPY_NUMBER_LOSS SO_TERM_SEQUENCE_ALTERATION SO_TERM_SUBSTITUTION SO_TERM_MOBILE_ELEMENT_DELETION SO_TERM_STRUCTURAL_VARIANT SO_TERM_TANDEM_DUPLICATION SO_TERM_INVERSION SO_TERM_COPY_NUMBER_VARIATION SO_TERM_COPY_NUMBER_GAIN SO_TERM_INTRACHROMOSOMAL_TRANSLOCATION SO_TERM_NOVEL_SEQUENCE_INSERTION SO_TERM_INSERTION SO_TERM_PROBE SO_TERM_TANDEM_REPEAT SO_TERM_LOSS_OF_HETEROZYGOSITY SO_TERM_SHORT_TANDEM_REPEAT_VARIATION SO_TERM_INDEL SO_TERM_DELETION SO_TERM_INTERCHROMOSOMAL_BREAKPOINT SO_TERM_INTRACHROMOSOMAL_BREAKPOINT SO_TERM_COMPLEX_STRUCTURAL_ALTERATION SO_TERM_TRANSLOCATION SO_TERM_COMPLEX_SUBSTITUTION SO_TERM_SNV)],  );

use Bio::EnsEMBL::Variation::OverlapConsequence;

use constant ATTRIB_TYPE_SO_ACCESSION => 'SO_accession';
use constant ATTRIB_TYPE_SO_TERM => 'SO_term';
use constant ATTRIB_TYPE_DISPLAY_TERM => 'display_term';
use constant ATTRIB_TYPE_NCBI_TERM => 'NCBI_term';
use constant ATTRIB_TYPE_FEATURE_SO_TERM => 'feature_SO_term';
use constant ATTRIB_TYPE_RANK => 'rank';
use constant ATTRIB_TYPE_POLYPHEN_PREDICTION => 'polyphen_prediction';
use constant ATTRIB_TYPE_SIFT_PREDICTION => 'sift_prediction';
use constant ATTRIB_TYPE_SHORT_NAME => 'short_name';
use constant ATTRIB_TYPE_DBSNP_CLIN_SIG => 'dbsnp_clin_sig';
use constant ATTRIB_TYPE_DGVA_CLIN_SIG => 'dgva_clin_sig';
use constant ATTRIB_TYPE_CLINVAR_CLIN_SIG => 'clinvar_clin_sig';
use constant ATTRIB_TYPE_PROT_FUNC_ANALYSIS => 'prot_func_analysis';
use constant ATTRIB_TYPE_ASSOCIATED_GENE => 'associated_gene';
use constant ATTRIB_TYPE_RISK_ALLELE => 'risk_allele';
use constant ATTRIB_TYPE_P_VALUE => 'p_value';
use constant ATTRIB_TYPE_VARIATION_NAMES => 'variation_names';
use constant ATTRIB_TYPE_SAMPLE_ID => 'sample_id';
use constant ATTRIB_TYPE_STRAIN_ID => 'strain_id';
use constant ATTRIB_TYPE_LOD_SCORE => 'lod_score';
use constant ATTRIB_TYPE_VARIANCE => 'variance';
use constant ATTRIB_TYPE_INHERITANCE_TYPE => 'inheritance_type';
use constant ATTRIB_TYPE_EXTERNAL_ID => 'external_id';
use constant ATTRIB_TYPE_ODDS_RATIO => 'odds_ratio';
use constant ATTRIB_TYPE_BETA_COEF => 'beta_coef';
use constant ATTRIB_TYPE_ALLELE_SYMBOL => 'allele_symbol';
use constant ATTRIB_TYPE_ALLELE_ACCESSION_ID => 'allele_accession_id';
use constant ATTRIB_TYPE_MARKER_ACCESSION_ID => 'marker_accession_id';
use constant ATTRIB_TYPE_EVIDENCE => 'evidence';
use constant ATTRIB_TYPE_SEQUENCE_NUMBER => 'sequence_number';
use constant ATTRIB_TYPE_BASED_ON => 'based_on';
use constant ATTRIB_TYPE_CONSERVATION_SCORE => 'conservation_score';
use constant ATTRIB_TYPE_REVIEW_STATUS => 'review_status';
use constant ATTRIB_TYPE_PHENOTYPE_TYPE => 'phenotype_type';

use constant ATTRIB_TRAIT => 'trait';

use constant SO_TERM_SNV => 'SNV';
use constant SO_TERM_SUBSTITUTION => 'substitution';
use constant SO_TERM_INSERTION => 'insertion';
use constant SO_TERM_DELETION => 'deletion';
use constant SO_TERM_INDEL => 'indel';
use constant SO_TERM_TANDEM_REPEAT => 'tandem_repeat';
use constant SO_TERM_SEQUENCE_ALTERATION => 'sequence_alteration';
use constant SO_TERM_GENETIC_MARKER => 'genetic_marker';
use constant SO_TERM_STRUCTURAL_VARIANT => 'structural_variant';
use constant SO_TERM_COPY_NUMBER_VARIATION => 'copy_number_variation';
use constant SO_TERM_PROBE => 'probe';
use constant SO_TERM_COPY_NUMBER_GAIN => 'copy_number_gain';
use constant SO_TERM_COPY_NUMBER_LOSS => 'copy_number_loss';
use constant SO_TERM_INVERSION => 'inversion';
use constant SO_TERM_COMPLEX_STRUCTURAL_ALTERATION => 'complex_structural_alteration';
use constant SO_TERM_TANDEM_DUPLICATION => 'tandem_duplication';
use constant SO_TERM_MOBILE_ELEMENT_INSERTION => 'mobile_element_insertion';
use constant SO_TERM_MOBILE_ELEMENT_DELETION => 'mobile_element_deletion';
use constant SO_TERM_INTERCHROMOSOMAL_BREAKPOINT => 'interchromosomal_breakpoint';
use constant SO_TERM_INTRACHROMOSOMAL_BREAKPOINT => 'intrachromosomal_breakpoint';
use constant SO_TERM_TRANSLOCATION => 'translocation';
use constant SO_TERM_DUPLICATION => 'duplication';
use constant SO_TERM_NOVEL_SEQUENCE_INSERTION => 'novel_sequence_insertion';
use constant SO_TERM_INTERCHROMOSOMAL_TRANSLOCATION => 'interchromosomal_translocation';
use constant SO_TERM_INTRACHROMOSOMAL_TRANSLOCATION => 'intrachromosomal_translocation';
use constant SO_TERM_ALU_INSERTION => 'Alu_insertion';
use constant SO_TERM_ALU_DELETION => 'Alu_deletion';
use constant SO_TERM_HERV_DELETION => 'HERV_deletion';
use constant SO_TERM_HERV_INSERTION => 'HERV_insertion';
use constant SO_TERM_LINE1_DELETION => 'LINE1_deletion';
use constant SO_TERM_LINE1_INSERTION => 'LINE1_insertion';
use constant SO_TERM_SVA_DELETION => 'SVA_deletion';
use constant SO_TERM_SVA_INSERTION => 'SVA_insertion';
use constant SO_TERM_COMPLEX_SUBSTITUTION => 'complex_substitution';
use constant SO_TERM_SHORT_TANDEM_REPEAT_VARIATION => 'short_tandem_repeat_variation';
use constant SO_TERM_LOSS_OF_HETEROZYGOSITY => 'loss_of_heterozygosity';
use constant SO_TERM_SEQUENCE_VARIANT => 'sequence_variant';
use constant SO_TERM_INTERGENIC_VARIANT => 'intergenic_variant';
use constant SO_TERM_UPSTREAM_GENE_VARIANT => 'upstream_gene_variant';
use constant SO_TERM_DOWNSTREAM_GENE_VARIANT => 'downstream_gene_variant';
use constant SO_TERM_SPLICE_DONOR_VARIANT => 'splice_donor_variant';
use constant SO_TERM_SPLICE_ACCEPTOR_VARIANT => 'splice_acceptor_variant';
use constant SO_TERM_SPLICE_REGION_VARIANT => 'splice_region_variant';
use constant SO_TERM_INTRON_VARIANT => 'intron_variant';
use constant SO_TERM_5_PRIME_UTR_VARIANT => '5_prime_UTR_variant';
use constant SO_TERM_3_PRIME_UTR_VARIANT => '3_prime_UTR_variant';
use constant SO_TERM_SYNONYMOUS_VARIANT => 'synonymous_variant';
use constant SO_TERM_MISSENSE_VARIANT => 'missense_variant';
use constant SO_TERM_INFRAME_INSERTION => 'inframe_insertion';
use constant SO_TERM_INFRAME_DELETION => 'inframe_deletion';
use constant SO_TERM_STOP_GAINED => 'stop_gained';
use constant SO_TERM_STOP_LOST => 'stop_lost';
use constant SO_TERM_STOP_RETAINED_VARIANT => 'stop_retained_variant';
use constant SO_TERM_START_LOST => 'start_lost';
use constant SO_TERM_START_RETAINED_VARIANT => 'start_retained_variant';
use constant SO_TERM_FRAMESHIFT_VARIANT => 'frameshift_variant';
use constant SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT => 'incomplete_terminal_codon_variant';
use constant SO_TERM_NMD_TRANSCRIPT_VARIANT => 'NMD_transcript_variant';
use constant SO_TERM_NON_CODING_TRANSCRIPT_VARIANT => 'non_coding_transcript_variant';
use constant SO_TERM_NON_CODING_TRANSCRIPT_EXON_VARIANT => 'non_coding_transcript_exon_variant';
use constant SO_TERM_MATURE_MIRNA_VARIANT => 'mature_miRNA_variant';
use constant SO_TERM_CODING_SEQUENCE_VARIANT => 'coding_sequence_variant';
use constant SO_TERM_REGULATORY_REGION_VARIANT => 'regulatory_region_variant';
use constant SO_TERM_TF_BINDING_SITE_VARIANT => 'TF_binding_site_variant';
use constant SO_TERM_TRANSCRIPT_ABLATION => 'transcript_ablation';
use constant SO_TERM_TRANSCRIPT_AMPLIFICATION => 'transcript_amplification';
use constant SO_TERM_TFBS_ABLATION => 'TFBS_ablation';
use constant SO_TERM_TFBS_AMPLIFICATION => 'TFBS_amplification';
use constant SO_TERM_REGULATORY_REGION_ABLATION => 'regulatory_region_ablation';
use constant SO_TERM_REGULATORY_REGION_AMPLIFICATION => 'regulatory_region_amplification';
use constant SO_TERM_FEATURE_ELONGATION => 'feature_elongation';
use constant SO_TERM_FEATURE_TRUNCATION => 'feature_truncation';
use constant SO_TERM_PROTEIN_ALTERING_VARIANT => 'protein_altering_variant';

our %VARIATION_CLASSES = (
'SNV' => {
  'SO_accession' => 'SO:0001483',
  'display_term' => 'SNP',
  'somatic_display_term' => 'somatic SNV'
}
,
'substitution' => {
  'SO_accession' => 'SO:1000002',
  'display_term' => 'substitution',
  'somatic_display_term' => 'somatic substitution'
}
,
'insertion' => {
  'SO_accession' => 'SO:0000667',
  'display_term' => 'insertion',
  'somatic_display_term' => 'somatic insertion'
}
,
'deletion' => {
  'SO_accession' => 'SO:0000159',
  'display_term' => 'deletion',
  'somatic_display_term' => 'somatic deletion'
}
,
'indel' => {
  'SO_accession' => 'SO:1000032',
  'display_term' => 'indel',
  'somatic_display_term' => 'somatic indel'
}
,
'tandem_repeat' => {
  'SO_accession' => 'SO:0000705',
  'display_term' => 'tandem repeat',
  'somatic_display_term' => 'somatic tandem repeat'
}
,
'sequence_alteration' => {
  'SO_accession' => 'SO:0001059',
  'display_term' => 'sequence alteration',
  'somatic_display_term' => 'somatic sequence alteration'
}
,
'genetic_marker' => {
  'SO_accession' => 'SO:0001645',
  'display_term' => 'genetic marker',
  'somatic_display_term' => 'somatic genetic marker'
}
,
'structural_variant' => {
  'SO_accession' => 'SO:0001537',
  'display_term' => 'SV',
  'somatic_display_term' => 'somatic SV',
  'type' => 'sv'
}
,
'copy_number_variation' => {
  'SO_accession' => 'SO:0001019',
  'display_term' => 'CNV',
  'somatic_display_term' => 'somatic CNV',
  'type' => 'sv'
}
,
'probe' => {
  'SO_accession' => 'SO:0000051',
  'display_term' => 'CNV_PROBE',
  'somatic_display_term' => 'somatic CNV_PROBE',
  'type' => 'sv'
}
,
'copy_number_gain' => {
  'SO_accession' => 'SO:0001742',
  'display_term' => 'gain',
  'somatic_display_term' => 'somatic gain',
  'type' => 'sv'
}
,
'copy_number_loss' => {
  'SO_accession' => 'SO:0001743',
  'display_term' => 'loss',
  'somatic_display_term' => 'somatic loss',
  'type' => 'sv'
}
,
'inversion' => {
  'SO_accession' => 'SO:1000036',
  'display_term' => 'inversion',
  'somatic_display_term' => 'somatic inversion',
  'type' => 'sv'
}
,
'complex_structural_alteration' => {
  'SO_accession' => 'SO:0001784',
  'display_term' => 'complex alteration',
  'somatic_display_term' => 'somatic complex alteration',
  'type' => 'sv'
}
,
'tandem_duplication' => {
  'SO_accession' => 'SO:1000173',
  'display_term' => 'tandem duplication',
  'somatic_display_term' => 'somatic tandem duplication',
  'type' => 'sv'
}
,
'mobile_element_insertion' => {
  'SO_accession' => 'SO:0001837',
  'display_term' => 'mobile element insertion',
  'somatic_display_term' => 'somatic mobile element insertion',
  'type' => 'sv'
}
,
'mobile_element_deletion' => {
  'SO_accession' => 'SO:0002066',
  'display_term' => 'mobile element deletion',
  'somatic_display_term' => 'somatic mobile element deletion',
  'type' => 'sv'
}
,
'interchromosomal_breakpoint' => {
  'SO_accession' => 'SO:0001873',
  'display_term' => 'interchromosomal breakpoint',
  'somatic_display_term' => 'somatic interchromosomal breakpoint',
  'type' => 'sv'
}
,
'intrachromosomal_breakpoint' => {
  'SO_accession' => 'SO:0001874',
  'display_term' => 'intrachromosomal breakpoint',
  'somatic_display_term' => 'somatic intrachromosomal breakpoint',
  'type' => 'sv'
}
,
'translocation' => {
  'SO_accession' => 'SO:0000199',
  'display_term' => 'translocation',
  'somatic_display_term' => 'somatic translocation',
  'type' => 'sv'
}
,
'duplication' => {
  'SO_accession' => 'SO:1000035',
  'display_term' => 'duplication',
  'somatic_display_term' => 'somatic duplication',
  'type' => 'sv'
}
,
'novel_sequence_insertion' => {
  'SO_accession' => 'SO:0001838',
  'display_term' => 'novel sequence insertion',
  'somatic_display_term' => 'somatic novel sequence insertion',
  'type' => 'sv'
}
,
'interchromosomal_translocation' => {
  'SO_accession' => 'SO:0002060',
  'display_term' => 'interchromosomal translocation',
  'somatic_display_term' => 'somatic interchromosomal translocation',
  'type' => 'sv'
}
,
'intrachromosomal_translocation' => {
  'SO_accession' => 'SO:0002061',
  'display_term' => 'intrachromosomal translocation',
  'somatic_display_term' => 'somatic intrachromosomal translocation',
  'type' => 'sv'
}
,
'Alu_insertion' => {
  'SO_accession' => 'SO:0002063',
  'display_term' => 'Alu insertion',
  'somatic_display_term' => 'somatic alu insertion',
  'type' => 'sv'
}
,
'Alu_deletion' => {
  'SO_accession' => 'SO:0002070',
  'display_term' => 'Alu deletion',
  'somatic_display_term' => 'somatic alu deletion',
  'type' => 'sv'
}
,
'HERV_deletion' => {
  'SO_accession' => 'SO:0002067',
  'display_term' => 'HERV deletion',
  'somatic_display_term' => 'somatic herv deletion',
  'type' => 'sv'
}
,
'HERV_insertion' => {
  'SO_accession' => 'SO:0002187',
  'display_term' => 'HERV insertion',
  'somatic_display_term' => 'somatic herv insertion',
  'type' => 'sv'
}
,
'LINE1_deletion' => {
  'SO_accession' => 'SO:0002069',
  'display_term' => 'LINE1 deletion',
  'somatic_display_term' => 'somatic line1 deletion',
  'type' => 'sv'
}
,
'LINE1_insertion' => {
  'SO_accession' => 'SO:0002064',
  'display_term' => 'LINE1 insertion',
  'somatic_display_term' => 'somatic line1 insertion',
  'type' => 'sv'
}
,
'SVA_deletion' => {
  'SO_accession' => 'SO:0002068',
  'display_term' => 'SVA deletion',
  'somatic_display_term' => 'somatic sva deletion',
  'type' => 'sv'
}
,
'SVA_insertion' => {
  'SO_accession' => 'SO:0002065',
  'display_term' => 'SVA insertion',
  'somatic_display_term' => 'somatic sva insertion',
  'type' => 'sv'
}
,
'complex_substitution' => {
  'SO_accession' => 'SO:1000005',
  'display_term' => 'complex substitution',
  'somatic_display_term' => 'somatic complex substitution',
  'type' => 'sv'
}
,
'short_tandem_repeat_variation' => {
  'SO_accession' => 'SO:0002096',
  'display_term' => 'short tandem repeat variation',
  'somatic_display_term' => 'somatic short tandem repeat variation',
  'type' => 'sv'
}
,
'loss_of_heterozygosity' => {
  'SO_accession' => 'SO:0001786',
  'display_term' => 'loss of heterozygosity',
  'somatic_display_term' => 'somatic loss of heterozygosity',
  'type' => 'sv'
}
,
);

our $DEFAULT_OVERLAP_CONSEQUENCE = Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001628',
  'SO_term' => 'intergenic_variant',
  'description' => 'A sequence variant located in the intergenic region, between genes',
  'display_term' => 'INTERGENIC',
  'impact' => 'MODIFIER',
  'include' => {
                 'within_feature' => 0
               },
  'is_default' => 1,
  'label' => 'intergenic variant',
  'rank' => '38',
  'tier' => '4'
}
);


our %OVERLAP_CONSEQUENCES = (
'sequence_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001060',
  'SO_term' => 'sequence_variant',
  'description' => 'A sequence_variant is a non exact copy of a sequence_feature or genome exhibiting one or more sequence_alteration',
  'display_term' => 'SEQUENCE_VARIANT',
  'impact' => 'MODIFIER',
  'include' => {
                 'within_feature' => 0
               },
  'label' => 'sequence variant',
  'rank' => '39',
  'tier' => '4'
}
),
'intergenic_variant' => $DEFAULT_OVERLAP_CONSEQUENCE,
'upstream_gene_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001631',
  'SO_term' => 'upstream_gene_variant',
  'description' => 'A sequence variant located 5\' of a gene',
  'display_term' => 'UPSTREAM',
  'feature_SO_term' => 'transcript',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODIFIER',
  'include' => {
                 'within_feature' => 0
               },
  'label' => 'upstream gene variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream',
  'rank' => '24',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'downstream_gene_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001632',
  'SO_term' => 'downstream_gene_variant',
  'description' => 'A sequence variant located 3\' of a gene',
  'display_term' => 'DOWNSTREAM',
  'feature_SO_term' => 'transcript',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODIFIER',
  'include' => {
                 'within_feature' => 0
               },
  'label' => 'downstream gene variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream',
  'rank' => '25',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'splice_donor_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'NCBI_term' => 'splice-5',
  'SO_accession' => 'SO:0001575',
  'SO_term' => 'splice_donor_variant',
  'description' => 'A splice variant that changes the 2 base region at the 5\' end of an intron',
  'display_term' => 'ESSENTIAL_SPLICE_SITE',
  'feature_SO_term' => 'primary_transcript',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'HIGH',
  'include' => {
                 'intron_boundary' => 1
               },
  'label' => 'splice donor variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::donor_splice_site',
  'rank' => '3',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature'
}
),
'splice_acceptor_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'NCBI_term' => 'splice-3',
  'SO_accession' => 'SO:0001574',
  'SO_term' => 'splice_acceptor_variant',
  'description' => 'A splice variant that changes the 2 base region at the 3\' end of an intron',
  'display_term' => 'ESSENTIAL_SPLICE_SITE',
  'feature_SO_term' => 'primary_transcript',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'HIGH',
  'include' => {
                 'intron_boundary' => 1
               },
  'label' => 'splice acceptor variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::acceptor_splice_site',
  'rank' => '3',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature'
}
),
'splice_region_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001630',
  'SO_term' => 'splice_region_variant',
  'description' => 'A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron',
  'display_term' => 'SPLICE_SITE',
  'feature_SO_term' => 'primary_transcript',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'LOW',
  'include' => {
                 'intron_boundary' => 1
               },
  'label' => 'splice region variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::splice_region',
  'rank' => '13',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature'
}
),
'intron_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'NCBI_term' => 'intron',
  'SO_accession' => 'SO:0001627',
  'SO_term' => 'intron_variant',
  'description' => 'A transcript variant occurring within an intron',
  'display_term' => 'INTRONIC',
  'feature_SO_term' => 'primary_transcript',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODIFIER',
  'include' => {
                 'intron' => 1
               },
  'label' => 'intron variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_intron',
  'rank' => '21',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'5_prime_UTR_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'NCBI_term' => 'untranslated_5',
  'SO_accession' => 'SO:0001623',
  'SO_term' => '5_prime_UTR_variant',
  'description' => 'A UTR variant of the 5\' UTR',
  'display_term' => '5PRIME_UTR',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODIFIER',
  'include' => {
                 'exon' => 1,
                 'utr' => 1
               },
  'label' => '5 prime UTR variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_5_prime_utr',
  'rank' => '18',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'3_prime_UTR_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'NCBI_term' => 'untranslated_3',
  'SO_accession' => 'SO:0001624',
  'SO_term' => '3_prime_UTR_variant',
  'description' => 'A UTR variant of the 3\' UTR',
  'display_term' => '3PRIME_UTR',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODIFIER',
  'include' => {
                 'exon' => 1,
                 'utr' => 1
               },
  'label' => '3 prime UTR variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_3_prime_utr',
  'rank' => '19',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'synonymous_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'NCBI_term' => 'cds-synon',
  'SO_accession' => 'SO:0001819',
  'SO_term' => 'synonymous_variant',
  'description' => 'A sequence variant where there is no resulting change to the encoded amino acid',
  'display_term' => 'SYNONYMOUS_CODING',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'LOW',
  'include' => {
                 'coding' => 1
               },
  'label' => 'synonymous variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::synonymous_variant',
  'rank' => '15',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature'
}
),
'missense_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'NCBI_term' => 'missense',
  'SO_accession' => 'SO:0001583',
  'SO_term' => 'missense_variant',
  'description' => 'A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODERATE',
  'include' => {
                 'coding' => 1,
                 'decrease_length' => 0,
                 'increase_length' => 0
               },
  'label' => 'missense variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::missense_variant',
  'rank' => '12',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature'
}
),
'inframe_insertion' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001821',
  'SO_term' => 'inframe_insertion',
  'description' => 'An inframe non synonymous variant that inserts bases into in the coding sequence',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODERATE',
  'include' => {
                 'coding' => 1,
                 'insertion' => 1
               },
  'label' => 'inframe insertion',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_insertion',
  'rank' => '10',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'inframe_deletion' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001822',
  'SO_term' => 'inframe_deletion',
  'description' => 'An inframe non synonymous variant that deletes bases from the coding sequence',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODERATE',
  'include' => {
                 'coding' => 1,
                 'deletion' => 1
               },
  'label' => 'inframe deletion',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_deletion',
  'rank' => '11',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'stop_gained' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'NCBI_term' => 'nonsense',
  'SO_accession' => 'SO:0001587',
  'SO_term' => 'stop_gained',
  'description' => 'A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript',
  'display_term' => 'STOP_GAINED',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'HIGH',
  'include' => {
                 'coding' => 1
               },
  'label' => 'stop gained',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_gained',
  'rank' => '4',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature'
}
),
'stop_lost' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001578',
  'SO_term' => 'stop_lost',
  'description' => 'A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript',
  'display_term' => 'STOP_LOST',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'HIGH',
  'include' => {
                 'coding' => 1
               },
  'label' => 'stop lost',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_lost',
  'rank' => '6',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'stop_retained_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001567',
  'SO_term' => 'stop_retained_variant',
  'description' => 'A sequence variant where at least one base in the terminator codon is changed, but the terminator remains',
  'display_term' => 'SYNONYMOUS_CODING',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'LOW',
  'include' => {
                 'coding' => 1
               },
  'label' => 'stop retained variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_retained',
  'rank' => '15',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature'
}
),
'start_lost' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0002012',
  'SO_term' => 'start_lost',
  'description' => 'A codon variant that changes at least one base of the canonical start codon',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'HIGH',
  'include' => {
                 'coding' => 1
               },
  'label' => 'start lost',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::start_lost',
  'rank' => '7',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'start_retained_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0002019',
  'SO_term' => 'start_retained_variant',
  'description' => 'A sequence variant where at least one base in the start codon is changed, but the start remains',
  'display_term' => 'SYNONYMOUS_CODING',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'LOW',
  'include' => {
                 'coding' => 1
               },
  'label' => 'start retained variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::start_retained_variant',
  'rank' => '15',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'frameshift_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'NCBI_term' => 'frameshift',
  'SO_accession' => 'SO:0001589',
  'SO_term' => 'frameshift_variant',
  'description' => 'A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three',
  'display_term' => 'FRAMESHIFT_CODING',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'HIGH',
  'include' => {
                 'coding' => 1,
                 'snp' => 0
               },
  'label' => 'frameshift variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::frameshift',
  'rank' => '5',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'incomplete_terminal_codon_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001626',
  'SO_term' => 'incomplete_terminal_codon_variant',
  'description' => 'A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed',
  'display_term' => 'PARTIAL_CODON',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'LOW',
  'include' => {
                 'coding' => 1
               },
  'label' => 'incomplete terminal codon variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::partial_codon',
  'rank' => '14',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature'
}
),
'NMD_transcript_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001621',
  'SO_term' => 'NMD_transcript_variant',
  'description' => 'A variant in a transcript that is the target of NMD',
  'display_term' => 'NMD_TRANSCRIPT',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODIFIER',
  'include' => {
                 'nonsense_mediated_decay' => 1,
                 'within_feature' => 1
               },
  'label' => 'NMD transcript variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_nmd_transcript',
  'rank' => '22',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'non_coding_transcript_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001619',
  'SO_term' => 'non_coding_transcript_variant',
  'description' => 'A transcript variant of a non coding RNA gene',
  'display_term' => 'WITHIN_NON_CODING_GENE',
  'feature_SO_term' => 'ncRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODIFIER',
  'include' => {
                 'protein_coding' => 0,
                 'within_feature' => 1
               },
  'label' => 'non coding transcript variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_non_coding_gene',
  'rank' => '23',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'non_coding_transcript_exon_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001792',
  'SO_term' => 'non_coding_transcript_exon_variant',
  'description' => 'A sequence variant that changes non-coding exon sequence in a non-coding transcript',
  'display_term' => 'WITHIN_NON_CODING_GENE',
  'feature_SO_term' => 'ncRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODIFIER',
  'include' => {
                 'exon' => 1,
                 'protein_coding' => 0,
                 'within_feature' => 1
               },
  'label' => 'non coding transcript exon variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::non_coding_exon_variant',
  'rank' => '20',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'mature_miRNA_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001620',
  'SO_term' => 'mature_miRNA_variant',
  'description' => 'A transcript variant located with the sequence of the mature miRNA',
  'display_term' => 'WITHIN_MATURE_miRNA',
  'feature_SO_term' => 'miRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODIFIER',
  'include' => {
                 'nonsense_mediated_decay' => 0,
                 'protein_coding' => 0,
                 'within_feature' => 1
               },
  'label' => 'mature miRNA variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_mature_miRNA',
  'rank' => '17',
  'tier' => '2',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'coding_sequence_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001580',
  'SO_term' => 'coding_sequence_variant',
  'description' => 'A sequence variant that changes the coding sequence',
  'display_term' => 'CODING_UNKNOWN',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODIFIER',
  'include' => {
                 'coding' => 1
               },
  'label' => 'coding sequence variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::coding_unknown',
  'rank' => '16',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'regulatory_region_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001566',
  'SO_term' => 'regulatory_region_variant',
  'description' => 'A sequence variant located within a regulatory region',
  'display_term' => 'REGULATORY_REGION',
  'feature_SO_term' => 'regulatory_region',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
  'impact' => 'MODIFIER',
  'label' => 'regulatory region variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_regulatory_feature',
  'rank' => '36',
  'tier' => '2',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'TF_binding_site_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001782',
  'SO_term' => 'TF_binding_site_variant',
  'description' => 'A sequence variant located within a transcription factor binding site',
  'display_term' => 'REGULATORY_REGION',
  'feature_SO_term' => 'TF_binding_site',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::MotifFeature',
  'impact' => 'MODIFIER',
  'label' => 'TF binding site',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_motif_feature',
  'rank' => '30',
  'tier' => '2',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'transcript_ablation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001893',
  'SO_term' => 'transcript_ablation',
  'description' => 'A feature ablation whereby the deleted region includes a transcript feature',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'HIGH',
  'include' => {
                 'complete_overlap' => 1,
                 'deletion' => 1
               },
  'label' => 'transcript ablation',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
  'rank' => '1',
  'tier' => '1',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'transcript_amplification' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001889',
  'SO_term' => 'transcript_amplification',
  'description' => 'A feature amplification of a region containing a transcript',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'HIGH',
  'include' => {
                 'complete_overlap' => 1,
                 'increase_length' => 1
               },
  'label' => 'transcript amplification',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
  'rank' => '8',
  'tier' => '1',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'TFBS_ablation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001895',
  'SO_term' => 'TFBS_ablation',
  'description' => 'A feature ablation whereby the deleted region includes a transcription factor binding site',
  'feature_SO_term' => 'TF_binding_site',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::MotifFeature',
  'impact' => 'MODERATE',
  'include' => {
                 'complete_overlap' => 1,
                 'deletion' => 1
               },
  'label' => 'TFBS ablation',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
  'rank' => '26',
  'tier' => '2',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'TFBS_amplification' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001892',
  'SO_term' => 'TFBS_amplification',
  'description' => 'A feature amplification of a region containing a transcription factor binding site',
  'feature_SO_term' => 'TF_binding_site',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::MotifFeature',
  'impact' => 'MODIFIER',
  'include' => {
                 'complete_overlap' => 1,
                 'increase_length' => 1
               },
  'label' => 'TFBS amplification',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
  'rank' => '28',
  'tier' => '2',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'regulatory_region_ablation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001894',
  'SO_term' => 'regulatory_region_ablation',
  'description' => 'A feature ablation whereby the deleted region includes a regulatory region',
  'feature_SO_term' => 'TF_binding_site',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
  'impact' => 'MODIFIER',
  'include' => {
                 'complete_overlap' => 1,
                 'deletion' => 1
               },
  'label' => 'regulatory region ablation',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
  'rank' => '31',
  'tier' => '2',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'regulatory_region_amplification' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001891',
  'SO_term' => 'regulatory_region_amplification',
  'description' => 'A feature amplification of a region containing a regulatory region',
  'feature_SO_term' => 'TF_binding_site',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
  'impact' => 'MODIFIER',
  'include' => {
                 'complete_overlap' => 1,
                 'increase_length' => 1
               },
  'label' => 'regulatory region amplification',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
  'rank' => '33',
  'tier' => '2',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'feature_elongation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001907',
  'SO_term' => 'feature_elongation',
  'description' => 'A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence',
  'feature_SO_term' => 'sequence_feature',
  'feature_class' => 'Bio::EnsEMBL::Feature',
  'impact' => 'MODIFIER',
  'include' => {
                 'increase_length' => 1,
                 'sv' => 1
               },
  'label' => 'feature elongation',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_elongation',
  'rank' => '36',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'feature_truncation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001906',
  'SO_term' => 'feature_truncation',
  'description' => 'A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence',
  'feature_SO_term' => 'sequence_feature',
  'feature_class' => 'Bio::EnsEMBL::Feature',
  'impact' => 'MODIFIER',
  'include' => {
                 'decrease_length' => 1,
                 'sv' => 1
               },
  'label' => 'feature truncation',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_truncation',
  'rank' => '37',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature'
}
),
'protein_altering_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_accession' => 'SO:0001818',
  'SO_term' => 'protein_altering_variant',
  'description' => 'A sequence_variant which is predicted to change the protein encoded in the coding sequence',
  'feature_SO_term' => 'mRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript',
  'impact' => 'MODERATE',
  'include' => {
                 'coding' => 1
               },
  'label' => 'protein altering variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::protein_altering_variant',
  'rank' => '12',
  'tier' => '3',
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature'
}
),
);

our $SO_ACC_MAPPER = {
  'Bio::EnsEMBL::Variation::StructuralVariationFeature' => {
                                                             'acc' => 'SO:0001537',
                                                             'term' => 'structural_variant'
                                                           },
  'Bio::EnsEMBL::Variation::VariationFeature' => {
                                                   'acc' => 'SO:0001060',
                                                   'term' => 'sequence_variant'
                                                 }
}
;

1;

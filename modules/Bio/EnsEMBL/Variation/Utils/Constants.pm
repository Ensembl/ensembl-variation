package Bio::EnsEMBL::Variation::Utils::Constants;

#####################################################################
# NB: THIS FILE HAS BEEN AUTOMATICALLY GENERATED, EDIT WITH CAUTION #
#####################################################################

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(%OVERLAP_CONSEQUENCES %VARIATION_CLASSES $DEFAULT_OVERLAP_CONSEQUENCE SO_TERM_INDEL ATTRIB_TYPE_ALLELE_ACCESSION_ID SO_TERM_CODING_SEQUENCE_VARIANT SO_TERM_STOP_RETAINED_VARIANT SO_TERM_COMPLEX_STRUCTURAL_ALTERATION SO_TERM_ALU_INSERTION SO_TERM_SPLICE_ACCEPTOR_VARIANT SO_TERM_PROTEIN_ALTERING_VARIANT SO_TERM_TANDEM_REPEAT SO_TERM_INTERCHROMOSOMAL_BREAKPOINT ATTRIB_TYPE_VARIATION_NAMES SO_TERM_FEATURE_ELONGATION SO_TERM_TRANSCRIPT_ABLATION SO_TERM_TF_BINDING_SITE_VARIANT ATTRIB_TYPE_DISPLAY_TERM SO_TERM_INTRON_VARIANT SO_TERM_STOP_GAINED ATTRIB_TYPE_REVIEW_STATUS SO_TERM_SPLICE_DONOR_VARIANT ATTRIB_TYPE_SIFT_PREDICTION SO_TERM_START_LOST SO_TERM_FEATURE_TRUNCATION ATTRIB_TYPE_VARIANCE SO_TERM_COMPLEX_SUBSTITUTION ATTRIB_TYPE_BASED_ON SO_TERM_INTRACHROMOSOMAL_BREAKPOINT SO_TERM_SPLICE_REGION_VARIANT SO_TERM_TFBS_ABLATION ATTRIB_TYPE_ODDS_RATIO SO_TERM_REGULATORY_REGION_AMPLIFICATION SO_TERM_MISSENSE_VARIANT SO_TERM_SNV ATTRIB_TYPE_EXTERNAL_ID SO_TERM_STRUCTURAL_VARIANT SO_TERM_PROBE ATTRIB_TYPE_DBSNP_CLIN_SIG ATTRIB_TYPE_LOD_SCORE SO_TERM_SEQUENCE_ALTERATION SO_TERM_REGULATORY_REGION_VARIANT SO_TERM_START_RETAINED_VARIANT ATTRIB_TYPE_MARKER_ACCESSION_ID SO_TERM_SUBSTITUTION ATTRIB_TYPE_DGVA_CLIN_SIG SO_TERM_NON_CODING_TRANSCRIPT_EXON_VARIANT SO_TERM_NOVEL_SEQUENCE_INSERTION SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT SO_TERM_DUPLICATION SO_TERM_TFBS_AMPLIFICATION SO_TERM_5_PRIME_UTR_VARIANT ATTRIB_TYPE_CONSERVATION_SCORE ATTRIB_TYPE_PROT_FUNC_ANALYSIS ATTRIB_TYPE_STRAIN_ID SO_TERM_TANDEM_DUPLICATION ATTRIB_TYPE_P_VALUE ATTRIB_TYPE_ASSOCIATED_GENE SO_TERM_MOBILE_ELEMENT_INSERTION SO_TERM_MATURE_MIRNA_VARIANT SO_TERM_NON_CODING_TRANSCRIPT_VARIANT SO_TERM_DOWNSTREAM_GENE_VARIANT SO_TERM_INFRAME_INSERTION SO_TERM_INSERTION SO_TERM_SHORT_TANDEM_REPEAT_VARIATION SO_TERM_NMD_TRANSCRIPT_VARIANT ATTRIB_TYPE_SO_TERM ATTRIB_TYPE_SAMPLE_ID SO_TERM_INTERGENIC_VARIANT SO_TERM_SYNONYMOUS_VARIANT ATTRIB_TYPE_SHORT_NAME SO_TERM_TRANSLOCATION ATTRIB_TYPE_INHERITANCE_TYPE SO_TERM_COPY_NUMBER_VARIATION ATTRIB_TYPE_POLYPHEN_PREDICTION SO_TERM_REGULATORY_REGION_ABLATION SO_TERM_INVERSION SO_TERM_LOSS_OF_HETEROZYGOSITY SO_TERM_STOP_LOST ATTRIB_TYPE_RANK SO_TERM_INTERCHROMOSOMAL_TRANSLOCATION SO_TERM_3_PRIME_UTR_VARIANT SO_TERM_UPSTREAM_GENE_VARIANT SO_TERM_INTRACHROMOSOMAL_TRANSLOCATION ATTRIB_TYPE_CLINVAR_CLIN_SIG ATTRIB_TYPE_SO_ACCESSION ATTRIB_TYPE_ALLELE_SYMBOL SO_TERM_COPY_NUMBER_LOSS SO_TERM_MOBILE_ELEMENT_DELETION ATTRIB_TYPE_SEQUENCE_NUMBER ATTRIB_TYPE_NCBI_TERM SO_TERM_DELETION ATTRIB_TYPE_FEATURE_SO_TERM SO_TERM_FRAMESHIFT_VARIANT SO_TERM_COPY_NUMBER_GAIN ATTRIB_TYPE_RISK_ALLELE ATTRIB_TYPE_BETA_COEF ATTRIB_TYPE_EVIDENCE SO_TERM_INFRAME_DELETION SO_TERM_GENETIC_MARKER SO_TERM_TRANSCRIPT_AMPLIFICATION);

our %EXPORT_TAGS = ( attrib_types => [qw(ATTRIB_TYPE_ALLELE_ACCESSION_ID ATTRIB_TYPE_SO_TERM ATTRIB_TYPE_SAMPLE_ID ATTRIB_TYPE_VARIATION_NAMES ATTRIB_TYPE_SHORT_NAME ATTRIB_TYPE_DISPLAY_TERM ATTRIB_TYPE_INHERITANCE_TYPE ATTRIB_TYPE_REVIEW_STATUS ATTRIB_TYPE_SIFT_PREDICTION ATTRIB_TYPE_VARIANCE ATTRIB_TYPE_POLYPHEN_PREDICTION ATTRIB_TYPE_BASED_ON ATTRIB_TYPE_RANK ATTRIB_TYPE_ODDS_RATIO ATTRIB_TYPE_EXTERNAL_ID ATTRIB_TYPE_SO_ACCESSION ATTRIB_TYPE_CLINVAR_CLIN_SIG ATTRIB_TYPE_DBSNP_CLIN_SIG ATTRIB_TYPE_LOD_SCORE ATTRIB_TYPE_ALLELE_SYMBOL ATTRIB_TYPE_MARKER_ACCESSION_ID ATTRIB_TYPE_SEQUENCE_NUMBER ATTRIB_TYPE_DGVA_CLIN_SIG ATTRIB_TYPE_NCBI_TERM ATTRIB_TYPE_FEATURE_SO_TERM ATTRIB_TYPE_RISK_ALLELE ATTRIB_TYPE_CONSERVATION_SCORE ATTRIB_TYPE_BETA_COEF ATTRIB_TYPE_PROT_FUNC_ANALYSIS ATTRIB_TYPE_STRAIN_ID ATTRIB_TYPE_EVIDENCE ATTRIB_TYPE_P_VALUE ATTRIB_TYPE_ASSOCIATED_GENE)], SO_consequence_terms => [qw(SO_TERM_NON_CODING_TRANSCRIPT_VARIANT SO_TERM_DOWNSTREAM_GENE_VARIANT SO_TERM_CODING_SEQUENCE_VARIANT SO_TERM_STOP_RETAINED_VARIANT SO_TERM_INFRAME_INSERTION SO_TERM_NMD_TRANSCRIPT_VARIANT SO_TERM_PROTEIN_ALTERING_VARIANT SO_TERM_SPLICE_ACCEPTOR_VARIANT SO_TERM_INTERGENIC_VARIANT SO_TERM_FEATURE_ELONGATION SO_TERM_TRANSCRIPT_ABLATION SO_TERM_SYNONYMOUS_VARIANT SO_TERM_TF_BINDING_SITE_VARIANT SO_TERM_INTRON_VARIANT SO_TERM_STOP_GAINED SO_TERM_SPLICE_DONOR_VARIANT SO_TERM_FEATURE_TRUNCATION SO_TERM_REGULATORY_REGION_ABLATION SO_TERM_START_LOST SO_TERM_TFBS_ABLATION SO_TERM_SPLICE_REGION_VARIANT SO_TERM_STOP_LOST SO_TERM_REGULATORY_REGION_AMPLIFICATION SO_TERM_3_PRIME_UTR_VARIANT SO_TERM_UPSTREAM_GENE_VARIANT SO_TERM_MISSENSE_VARIANT SO_TERM_REGULATORY_REGION_VARIANT SO_TERM_START_RETAINED_VARIANT SO_TERM_NON_CODING_TRANSCRIPT_EXON_VARIANT SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT SO_TERM_FRAMESHIFT_VARIANT SO_TERM_TFBS_AMPLIFICATION SO_TERM_5_PRIME_UTR_VARIANT SO_TERM_INFRAME_DELETION SO_TERM_TRANSCRIPT_AMPLIFICATION SO_TERM_MATURE_MIRNA_VARIANT)], SO_class_terms => [qw(SO_TERM_INTERCHROMOSOMAL_TRANSLOCATION SO_TERM_INDEL SO_TERM_SNV SO_TERM_STRUCTURAL_VARIANT SO_TERM_ALU_INSERTION SO_TERM_INTRACHROMOSOMAL_TRANSLOCATION SO_TERM_COMPLEX_STRUCTURAL_ALTERATION SO_TERM_INSERTION SO_TERM_SHORT_TANDEM_REPEAT_VARIATION SO_TERM_PROBE SO_TERM_COPY_NUMBER_LOSS SO_TERM_SEQUENCE_ALTERATION SO_TERM_TANDEM_REPEAT SO_TERM_INTERCHROMOSOMAL_BREAKPOINT SO_TERM_MOBILE_ELEMENT_DELETION SO_TERM_SUBSTITUTION SO_TERM_NOVEL_SEQUENCE_INSERTION SO_TERM_DELETION SO_TERM_COPY_NUMBER_GAIN SO_TERM_DUPLICATION SO_TERM_TRANSLOCATION SO_TERM_COPY_NUMBER_VARIATION SO_TERM_TANDEM_DUPLICATION SO_TERM_GENETIC_MARKER SO_TERM_COMPLEX_SUBSTITUTION SO_TERM_INVERSION SO_TERM_LOSS_OF_HETEROZYGOSITY SO_TERM_INTRACHROMOSOMAL_BREAKPOINT SO_TERM_MOBILE_ELEMENT_INSERTION)],  );

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
use constant SO_TERM_COMPLEX_SUBSTITUTION => 'complex_substitution';
use constant SO_TERM_SHORT_TANDEM_REPEAT_VARIATION => 'short_tandem_repeat_variation';
use constant SO_TERM_LOSS_OF_HETEROZYGOSITY => 'loss_of_heterozygosity';
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
  'somatic_display_term' => 'somatic SNV',
  'SO_accession' => 'SO:0001483',
  'display_term' => 'SNP'
}
,
'substitution' => {
  'somatic_display_term' => 'somatic substitution',
  'SO_accession' => 'SO:1000002',
  'display_term' => 'substitution'
}
,
'insertion' => {
  'somatic_display_term' => 'somatic insertion',
  'SO_accession' => 'SO:0000667',
  'display_term' => 'insertion'
}
,
'deletion' => {
  'somatic_display_term' => 'somatic deletion',
  'SO_accession' => 'SO:0000159',
  'display_term' => 'deletion'
}
,
'indel' => {
  'somatic_display_term' => 'somatic indel',
  'SO_accession' => 'SO:1000032',
  'display_term' => 'indel'
}
,
'tandem_repeat' => {
  'somatic_display_term' => 'somatic tandem repeat',
  'SO_accession' => 'SO:0000705',
  'display_term' => 'tandem repeat'
}
,
'sequence_alteration' => {
  'somatic_display_term' => 'somatic sequence alteration',
  'SO_accession' => 'SO:0001059',
  'display_term' => 'sequence alteration'
}
,
'genetic_marker' => {
  'somatic_display_term' => 'somatic genetic marker',
  'SO_accession' => 'SO:0001645',
  'display_term' => 'genetic marker'
}
,
'structural_variant' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic SV',
  'SO_accession' => 'SO:0001537',
  'display_term' => 'SV'
}
,
'copy_number_variation' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic CNV',
  'SO_accession' => 'SO:0001019',
  'display_term' => 'CNV'
}
,
'probe' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic CNV_PROBE',
  'SO_accession' => 'SO:0000051',
  'display_term' => 'CNV_PROBE'
}
,
'copy_number_gain' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic gain',
  'SO_accession' => 'SO:0001742',
  'display_term' => 'gain'
}
,
'copy_number_loss' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic loss',
  'SO_accession' => 'SO:0001743',
  'display_term' => 'loss'
}
,
'inversion' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic inversion',
  'SO_accession' => 'SO:1000036',
  'display_term' => 'inversion'
}
,
'complex_structural_alteration' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic complex alteration',
  'SO_accession' => 'SO:0001784',
  'display_term' => 'complex alteration'
}
,
'tandem_duplication' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic tandem duplication',
  'SO_accession' => 'SO:1000173',
  'display_term' => 'tandem duplication'
}
,
'mobile_element_insertion' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic mobile element insertion',
  'SO_accession' => 'SO:0001837',
  'display_term' => 'mobile element insertion'
}
,
'mobile_element_deletion' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic mobile element deletion',
  'SO_accession' => 'SO:0002066',
  'display_term' => 'mobile element deletion'
}
,
'interchromosomal_breakpoint' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic interchromosomal breakpoint',
  'SO_accession' => 'SO:0001873',
  'display_term' => 'interchromosomal breakpoint'
}
,
'intrachromosomal_breakpoint' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic intrachromosomal breakpoint',
  'SO_accession' => 'SO:0001874',
  'display_term' => 'intrachromosomal breakpoint'
}
,
'translocation' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic translocation',
  'SO_accession' => 'SO:0000199',
  'display_term' => 'translocation'
}
,
'duplication' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic duplication',
  'SO_accession' => 'SO:1000035',
  'display_term' => 'duplication'
}
,
'novel_sequence_insertion' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic novel sequence insertion',
  'SO_accession' => 'SO:0001838',
  'display_term' => 'novel sequence insertion'
}
,
'interchromosomal_translocation' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic interchromosomal translocation',
  'SO_accession' => 'SO:0002060',
  'display_term' => 'interchromosomal translocation'
}
,
'intrachromosomal_translocation' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic intrachromosomal translocation',
  'SO_accession' => 'SO:0002061',
  'display_term' => 'intrachromosomal translocation'
}
,
'Alu_insertion' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic alu insertion',
  'SO_accession' => 'SO:0002063',
  'display_term' => 'Alu insertion'
}
,
'complex_substitution' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic complex substitution',
  'SO_accession' => 'SO:1000005',
  'display_term' => 'complex substitution'
}
,
'short_tandem_repeat_variation' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic short tandem repeat variation',
  'SO_accession' => 'SO:0002096',
  'display_term' => 'short tandem repeat variation'
}
,
'loss_of_heterozygosity' => {
  'type' => 'sv',
  'somatic_display_term' => 'somatic loss of heterozygosity',
  'SO_accession' => 'SO:0001786',
  'display_term' => 'loss of heterozygosity'
}
,
);

our $DEFAULT_OVERLAP_CONSEQUENCE = Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'is_default' => 1,
  'include' => {
                 'within_feature' => 0
               },
  'description' => 'A sequence variant located in the intergenic region, between genes',
  'SO_accession' => 'SO:0001628',
  'SO_term' => 'intergenic_variant',
  'tier' => '4',
  'label' => 'intergenic variant',
  'rank' => '38',
  'impact' => 'MODIFIER',
  'display_term' => 'INTERGENIC'
}
);


our %OVERLAP_CONSEQUENCES = (
'intergenic_variant' => $DEFAULT_OVERLAP_CONSEQUENCE,
'sequence_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'include' => {
                 'within_feature' => 0
               },
  'description' => 'A sequence_variant is a non exact copy of a sequence_feature or genome exhibiting one or more sequence_alteration',
  'SO_accession' => 'SO:0001060',
  'SO_term' => 'sequence_variant',
  'tier' => '4',
  'label' => 'sequence variant',
  'rank' => '39',
  'impact' => 'MODIFIER',
  'display_term' => 'SEQUENCE_VARIANT',
}
),
'upstream_gene_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'within_feature' => 0
               },
  'feature_SO_term' => 'transcript',
  'description' => 'A sequence variant located 5\' of a gene',
  'SO_accession' => 'SO:0001631',
  'SO_term' => 'upstream_gene_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream',
  'label' => 'upstream gene variant',
  'rank' => '24',
  'impact' => 'MODIFIER',
  'display_term' => 'UPSTREAM',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'downstream_gene_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'within_feature' => 0
               },
  'feature_SO_term' => 'transcript',
  'description' => 'A sequence variant located 3\' of a gene',
  'SO_accession' => 'SO:0001632',
  'SO_term' => 'downstream_gene_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream',
  'label' => 'downstream gene variant',
  'rank' => '25',
  'impact' => 'MODIFIER',
  'display_term' => 'DOWNSTREAM',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'splice_donor_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'include' => {
                 'intron_boundary' => 1
               },
  'NCBI_term' => 'splice-5',
  'feature_SO_term' => 'primary_transcript',
  'description' => 'A splice variant that changes the 2 base region at the 5\' end of an intron',
  'SO_accession' => 'SO:0001575',
  'tier' => '3',
  'SO_term' => 'splice_donor_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::donor_splice_site',
  'label' => 'splice donor variant',
  'rank' => '3',
  'impact' => 'HIGH',
  'display_term' => 'ESSENTIAL_SPLICE_SITE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'splice_acceptor_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'include' => {
                 'intron_boundary' => 1
               },
  'NCBI_term' => 'splice-3',
  'feature_SO_term' => 'primary_transcript',
  'description' => 'A splice variant that changes the 2 base region at the 3\' end of an intron',
  'SO_accession' => 'SO:0001574',
  'tier' => '3',
  'SO_term' => 'splice_acceptor_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::acceptor_splice_site',
  'label' => 'splice acceptor variant',
  'rank' => '3',
  'impact' => 'HIGH',
  'display_term' => 'ESSENTIAL_SPLICE_SITE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'splice_region_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'include' => {
                 'intron_boundary' => 1
               },
  'feature_SO_term' => 'primary_transcript',
  'description' => 'A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron',
  'SO_accession' => 'SO:0001630',
  'SO_term' => 'splice_region_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::splice_region',
  'label' => 'splice region variant',
  'rank' => '13',
  'impact' => 'LOW',
  'display_term' => 'SPLICE_SITE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'intron_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'intron' => 1
               },
  'NCBI_term' => 'intron',
  'feature_SO_term' => 'primary_transcript',
  'description' => 'A transcript variant occurring within an intron',
  'SO_accession' => 'SO:0001627',
  'tier' => '3',
  'SO_term' => 'intron_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_intron',
  'label' => 'intron variant',
  'rank' => '21',
  'impact' => 'MODIFIER',
  'display_term' => 'INTRONIC',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'5_prime_UTR_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'utr' => 1,
                 'exon' => 1
               },
  'NCBI_term' => 'untranslated_5',
  'feature_SO_term' => 'mRNA',
  'description' => 'A UTR variant of the 5\' UTR',
  'SO_accession' => 'SO:0001623',
  'tier' => '3',
  'SO_term' => '5_prime_UTR_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_5_prime_utr',
  'label' => '5 prime UTR variant',
  'rank' => '18',
  'impact' => 'MODIFIER',
  'display_term' => '5PRIME_UTR',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'3_prime_UTR_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'utr' => 1,
                 'exon' => 1
               },
  'NCBI_term' => 'untranslated_3',
  'feature_SO_term' => 'mRNA',
  'description' => 'A UTR variant of the 3\' UTR',
  'SO_accession' => 'SO:0001624',
  'tier' => '3',
  'SO_term' => '3_prime_UTR_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_3_prime_utr',
  'label' => '3 prime UTR variant',
  'rank' => '19',
  'impact' => 'MODIFIER',
  'display_term' => '3PRIME_UTR',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'synonymous_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'include' => {
                 'coding' => 1
               },
  'NCBI_term' => 'cds-synon',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant where there is no resulting change to the encoded amino acid',
  'SO_accession' => 'SO:0001819',
  'tier' => '3',
  'SO_term' => 'synonymous_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::synonymous_variant',
  'label' => 'synonymous variant',
  'rank' => '15',
  'impact' => 'LOW',
  'display_term' => 'SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'missense_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'include' => {
                 'increase_length' => 0,
                 'decrease_length' => 0,
                 'coding' => 1
               },
  'NCBI_term' => 'missense',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved',
  'SO_accession' => 'SO:0001583',
  'tier' => '3',
  'SO_term' => 'missense_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::missense_variant',
  'label' => 'missense variant',
  'rank' => '12',
  'impact' => 'MODERATE',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'inframe_insertion' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'insertion' => 1,
                 'coding' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'An inframe non synonymous variant that inserts bases into in the coding sequence',
  'SO_accession' => 'SO:0001821',
  'SO_term' => 'inframe_insertion',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_insertion',
  'label' => 'inframe insertion',
  'rank' => '10',
  'impact' => 'MODERATE',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'inframe_deletion' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'deletion' => 1,
                 'coding' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'An inframe non synonymous variant that deletes bases from the coding sequence',
  'SO_accession' => 'SO:0001822',
  'SO_term' => 'inframe_deletion',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_deletion',
  'label' => 'inframe deletion',
  'rank' => '11',
  'impact' => 'MODERATE',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'stop_gained' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'include' => {
                 'coding' => 1
               },
  'NCBI_term' => 'nonsense',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript',
  'SO_accession' => 'SO:0001587',
  'tier' => '3',
  'SO_term' => 'stop_gained',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_gained',
  'label' => 'stop gained',
  'rank' => '4',
  'impact' => 'HIGH',
  'display_term' => 'STOP_GAINED',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'stop_lost' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'coding' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript',
  'SO_accession' => 'SO:0001578',
  'SO_term' => 'stop_lost',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_lost',
  'label' => 'stop lost',
  'rank' => '6',
  'impact' => 'HIGH',
  'display_term' => 'STOP_LOST',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'stop_retained_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'include' => {
                 'coding' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant where at least one base in the terminator codon is changed, but the terminator remains',
  'SO_accession' => 'SO:0001567',
  'SO_term' => 'stop_retained_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_retained',
  'label' => 'stop retained variant',
  'rank' => '15',
  'impact' => 'LOW',
  'display_term' => 'SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'start_lost' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'coding' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'A codon variant that changes at least one base of the canonical start codon',
  'SO_accession' => 'SO:0002012',
  'SO_term' => 'start_lost',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::start_lost',
  'label' => 'start lost',
  'rank' => '7',
  'impact' => 'HIGH',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'start_retained_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'coding' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant where at least one base in the start codon is changed, but the start remains',
  'SO_accession' => 'SO:0002019',
  'SO_term' => 'start_retained_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::start_retained_variant',
  'label' => 'start retained variant',
  'rank' => '15',
  'impact' => 'LOW',
  'display_term' => 'SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'frameshift_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'snp' => 0,
                 'coding' => 1
               },
  'NCBI_term' => 'frameshift',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three',
  'SO_accession' => 'SO:0001589',
  'tier' => '3',
  'SO_term' => 'frameshift_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::frameshift',
  'label' => 'frameshift variant',
  'rank' => '5',
  'impact' => 'HIGH',
  'display_term' => 'FRAMESHIFT_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'incomplete_terminal_codon_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'include' => {
                 'coding' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed',
  'SO_accession' => 'SO:0001626',
  'SO_term' => 'incomplete_terminal_codon_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::partial_codon',
  'label' => 'incomplete terminal codon variant',
  'rank' => '14',
  'impact' => 'LOW',
  'display_term' => 'PARTIAL_CODON',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'NMD_transcript_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'within_feature' => 1,
                 'nonsense_mediated_decay' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'A variant in a transcript that is the target of NMD',
  'SO_accession' => 'SO:0001621',
  'SO_term' => 'NMD_transcript_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_nmd_transcript',
  'label' => 'NMD transcript variant',
  'rank' => '22',
  'impact' => 'MODIFIER',
  'display_term' => 'NMD_TRANSCRIPT',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'non_coding_transcript_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'within_feature' => 1,
                 'protein_coding' => 0
               },
  'feature_SO_term' => 'ncRNA',
  'description' => 'A transcript variant of a non coding RNA gene',
  'SO_accession' => 'SO:0001619',
  'SO_term' => 'non_coding_transcript_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_non_coding_gene',
  'label' => 'non coding transcript variant',
  'rank' => '23',
  'impact' => 'MODIFIER',
  'display_term' => 'WITHIN_NON_CODING_GENE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'non_coding_transcript_exon_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'within_feature' => 1,
                 'exon' => 1,
                 'protein_coding' => 0
               },
  'feature_SO_term' => 'ncRNA',
  'description' => 'A sequence variant that changes non-coding exon sequence in a non-coding transcript',
  'SO_accession' => 'SO:0001792',
  'SO_term' => 'non_coding_transcript_exon_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::non_coding_exon_variant',
  'label' => 'non coding transcript exon variant',
  'rank' => '20',
  'impact' => 'MODIFIER',
  'display_term' => 'WITHIN_NON_CODING_GENE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'mature_miRNA_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'within_feature' => 1,
                 'nonsense_mediated_decay' => 0,
                 'protein_coding' => 0
               },
  'feature_SO_term' => 'miRNA',
  'description' => 'A transcript variant located with the sequence of the mature miRNA',
  'SO_accession' => 'SO:0001620',
  'SO_term' => 'mature_miRNA_variant',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_mature_miRNA',
  'label' => 'mature miRNA variant',
  'rank' => '17',
  'impact' => 'MODIFIER',
  'display_term' => 'WITHIN_MATURE_miRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'coding_sequence_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'coding' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant that changes the coding sequence',
  'SO_accession' => 'SO:0001580',
  'SO_term' => 'coding_sequence_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::coding_unknown',
  'label' => 'coding sequence variant',
  'rank' => '16',
  'impact' => 'MODIFIER',
  'display_term' => 'CODING_UNKNOWN',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'regulatory_region_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'regulatory_region',
  'description' => 'A sequence variant located within a regulatory region',
  'SO_accession' => 'SO:0001566',
  'SO_term' => 'regulatory_region_variant',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_regulatory_feature',
  'label' => 'regulatory region variant',
  'rank' => '36',
  'impact' => 'MODIFIER',
  'display_term' => 'REGULATORY_REGION',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature'
}
),
'TF_binding_site_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'TF_binding_site',
  'description' => 'A sequence variant located within a transcription factor binding site',
  'SO_accession' => 'SO:0001782',
  'SO_term' => 'TF_binding_site_variant',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_motif_feature',
  'label' => 'TF binding site',
  'rank' => '30',
  'impact' => 'MODIFIER',
  'display_term' => 'REGULATORY_REGION',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::MotifFeature'
}
),
'transcript_ablation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'complete_overlap' => 1,
                 'deletion' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'A feature ablation whereby the deleted region includes a transcript feature',
  'SO_accession' => 'SO:0001893',
  'SO_term' => 'transcript_ablation',
  'tier' => '1',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
  'label' => 'transcript ablation',
  'rank' => '1',
  'impact' => 'HIGH',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'transcript_amplification' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'complete_overlap' => 1,
                 'increase_length' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'A feature amplification of a region containing a transcript',
  'SO_accession' => 'SO:0001889',
  'SO_term' => 'transcript_amplification',
  'tier' => '1',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
  'label' => 'transcript amplification',
  'rank' => '8',
  'impact' => 'HIGH',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'TFBS_ablation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'complete_overlap' => 1,
                 'deletion' => 1
               },
  'feature_SO_term' => 'TF_binding_site',
  'description' => 'A feature ablation whereby the deleted region includes a transcription factor binding site',
  'SO_accession' => 'SO:0001895',
  'SO_term' => 'TFBS_ablation',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
  'label' => 'TFBS ablation',
  'rank' => '26',
  'impact' => 'MODERATE',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::MotifFeature'
}
),
'TFBS_amplification' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'complete_overlap' => 1,
                 'increase_length' => 1
               },
  'feature_SO_term' => 'TF_binding_site',
  'description' => 'A feature amplification of a region containing a transcription factor binding site',
  'SO_accession' => 'SO:0001892',
  'SO_term' => 'TFBS_amplification',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
  'label' => 'TFBS amplification',
  'rank' => '28',
  'impact' => 'MODIFIER',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::MotifFeature'
}
),
'regulatory_region_ablation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'complete_overlap' => 1,
                 'deletion' => 1
               },
  'feature_SO_term' => 'TF_binding_site',
  'description' => 'A feature ablation whereby the deleted region includes a regulatory region',
  'SO_accession' => 'SO:0001894',
  'SO_term' => 'regulatory_region_ablation',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
  'label' => 'regulatory region ablation',
  'rank' => '31',
  'impact' => 'MODERATE',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature'
}
),
'regulatory_region_amplification' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'complete_overlap' => 1,
                 'increase_length' => 1
               },
  'feature_SO_term' => 'TF_binding_site',
  'description' => 'A feature amplification of a region containing a regulatory region',
  'SO_accession' => 'SO:0001891',
  'SO_term' => 'regulatory_region_amplification',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
  'label' => 'regulatory region amplification',
  'rank' => '33',
  'impact' => 'MODIFIER',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature'
}
),
'feature_elongation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'sv' => 1,
                 'increase_length' => 1
               },
  'feature_SO_term' => 'sequence_feature',
  'description' => 'A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence',
  'SO_accession' => 'SO:0001907',
  'SO_term' => 'feature_elongation',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_elongation',
  'label' => 'feature elongation',
  'rank' => '36',
  'impact' => 'MODIFIER',
  'feature_class' => 'Bio::EnsEMBL::Feature'
}
),
'feature_truncation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'include' => {
                 'sv' => 1,
                 'decrease_length' => 1
               },
  'feature_SO_term' => 'sequence_feature',
  'description' => 'A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence',
  'SO_accession' => 'SO:0001906',
  'SO_term' => 'feature_truncation',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_truncation',
  'label' => 'feature truncation',
  'rank' => '37',
  'impact' => 'MODIFIER',
  'feature_class' => 'Bio::EnsEMBL::Feature'
}
),
'protein_altering_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'include' => {
                 'coding' => 1
               },
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence_variant which is predicted to change the protein encoded in the coding sequence',
  'SO_accession' => 'SO:0001818',
  'SO_term' => 'protein_altering_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::protein_altering_variant',
  'label' => 'protein altering variant',
  'rank' => '12',
  'impact' => 'MODERATE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
);

1;

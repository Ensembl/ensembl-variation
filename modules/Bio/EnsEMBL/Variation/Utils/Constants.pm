package Bio::EnsEMBL::Variation::Utils::Constants;

#####################################################################
# NB: THIS FILE HAS BEEN AUTOMATICALLY GENERATED, EDIT WITH CAUTION #
#####################################################################

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(%OVERLAP_CONSEQUENCES %VARIATION_CLASSES $DEFAULT_OVERLAP_CONSEQUENCE SO_TERM_NC_TRANSCRIPT_VARIANT SO_TERM_INDEL SO_TERM_PARTIAL_OVERLAP SO_TERM_CODING_SEQUENCE_VARIANT SO_TERM_STOP_RETAINED_VARIANT SO_TERM_COMPLEX_STRUCTURAL_ALTERATION SO_TERM_INFRAME_CODON_GAIN SO_TERM_INSERTION SO_TERM_NMD_TRANSCRIPT_VARIANT SO_TERM_SPLICE_ACCEPTOR_VARIANT SO_TERM_TANDEM_REPEAT ATTRIB_TYPE_SO_TERM SO_TERM_INTERGENIC_VARIANT SO_TERM_COMPLETE_OVERLAP_FEATURE SO_TERM_TF_BINDING_SITE_VARIANT ATTRIB_TYPE_SHORT_NAME ATTRIB_TYPE_DISPLAY_TERM SO_TERM_INTRON_VARIANT SO_TERM_ENTIRELY_WITHIN_FEATURE SO_TERM_STOP_GAINED SO_TERM_NON_SYNONYMOUS_CODON SO_TERM_COPY_NUMBER_VARIATION SO_TERM_SPLICE_DONOR_VARIANT ATTRIB_TYPE_SIFT_PREDICTION ATTRIB_TYPE_POLYPHEN_PREDICTION SO_TERM_INVERSION SO_TERM_500B_DOWNSTREAM_VARIANT SO_TERM_SPLICE_REGION_VARIANT SO_TERM_STOP_LOST SO_TERM_INITIATOR_CODON_CHANGE ATTRIB_TYPE_RANK SO_TERM_3_PRIME_UTR_VARIANT SO_TERM_SNV SO_TERM_INFRAME_CODON_LOSS SO_TERM_STRUCTURAL_VARIANT ATTRIB_TYPE_SO_ACCESSION SO_TERM_PROBE ATTRIB_TYPE_DBSNP_CLIN_SIG SO_TERM_SEQUENCE_ALTERATION SO_TERM_COPY_NUMBER_LOSS SO_TERM_REGULATORY_REGION_VARIANT SO_TERM_5KB_UPSTREAM_VARIANT SO_TERM_SUBSTITUTION ATTRIB_TYPE_DGVA_CLIN_SIG SO_TERM_2KB_UPSTREAM_VARIANT ATTRIB_TYPE_NCBI_TERM SO_TERM_DELETION SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT SO_TERM_COMPLEX_CHANGE_IN_TRANSCRIPT SO_TERM_SYNONYMOUS_CODON SO_TERM_FRAMESHIFT_VARIANT ATTRIB_TYPE_FEATURE_SO_TERM SO_TERM_COPY_NUMBER_GAIN SO_TERM_DUPLICATION SO_TERM_5_PRIME_UTR_VARIANT SO_TERM_TANDEM_DUPLICATION SO_TERM_5KB_DOWNSTREAM_VARIANT SO_TERM_MOBILE_ELEMENT_INSERTION SO_TERM_MATURE_MIRNA_VARIANT);

our %EXPORT_TAGS = ( attrib_types => [qw(ATTRIB_TYPE_RANK ATTRIB_TYPE_NCBI_TERM ATTRIB_TYPE_SHORT_NAME ATTRIB_TYPE_FEATURE_SO_TERM ATTRIB_TYPE_DISPLAY_TERM ATTRIB_TYPE_SO_ACCESSION ATTRIB_TYPE_DBSNP_CLIN_SIG ATTRIB_TYPE_SO_TERM ATTRIB_TYPE_SIFT_PREDICTION ATTRIB_TYPE_POLYPHEN_PREDICTION ATTRIB_TYPE_DGVA_CLIN_SIG)], SO_consequence_terms => [qw(SO_TERM_NC_TRANSCRIPT_VARIANT SO_TERM_PARTIAL_OVERLAP SO_TERM_CODING_SEQUENCE_VARIANT SO_TERM_STOP_RETAINED_VARIANT SO_TERM_INSERTION SO_TERM_INFRAME_CODON_GAIN SO_TERM_NMD_TRANSCRIPT_VARIANT SO_TERM_SPLICE_ACCEPTOR_VARIANT SO_TERM_COMPLETE_OVERLAP_FEATURE SO_TERM_INTERGENIC_VARIANT SO_TERM_TF_BINDING_SITE_VARIANT SO_TERM_ENTIRELY_WITHIN_FEATURE SO_TERM_INTRON_VARIANT SO_TERM_NON_SYNONYMOUS_CODON SO_TERM_STOP_GAINED SO_TERM_SPLICE_DONOR_VARIANT SO_TERM_500B_DOWNSTREAM_VARIANT SO_TERM_INITIATOR_CODON_CHANGE SO_TERM_STOP_LOST SO_TERM_SPLICE_REGION_VARIANT SO_TERM_3_PRIME_UTR_VARIANT SO_TERM_INFRAME_CODON_LOSS SO_TERM_5KB_UPSTREAM_VARIANT SO_TERM_REGULATORY_REGION_VARIANT SO_TERM_2KB_UPSTREAM_VARIANT SO_TERM_COMPLEX_CHANGE_IN_TRANSCRIPT SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT SO_TERM_DELETION SO_TERM_SYNONYMOUS_CODON SO_TERM_FRAMESHIFT_VARIANT SO_TERM_DUPLICATION SO_TERM_5_PRIME_UTR_VARIANT SO_TERM_TANDEM_DUPLICATION SO_TERM_5KB_DOWNSTREAM_VARIANT SO_TERM_MATURE_MIRNA_VARIANT)], SO_class_terms => [qw(SO_TERM_INDEL SO_TERM_DELETION SO_TERM_SNV SO_TERM_COPY_NUMBER_GAIN SO_TERM_STRUCTURAL_VARIANT SO_TERM_COMPLEX_STRUCTURAL_ALTERATION SO_TERM_INSERTION SO_TERM_PROBE SO_TERM_COPY_NUMBER_LOSS SO_TERM_SEQUENCE_ALTERATION SO_TERM_COPY_NUMBER_VARIATION SO_TERM_TANDEM_DUPLICATION SO_TERM_TANDEM_REPEAT SO_TERM_SUBSTITUTION SO_TERM_INVERSION SO_TERM_MOBILE_ELEMENT_INSERTION)],  );

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

use constant SO_TERM_SNV => 'SNV';
use constant SO_TERM_SUBSTITUTION => 'substitution';
use constant SO_TERM_COPY_NUMBER_VARIATION => 'copy_number_variation';
use constant SO_TERM_INSERTION => 'insertion';
use constant SO_TERM_DELETION => 'deletion';
use constant SO_TERM_INDEL => 'indel';
use constant SO_TERM_TANDEM_REPEAT => 'tandem_repeat';
use constant SO_TERM_SEQUENCE_ALTERATION => 'sequence_alteration';
use constant SO_TERM_STRUCTURAL_VARIANT => 'structural_variant';
use constant SO_TERM_PROBE => 'probe';
use constant SO_TERM_COPY_NUMBER_GAIN => 'copy_number_gain';
use constant SO_TERM_COPY_NUMBER_LOSS => 'copy_number_loss';
use constant SO_TERM_INVERSION => 'inversion';
use constant SO_TERM_COMPLEX_STRUCTURAL_ALTERATION => 'complex_structural_alteration';
use constant SO_TERM_TANDEM_DUPLICATION => 'tandem_duplication';
use constant SO_TERM_MOBILE_ELEMENT_INSERTION => 'mobile_element_insertion';
use constant SO_TERM_INTERGENIC_VARIANT => 'intergenic_variant';
use constant SO_TERM_5KB_UPSTREAM_VARIANT => '5KB_upstream_variant';
use constant SO_TERM_5KB_DOWNSTREAM_VARIANT => '5KB_downstream_variant';
use constant SO_TERM_2KB_UPSTREAM_VARIANT => '2KB_upstream_variant';
use constant SO_TERM_500B_DOWNSTREAM_VARIANT => '500B_downstream_variant';
use constant SO_TERM_SPLICE_DONOR_VARIANT => 'splice_donor_variant';
use constant SO_TERM_SPLICE_ACCEPTOR_VARIANT => 'splice_acceptor_variant';
use constant SO_TERM_SPLICE_REGION_VARIANT => 'splice_region_variant';
use constant SO_TERM_INTRON_VARIANT => 'intron_variant';
use constant SO_TERM_5_PRIME_UTR_VARIANT => '5_prime_UTR_variant';
use constant SO_TERM_3_PRIME_UTR_VARIANT => '3_prime_UTR_variant';
use constant SO_TERM_COMPLEX_CHANGE_IN_TRANSCRIPT => 'complex_change_in_transcript';
use constant SO_TERM_SYNONYMOUS_CODON => 'synonymous_codon';
use constant SO_TERM_NON_SYNONYMOUS_CODON => 'non_synonymous_codon';
use constant SO_TERM_INFRAME_CODON_GAIN => 'inframe_codon_gain';
use constant SO_TERM_INFRAME_CODON_LOSS => 'inframe_codon_loss';
use constant SO_TERM_STOP_GAINED => 'stop_gained';
use constant SO_TERM_STOP_LOST => 'stop_lost';
use constant SO_TERM_STOP_RETAINED_VARIANT => 'stop_retained_variant';
use constant SO_TERM_INITIATOR_CODON_CHANGE => 'initiator_codon_change';
use constant SO_TERM_FRAMESHIFT_VARIANT => 'frameshift_variant';
use constant SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT => 'incomplete_terminal_codon_variant';
use constant SO_TERM_NMD_TRANSCRIPT_VARIANT => 'NMD_transcript_variant';
use constant SO_TERM_NC_TRANSCRIPT_VARIANT => 'nc_transcript_variant';
use constant SO_TERM_MATURE_MIRNA_VARIANT => 'mature_miRNA_variant';
use constant SO_TERM_CODING_SEQUENCE_VARIANT => 'coding_sequence_variant';
use constant SO_TERM_REGULATORY_REGION_VARIANT => 'regulatory_region_variant';
use constant SO_TERM_TF_BINDING_SITE_VARIANT => 'TF_binding_site_variant';
use constant SO_TERM_PARTIAL_OVERLAP => 'partial_overlap';
use constant SO_TERM_ENTIRELY_WITHIN_FEATURE => 'entirely_within_feature';
use constant SO_TERM_COMPLETE_OVERLAP_FEATURE => 'complete_overlap_feature';
use constant SO_TERM_DUPLICATION => 'duplication';

our %VARIATION_CLASSES = (
'SNV' => {
  'somatic_display_term' => 'somatic_SNV',
  'SO_accession' => 'SO:0001483',
  'display_term' => 'SNP'
}
,
'substitution' => {
  'somatic_display_term' => 'somatic_substitution',
  'SO_accession' => 'SO:1000002',
  'display_term' => 'substitution'
}
,
'copy_number_variation' => {
  'somatic_display_term' => 'somatic_CNV',
  'SO_accession' => 'SO:0001019',
  'display_term' => 'CNV'
}
,
'insertion' => {
  'somatic_display_term' => 'somatic_insertion',
  'SO_accession' => 'SO:0000667',
  'display_term' => 'insertion'
}
,
'deletion' => {
  'somatic_display_term' => 'somatic_deletion',
  'SO_accession' => 'SO:0000159',
  'display_term' => 'deletion'
}
,
'indel' => {
  'somatic_display_term' => 'somatic_indel',
  'SO_accession' => 'SO:1000032',
  'display_term' => 'indel'
}
,
'tandem_repeat' => {
  'somatic_display_term' => 'somatic_tandem_repeat',
  'SO_accession' => 'SO:0000705',
  'display_term' => 'tandem_repeat'
}
,
'sequence_alteration' => {
  'somatic_display_term' => 'somatic_sequence_alteration',
  'SO_accession' => 'SO:0001059',
  'display_term' => 'sequence_alteration'
}
,
'structural_variant' => {
  'somatic_display_term' => 'somatic_SV',
  'SO_accession' => 'SO:0001537',
  'display_term' => 'SV'
}
,
'probe' => {
  'somatic_display_term' => 'somatic_CNV_PROBE',
  'SO_accession' => 'SO:0000051',
  'display_term' => 'CNV_PROBE'
}
,
'copy_number_gain' => {
  'somatic_display_term' => 'somatic_Gain',
  'SO_accession' => 'SO:0001742',
  'display_term' => 'Gain'
}
,
'copy_number_loss' => {
  'somatic_display_term' => 'somatic_Loss',
  'SO_accession' => 'SO:0001743',
  'display_term' => 'Loss'
}
,
'inversion' => {
  'somatic_display_term' => 'somatic_inversion',
  'SO_accession' => 'SO:1000036',
  'display_term' => 'inversion'
}
,
'complex_structural_alteration' => {
  'somatic_display_term' => 'somatic_Complex',
  'SO_accession' => 'SO:0001784',
  'display_term' => 'Complex'
}
,
'tandem_duplication' => {
  'somatic_display_term' => 'somatic_Tandem duplication',
  'SO_accession' => 'SO:1000173',
  'display_term' => 'Tandem duplication'
}
,
'mobile_element_insertion' => {
  'somatic_display_term' => 'somatic_Mobile element insertion',
  'SO_accession' => 'SO:0001837',
  'display_term' => 'Mobile element insertion'
}
,
);

our $DEFAULT_OVERLAP_CONSEQUENCE = Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'SO_term' => 'intergenic_variant',
  'is_default' => 1,
  'label' => 'Intergenic',
  'description' => 'More than 5 kb either upstream or downstream of a transcript',
  'rank' => '26',
  'SO_accession' => 'SO:0001628',
  'display_term' => 'INTERGENIC'
}
);


our %OVERLAP_CONSEQUENCES = (
'intergenic_variant' => $DEFAULT_OVERLAP_CONSEQUENCE,
'5KB_upstream_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'transcript',
  'description' => 'Within 5 kb upstream of the 5 prime end of a transcript',
  'SO_accession' => 'SO:0001635',
  'SO_term' => '5KB_upstream_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream_5KB',
  'label' => 'Upstream',
  'rank' => '22',
  'display_term' => 'UPSTREAM',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'5KB_downstream_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'transcript',
  'description' => 'Within 5 kb downstream of the 3 prime end of a transcript',
  'SO_accession' => 'SO:0001633',
  'SO_term' => '5KB_downstream_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream_5KB',
  'label' => 'Downstream',
  'rank' => '23',
  'display_term' => 'DOWNSTREAM',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'2KB_upstream_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'near-gene-5',
  'feature_SO_term' => 'transcript',
  'description' => 'Within 5 kb upstream of the 5 prime end of a transcript',
  'SO_accession' => 'SO:0001636',
  'SO_term' => '2KB_upstream_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream_2KB',
  'label' => 'Upstream',
  'rank' => '20',
  'display_term' => 'UPSTREAM',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'500B_downstream_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'near-gene-3',
  'feature_SO_term' => 'transcript',
  'description' => 'Within 5 kb downstream of the 3 prime end of a transcript',
  'SO_accession' => 'SO:0001634',
  'SO_term' => '500B_downstream_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream_500B',
  'label' => 'Downstream',
  'rank' => '21',
  'display_term' => 'DOWNSTREAM',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'splice_donor_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'splice-5',
  'feature_SO_term' => 'primary_transcript',
  'description' => 'In the first 2 or the last 2 basepairs of an intron',
  'SO_accession' => 'SO:0001575',
  'SO_term' => 'splice_donor_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::donor_splice_site',
  'label' => 'Essential splice site',
  'rank' => '1',
  'display_term' => 'ESSENTIAL_SPLICE_SITE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'splice_acceptor_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'splice-3',
  'feature_SO_term' => 'primary_transcript',
  'description' => 'In the first 2 or the last 2 basepairs of an intron',
  'SO_accession' => 'SO:0001574',
  'SO_term' => 'splice_acceptor_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::acceptor_splice_site',
  'label' => 'Essential splice site',
  'rank' => '1',
  'display_term' => 'ESSENTIAL_SPLICE_SITE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'splice_region_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'primary_transcript',
  'description' => '1-3 bps into an exon or 3-8 bps into an intron',
  'SO_accession' => 'SO:0001630',
  'SO_term' => 'splice_region_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::splice_region',
  'label' => 'Splice site',
  'rank' => '10',
  'display_term' => 'SPLICE_SITE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'intron_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'intron',
  'feature_SO_term' => 'primary_transcript',
  'description' => 'In intron',
  'SO_accession' => 'SO:0001627',
  'SO_term' => 'intron_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_intron',
  'label' => 'Intronic',
  'rank' => '17',
  'display_term' => 'INTRONIC',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'5_prime_UTR_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'untranslated_5',
  'feature_SO_term' => 'mRNA',
  'description' => 'In 5 prime untranslated region',
  'SO_accession' => 'SO:0001623',
  'SO_term' => '5_prime_UTR_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_5_prime_utr',
  'label' => '5 prime UTR',
  'rank' => '15',
  'display_term' => '5PRIME_UTR',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'3_prime_UTR_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'untranslated_3',
  'feature_SO_term' => 'mRNA',
  'description' => 'In 3 prime untranslated region',
  'SO_accession' => 'SO:0001624',
  'SO_term' => '3_prime_UTR_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_3_prime_utr',
  'label' => '3 prime UTR',
  'rank' => '16',
  'display_term' => '3PRIME_UTR',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'complex_change_in_transcript' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'primary_transcript',
  'description' => 'Insertion or deletion that spans an exon/intron or coding sequence/UTR border',
  'SO_accession' => 'SO:0001577',
  'SO_term' => 'complex_change_in_transcript',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::complex_indel',
  'label' => 'Complex in/del',
  'rank' => '4',
  'display_term' => 'COMPLEX_INDEL',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'synonymous_codon' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'cds-synon',
  'feature_SO_term' => 'mRNA',
  'description' => 'In coding sequence, not resulting in an amino acid change (silent mutation)',
  'SO_accession' => 'SO:0001588',
  'SO_term' => 'synonymous_codon',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::synonymous_codon',
  'label' => 'Synonymous coding',
  'rank' => '12',
  'display_term' => 'SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'non_synonymous_codon' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'missense',
  'feature_SO_term' => 'mRNA',
  'description' => 'In coding sequence and results in an amino acid change in the encoded peptide sequence',
  'SO_accession' => 'SO:0001583',
  'SO_term' => 'non_synonymous_codon',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::non_synonymous_codon',
  'label' => 'Non-synonymous coding',
  'rank' => '9',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'inframe_codon_gain' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'In coding sequence and results in an amino acid change in the encoded peptide sequence',
  'SO_accession' => 'SO:0001651',
  'SO_term' => 'inframe_codon_gain',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_codon_gain',
  'label' => 'Non-synonymous coding',
  'rank' => '8',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'inframe_codon_loss' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'In coding sequence and results in an amino acid change in the encoded peptide sequence',
  'SO_accession' => 'SO:0001652',
  'SO_term' => 'inframe_codon_loss',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_codon_loss',
  'label' => 'Non-synonymous coding',
  'rank' => '7',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'stop_gained' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'nonsense',
  'feature_SO_term' => 'mRNA',
  'description' => 'In coding sequence, resulting in the gain of a stop codon',
  'SO_accession' => 'SO:0001587',
  'SO_term' => 'stop_gained',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_gained',
  'label' => 'Stop gained',
  'rank' => '2',
  'display_term' => 'STOP_GAINED',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'stop_lost' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'In coding sequence, resulting in the loss of a stop codon',
  'SO_accession' => 'SO:0001578',
  'SO_term' => 'stop_lost',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_lost',
  'label' => 'Stop lost',
  'rank' => '3',
  'display_term' => 'STOP_LOST',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'stop_retained_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'In coding sequence, not resulting in an amino acid change (silent mutation)',
  'SO_accession' => 'SO:0001567',
  'SO_term' => 'stop_retained_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_retained',
  'label' => 'Synonymous coding',
  'rank' => '12',
  'display_term' => 'SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'initiator_codon_change' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'In coding sequence and results in an amino acid change in the encoded peptide sequence',
  'SO_accession' => 'SO:0001582',
  'SO_term' => 'initiator_codon_change',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::affects_start_codon',
  'label' => 'Non-synonymous coding',
  'rank' => '6',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'frameshift_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'frameshift',
  'feature_SO_term' => 'mRNA',
  'description' => 'In coding sequence, resulting in a frameshift',
  'SO_accession' => 'SO:0001589',
  'SO_term' => 'frameshift_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::frameshift',
  'label' => 'Frameshift coding',
  'rank' => '5',
  'display_term' => 'FRAMESHIFT_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'incomplete_terminal_codon_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'Located within the final, incomplete codon of a transcript whose end coordinate is unknown',
  'SO_accession' => 'SO:0001626',
  'SO_term' => 'incomplete_terminal_codon_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::partial_codon',
  'label' => 'Partial codon',
  'rank' => '11',
  'display_term' => 'PARTIAL_CODON',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'NMD_transcript_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'Located within a transcript predicted to undergo nonsense-mediated decay',
  'SO_accession' => 'SO:0001621',
  'SO_term' => 'NMD_transcript_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_nmd_transcript',
  'label' => 'NMD transcript',
  'rank' => '18',
  'display_term' => 'NMD_TRANSCRIPT',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'nc_transcript_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'ncRNA',
  'description' => 'Located within a gene that does not code for a protein',
  'SO_accession' => 'SO:0001619',
  'SO_term' => 'nc_transcript_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_non_coding_gene',
  'label' => 'Within non-coding gene',
  'rank' => '19',
  'display_term' => 'WITHIN_NON_CODING_GENE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'mature_miRNA_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'miRNA',
  'description' => 'Located within a microRNA',
  'SO_accession' => 'SO:0001620',
  'SO_term' => 'mature_miRNA_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_mature_miRNA',
  'label' => 'Within mature miRNA',
  'rank' => '14',
  'display_term' => 'WITHIN_MATURE_miRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'coding_sequence_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'In coding sequence with indeterminate effect',
  'SO_accession' => 'SO:0001580',
  'SO_term' => 'coding_sequence_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::coding_unknown',
  'label' => 'Coding unknown',
  'rank' => '13',
  'display_term' => 'CODING_UNKNOWN',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'regulatory_region_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'regulatory_region',
  'description' => 'In regulatory region annotated by Ensembl',
  'SO_accession' => 'SO:0001566',
  'SO_term' => 'regulatory_region_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_regulatory_feature',
  'label' => 'Regulatory region',
  'rank' => '25',
  'display_term' => 'REGULATORY_REGION',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature'
}
),
'TF_binding_site_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'TF_binding_site',
  'description' => 'In regulatory region annotated by Ensembl',
  'SO_accession' => 'SO:0001782',
  'SO_term' => 'TF_binding_site_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_motif_feature',
  'label' => 'Regulatory region',
  'rank' => '24',
  'display_term' => 'REGULATORY_REGION',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::MotifFeature'
}
),
'partial_overlap' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
  'SO_term' => 'partial_overlap',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::partial_overlap_feature',
  'feature_SO_term' => 'sequence_feature',
  'rank' => '999',
  'SO_accession' => 'SO:X000101',
  'feature_class' => 'Bio::EnsEMBL::Feature'
}
),
'entirely_within_feature' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
  'SO_term' => 'entirely_within_feature',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::complete_within_feature',
  'feature_SO_term' => 'sequence_feature',
  'rank' => '999',
  'SO_accession' => 'SO:X000102',
  'feature_class' => 'Bio::EnsEMBL::Feature'
}
),
'complete_overlap_feature' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
  'SO_term' => 'complete_overlap_feature',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::complete_overlap_feature',
  'feature_SO_term' => 'sequence_feature',
  'rank' => '999',
  'SO_accession' => 'SO:X000103',
  'feature_class' => 'Bio::EnsEMBL::Feature'
}
),
'deletion' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
  'SO_term' => 'deletion',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::deletion',
  'feature_SO_term' => 'sequence_feature',
  'rank' => '999',
  'SO_accession' => 'SO:0000159',
  'feature_class' => 'Bio::EnsEMBL::Feature'
}
),
'duplication' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
  'SO_term' => 'duplication',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::duplication',
  'feature_SO_term' => 'sequence_feature',
  'rank' => '999',
  'SO_accession' => 'SO:1000035',
  'feature_class' => 'Bio::EnsEMBL::Feature'
}
),
'tandem_duplication' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
  'SO_term' => 'tandem_duplication',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::tandem_duplication',
  'feature_SO_term' => 'sequence_feature',
  'rank' => '999',
  'SO_accession' => 'SO:1000173',
  'feature_class' => 'Bio::EnsEMBL::Feature'
}
),
'insertion' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::StructuralVariationFeature',
  'SO_term' => 'insertion',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::insertion',
  'feature_SO_term' => 'sequence_feature',
  'rank' => '999',
  'SO_accession' => 'SO:0000667',
  'feature_class' => 'Bio::EnsEMBL::Feature'
}
),
);

1;

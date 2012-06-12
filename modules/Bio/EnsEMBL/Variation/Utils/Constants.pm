package Bio::EnsEMBL::Variation::Utils::Constants;

#####################################################################
# NB: THIS FILE HAS BEEN AUTOMATICALLY GENERATED, EDIT WITH CAUTION #
#####################################################################

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(%OVERLAP_CONSEQUENCES %VARIATION_CLASSES $DEFAULT_OVERLAP_CONSEQUENCE SO_TERM_INDEL SO_TERM_CODING_SEQUENCE_VARIANT SO_TERM_STOP_RETAINED_VARIANT SO_TERM_COMPLEX_STRUCTURAL_ALTERATION SO_TERM_SPLICE_ACCEPTOR_VARIANT SO_TERM_TANDEM_REPEAT SO_TERM_INTERCHROMOSOMAL_BREAKPOINT SO_TERM_FEATURE_ELONGATION SO_TERM_TRANSCRIPT_ABLATION SO_TERM_TF_BINDING_SITE_VARIANT ATTRIB_TYPE_DISPLAY_TERM SO_TERM_INTRON_VARIANT SO_TERM_STOP_GAINED SO_TERM_SPLICE_DONOR_VARIANT ATTRIB_TYPE_SIFT_PREDICTION SO_TERM_FEATURE_TRUNCATION SO_TERM_INTRACHROMOSOMAL_BREAKPOINT SO_TERM_SPLICE_REGION_VARIANT SO_TERM_TFBS_ABLATION SO_TERM_REGULATORY_REGION_AMPLIFICATION SO_TERM_MISSENSE_VARIANT SO_TERM_SNV SO_TERM_STRUCTURAL_VARIANT SO_TERM_PROBE ATTRIB_TYPE_DBSNP_CLIN_SIG SO_TERM_SEQUENCE_ALTERATION SO_TERM_REGULATORY_REGION_VARIANT SO_TERM_SUBSTITUTION ATTRIB_TYPE_DGVA_CLIN_SIG SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT SO_TERM_TFBS_AMPLIFICATION SO_TERM_5_PRIME_UTR_VARIANT ATTRIB_TYPE_PROT_FUNC_ANALYSIS SO_TERM_TANDEM_DUPLICATION SO_TERM_MOBILE_ELEMENT_INSERTION SO_TERM_MATURE_MIRNA_VARIANT SO_TERM_NC_TRANSCRIPT_VARIANT SO_TERM_DOWNSTREAM_GENE_VARIANT SO_TERM_INFRAME_INSERTION SO_TERM_INSERTION SO_TERM_NMD_TRANSCRIPT_VARIANT ATTRIB_TYPE_SO_TERM SO_TERM_INTERGENIC_VARIANT SO_TERM_NON_CODING_EXON_VARIANT SO_TERM_SYNONYMOUS_VARIANT ATTRIB_TYPE_SHORT_NAME SO_TERM_TRANSLOCATION SO_TERM_COPY_NUMBER_VARIATION ATTRIB_TYPE_POLYPHEN_PREDICTION SO_TERM_REGULATORY_REGION_ABLATION SO_TERM_INVERSION SO_TERM_STOP_LOST ATTRIB_TYPE_RANK SO_TERM_INITIATOR_CODON_VARIANT SO_TERM_3_PRIME_UTR_VARIANT SO_TERM_UPSTREAM_GENE_VARIANT ATTRIB_TYPE_SO_ACCESSION SO_TERM_COPY_NUMBER_LOSS ATTRIB_TYPE_NCBI_TERM SO_TERM_DELETION ATTRIB_TYPE_FEATURE_SO_TERM SO_TERM_FRAMESHIFT_VARIANT SO_TERM_COPY_NUMBER_GAIN SO_TERM_INFRAME_DELETION SO_TERM_TRANSCRIPT_AMPLIFICATION);

our %EXPORT_TAGS = ( attrib_types => [qw(ATTRIB_TYPE_RANK ATTRIB_TYPE_NCBI_TERM ATTRIB_TYPE_SHORT_NAME ATTRIB_TYPE_FEATURE_SO_TERM ATTRIB_TYPE_DISPLAY_TERM ATTRIB_TYPE_SO_ACCESSION ATTRIB_TYPE_DBSNP_CLIN_SIG ATTRIB_TYPE_PROT_FUNC_ANALYSIS ATTRIB_TYPE_SO_TERM ATTRIB_TYPE_SIFT_PREDICTION ATTRIB_TYPE_POLYPHEN_PREDICTION ATTRIB_TYPE_DGVA_CLIN_SIG)], SO_consequence_terms => [qw(SO_TERM_NC_TRANSCRIPT_VARIANT SO_TERM_DOWNSTREAM_GENE_VARIANT SO_TERM_CODING_SEQUENCE_VARIANT SO_TERM_STOP_RETAINED_VARIANT SO_TERM_INFRAME_INSERTION SO_TERM_NMD_TRANSCRIPT_VARIANT SO_TERM_SPLICE_ACCEPTOR_VARIANT SO_TERM_INTERGENIC_VARIANT SO_TERM_FEATURE_ELONGATION SO_TERM_NON_CODING_EXON_VARIANT SO_TERM_TRANSCRIPT_ABLATION SO_TERM_SYNONYMOUS_VARIANT SO_TERM_TF_BINDING_SITE_VARIANT SO_TERM_INTRON_VARIANT SO_TERM_STOP_GAINED SO_TERM_SPLICE_DONOR_VARIANT SO_TERM_FEATURE_TRUNCATION SO_TERM_REGULATORY_REGION_ABLATION SO_TERM_TFBS_ABLATION SO_TERM_SPLICE_REGION_VARIANT SO_TERM_STOP_LOST SO_TERM_INITIATOR_CODON_VARIANT SO_TERM_3_PRIME_UTR_VARIANT SO_TERM_REGULATORY_REGION_AMPLIFICATION SO_TERM_UPSTREAM_GENE_VARIANT SO_TERM_MISSENSE_VARIANT SO_TERM_REGULATORY_REGION_VARIANT SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT SO_TERM_FRAMESHIFT_VARIANT SO_TERM_TFBS_AMPLIFICATION SO_TERM_5_PRIME_UTR_VARIANT SO_TERM_INFRAME_DELETION SO_TERM_TRANSCRIPT_AMPLIFICATION SO_TERM_MATURE_MIRNA_VARIANT)], SO_class_terms => [qw(SO_TERM_INDEL SO_TERM_SNV SO_TERM_STRUCTURAL_VARIANT SO_TERM_COMPLEX_STRUCTURAL_ALTERATION SO_TERM_INSERTION SO_TERM_PROBE SO_TERM_COPY_NUMBER_LOSS SO_TERM_SEQUENCE_ALTERATION SO_TERM_TANDEM_REPEAT SO_TERM_INTERCHROMOSOMAL_BREAKPOINT SO_TERM_SUBSTITUTION SO_TERM_DELETION SO_TERM_COPY_NUMBER_GAIN SO_TERM_TRANSLOCATION SO_TERM_COPY_NUMBER_VARIATION SO_TERM_TANDEM_DUPLICATION SO_TERM_INVERSION SO_TERM_INTRACHROMOSOMAL_BREAKPOINT SO_TERM_MOBILE_ELEMENT_INSERTION)],  );

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
use constant ATTRIB_TYPE_PROT_FUNC_ANALYSIS => 'prot_func_analysis';

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
use constant SO_TERM_INTERCHROMOSOMAL_BREAKPOINT => 'interchromosomal_breakpoint';
use constant SO_TERM_INTRACHROMOSOMAL_BREAKPOINT => 'intrachromosomal_breakpoint';
use constant SO_TERM_TRANSLOCATION => 'translocation';
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
use constant SO_TERM_INITIATOR_CODON_VARIANT => 'initiator_codon_variant';
use constant SO_TERM_FRAMESHIFT_VARIANT => 'frameshift_variant';
use constant SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT => 'incomplete_terminal_codon_variant';
use constant SO_TERM_NMD_TRANSCRIPT_VARIANT => 'NMD_transcript_variant';
use constant SO_TERM_NC_TRANSCRIPT_VARIANT => 'nc_transcript_variant';
use constant SO_TERM_NON_CODING_EXON_VARIANT => 'non_coding_exon_variant';
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
'interchromosomal_breakpoint' => {
  'somatic_display_term' => 'somatic_Interchromosomal breakpoint',
  'SO_accession' => 'SO:0001873',
  'display_term' => 'Interchromosomal breakpoint'
}
,
'intrachromosomal_breakpoint' => {
  'somatic_display_term' => 'somatic_Intrachromosomal breakpoint',
  'SO_accession' => 'SO:0001874',
  'display_term' => 'Intrachromosomal breakpoint'
}
,
'translocation' => {
  'somatic_display_term' => 'somatic_translocation',
  'SO_accession' => 'SO:0000199',
  'display_term' => 'translocation'
}
,
);

our $DEFAULT_OVERLAP_CONSEQUENCE = Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'tier' => '4',
  'SO_term' => 'intergenic_variant',
  'is_default' => 1,
  'label' => 'Intergenic variant',
  'description' => 'A sequence variant located in the intergenic region, between genes',
  'rank' => '38',
  'SO_accession' => 'SO:0001628',
  'display_term' => 'INTERGENIC'
}
);


our %OVERLAP_CONSEQUENCES = (
'intergenic_variant' => $DEFAULT_OVERLAP_CONSEQUENCE,
'upstream_gene_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'transcript',
  'description' => 'A sequence variant located 5\' of a gene',
  'SO_accession' => 'SO:0001631',
  'SO_term' => 'upstream_gene_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream',
  'label' => 'Upstream gene variant',
  'rank' => '24',
  'display_term' => 'UPSTREAM',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'downstream_gene_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'transcript',
  'description' => 'A sequence variant located 3\' of a gene',
  'SO_accession' => 'SO:0001632',
  'SO_term' => 'downstream_gene_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream',
  'label' => 'Downstream gene variant',
  'rank' => '25',
  'display_term' => 'DOWNSTREAM',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'splice_donor_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'splice-5',
  'feature_SO_term' => 'primary_transcript',
  'description' => 'A splice variant that changes the 2 base region at the 5\' end of an intron',
  'SO_accession' => 'SO:0001575',
  'tier' => '3',
  'SO_term' => 'splice_donor_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::donor_splice_site',
  'label' => 'Splice donor variant',
  'rank' => '3',
  'display_term' => 'ESSENTIAL_SPLICE_SITE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'splice_acceptor_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'splice-3',
  'feature_SO_term' => 'primary_transcript',
  'description' => 'A splice variant that changes the 2 base region at the 3\' end of an intron',
  'SO_accession' => 'SO:0001574',
  'tier' => '3',
  'SO_term' => 'splice_acceptor_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::acceptor_splice_site',
  'label' => 'Splice acceptor variant',
  'rank' => '3',
  'display_term' => 'ESSENTIAL_SPLICE_SITE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'splice_region_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'primary_transcript',
  'description' => 'A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron',
  'SO_accession' => 'SO:0001630',
  'SO_term' => 'splice_region_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::splice_region',
  'label' => 'Splice region variant',
  'rank' => '13',
  'display_term' => 'SPLICE_SITE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'intron_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'NCBI_term' => 'intron',
  'feature_SO_term' => 'primary_transcript',
  'description' => 'A transcript variant occurring within an intron',
  'SO_accession' => 'SO:0001627',
  'tier' => '3',
  'SO_term' => 'intron_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_intron',
  'label' => 'Intron variant',
  'rank' => '20',
  'display_term' => 'INTRONIC',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'5_prime_UTR_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'NCBI_term' => 'untranslated_5',
  'feature_SO_term' => 'mRNA',
  'description' => 'A UTR variant of the 5\' UTR',
  'SO_accession' => 'SO:0001623',
  'tier' => '3',
  'SO_term' => '5_prime_UTR_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_5_prime_utr',
  'label' => '5 prime UTR variant',
  'rank' => '18',
  'display_term' => '5PRIME_UTR',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'3_prime_UTR_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'NCBI_term' => 'untranslated_3',
  'feature_SO_term' => 'mRNA',
  'description' => 'A UTR variant of the 3\' UTR',
  'SO_accession' => 'SO:0001624',
  'tier' => '3',
  'SO_term' => '3_prime_UTR_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_3_prime_utr',
  'label' => '3 prime UTR variant',
  'rank' => '19',
  'display_term' => '3PRIME_UTR',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'synonymous_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'cds-synon',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant where there is no resulting change to the encoded amino acid',
  'SO_accession' => 'SO:0001819',
  'tier' => '3',
  'SO_term' => 'synonymous_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::synonymous_variant',
  'label' => 'Synonymous variant',
  'rank' => '15',
  'display_term' => 'SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'missense_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'missense',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant, where the change may be longer than 3 bases, and at least one base of a codon is changed resulting in a codon that encodes for a different amino acid',
  'SO_accession' => 'SO:0001583',
  'tier' => '3',
  'SO_term' => 'missense_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::missense_variant',
  'label' => 'Missense variant',
  'rank' => '12',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'inframe_insertion' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'An inframe non synonymous variant that inserts bases into in the coding sequence',
  'SO_accession' => 'SO:0001821',
  'SO_term' => 'inframe_insertion',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_insertion',
  'label' => 'Inframe insertion',
  'rank' => '10',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'inframe_deletion' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'An inframe non synonymous variant that deletes bases from the coding sequence',
  'SO_accession' => 'SO:0001822',
  'SO_term' => 'inframe_deletion',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_deletion',
  'label' => 'Inframe deletion',
  'rank' => '11',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'stop_gained' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'NCBI_term' => 'nonsense',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript',
  'SO_accession' => 'SO:0001587',
  'tier' => '3',
  'SO_term' => 'stop_gained',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_gained',
  'label' => 'Stop gained',
  'rank' => '4',
  'display_term' => 'STOP_GAINED',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'stop_lost' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript',
  'SO_accession' => 'SO:0001578',
  'SO_term' => 'stop_lost',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_lost',
  'label' => 'Stop lost',
  'rank' => '6',
  'display_term' => 'STOP_LOST',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'stop_retained_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant where at least one base in the terminator codon is changed, but the terminator remains',
  'SO_accession' => 'SO:0001567',
  'SO_term' => 'stop_retained_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_retained',
  'label' => 'Stop retained variant',
  'rank' => '15',
  'display_term' => 'SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'initiator_codon_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'A codon variant that changes at least one base of the first codon of a transcript',
  'SO_accession' => 'SO:0001582',
  'SO_term' => 'initiator_codon_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::affects_start_codon',
  'label' => 'Initiator codon variant',
  'rank' => '7',
  'display_term' => 'NON_SYNONYMOUS_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'frameshift_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'NCBI_term' => 'frameshift',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three',
  'SO_accession' => 'SO:0001589',
  'tier' => '3',
  'SO_term' => 'frameshift_variant',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::frameshift',
  'label' => 'Frameshift variant',
  'rank' => '5',
  'display_term' => 'FRAMESHIFT_CODING',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'incomplete_terminal_codon_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::VariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed',
  'SO_accession' => 'SO:0001626',
  'SO_term' => 'incomplete_terminal_codon_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::partial_codon',
  'label' => 'Incomplete terminal codon variant',
  'rank' => '14',
  'display_term' => 'PARTIAL_CODON',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'NMD_transcript_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'A variant in a transcript that is the target of NMD',
  'SO_accession' => 'SO:0001621',
  'SO_term' => 'NMD_transcript_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_nmd_transcript',
  'label' => 'NMD transcript variant',
  'rank' => '21',
  'display_term' => 'NMD_TRANSCRIPT',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'nc_transcript_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'ncRNA',
  'description' => 'A transcript variant of a non coding RNA',
  'SO_accession' => 'SO:0001619',
  'SO_term' => 'nc_transcript_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_non_coding_gene',
  'label' => 'NC transcript variant',
  'rank' => '23',
  'display_term' => 'WITHIN_NON_CODING_GENE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'non_coding_exon_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'ncRNA',
  'description' => 'A sequence variant that changes non-coding exon sequence',
  'SO_accession' => 'SO:0001792',
  'SO_term' => 'non_coding_exon_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::non_coding_exon_variant',
  'label' => 'Non coding exon variant',
  'rank' => '22',
  'display_term' => 'WITHIN_NON_CODING_GENE',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'mature_miRNA_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'miRNA',
  'description' => 'A transcript variant located with the sequence of the mature miRNA',
  'SO_accession' => 'SO:0001620',
  'SO_term' => 'mature_miRNA_variant',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_mature_miRNA',
  'label' => 'Mature miRNA variant',
  'rank' => '17',
  'display_term' => 'WITHIN_MATURE_miRNA',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'coding_sequence_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'A sequence variant that changes the coding sequence',
  'SO_accession' => 'SO:0001580',
  'SO_term' => 'coding_sequence_variant',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::coding_unknown',
  'label' => 'Coding sequence variant',
  'rank' => '16',
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
  'label' => 'Regulatory region variant',
  'rank' => '36',
  'display_term' => 'REGULATORY_REGION',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature'
}
),
'TF_binding_site_variant' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'TF_binding_site',
  'description' => 'In regulatory region annotated by Ensembl',
  'SO_accession' => 'SO:0001782',
  'SO_term' => 'TF_binding_site_variant',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_motif_feature',
  'label' => 'A sequence variant located within a transcription factor binding site',
  'rank' => '30',
  'display_term' => 'REGULATORY_REGION',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::MotifFeature'
}
),
'transcript_ablation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'A feature ablation whereby the deleted region includes a transcript feature',
  'SO_accession' => 'SO:0001893',
  'SO_term' => 'transcript_ablation',
  'tier' => '1',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
  'label' => 'Transcript ablation',
  'rank' => '1',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'transcript_amplification' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'mRNA',
  'description' => 'A feature amplification of a region containing a transcript',
  'SO_accession' => 'SO:0001889',
  'SO_term' => 'transcript_amplification',
  'tier' => '1',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
  'label' => 'Transcript amplification',
  'rank' => '8',
  'feature_class' => 'Bio::EnsEMBL::Transcript'
}
),
'TFBS_ablation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'TF_binding_site',
  'description' => 'A feature ablation whereby the deleted region includes a transcription factor binding site',
  'SO_accession' => 'SO:0001895',
  'SO_term' => 'TFBS_ablation',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
  'label' => 'TFBS ablation',
  'rank' => '26',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::MotifFeature'
}
),
'TFBS_amplification' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'TF_binding_site',
  'description' => 'A feature amplification of a region containing a transcription factor binding site',
  'SO_accession' => 'SO:0001892',
  'SO_term' => 'TFBS_amplification',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
  'label' => 'TFBS amplification',
  'rank' => '28',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::MotifFeature'
}
),
'regulatory_region_ablation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'TF_binding_site',
  'description' => 'A feature ablation whereby the deleted region includes a regulatory region',
  'SO_accession' => 'SO:0001894',
  'SO_term' => 'regulatory_region_ablation',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_ablation',
  'label' => 'Regulatory region ablation',
  'rank' => '31',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature'
}
),
'regulatory_region_amplification' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'feature_SO_term' => 'TF_binding_site',
  'description' => 'A feature amplification of a region containing a regulatory region',
  'SO_accession' => 'SO:0001891',
  'SO_term' => 'regulatory_region_amplification',
  'tier' => '2',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_amplification',
  'label' => 'Regulatory region amplification',
  'rank' => '33',
  'feature_class' => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature'
}
),
'feature_elongation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'description' => 'A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence',
  'SO_accession' => 'SO:0001907',
  'SO_term' => 'feature_elongation',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_elongation',
  'label' => 'Feature elongation',
  'rank' => '36',
  'feature_class' => 'Bio::EnsEMBL::Feature'
}
),
'feature_truncation' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
  'variant_feature_class' => 'Bio::EnsEMBL::Variation::BaseVariationFeature',
  'description' => 'A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence',
  'SO_accession' => 'SO:0001906',
  'SO_term' => 'feature_truncation',
  'tier' => '3',
  'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::feature_truncation',
  'label' => 'Feature truncation',
  'rank' => '37',
  'feature_class' => 'Bio::EnsEMBL::Feature'
}
),
);

1;

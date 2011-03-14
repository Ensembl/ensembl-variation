package Bio::EnsEMBL::Variation::Utils::Constants;

#####################################################################
# NB: THIS FILE HAS BEEN AUTOMATICALLY GENERATED, EDIT WITH CAUTION #
#####################################################################

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(@OVERLAP_CONSEQUENCES ATTRIB_TYPE_SO_ACCESSION ATTRIB_TYPE_SO_TERM ATTRIB_TYPE_DISPLAY_TERM ATTRIB_TYPE_NCBI_TERM ATTRIB_TYPE_FEATURE_SO_TERM ATTRIB_TYPE_RANK SO_TERM_INTERGENIC_VARIANT SO_TERM_5KB_UPSTREAM_VARIANT SO_TERM_5KB_DOWNSTREAM_VARIANT SO_TERM_2KB_UPSTREAM_VARIANT SO_TERM_500B_DOWNSTREAM_VARIANT SO_TERM_SPLICE_DONOR_VARIANT SO_TERM_SPLICE_ACCEPTOR_VARIANT SO_TERM_SPLICE_REGION_VARIANT SO_TERM_INTRON_VARIANT SO_TERM_5_PRIME_UTR_VARIANT SO_TERM_3_PRIME_UTR_VARIANT SO_TERM_COMPLEX_CHANGE_IN_TRANSCRIPT SO_TERM_SYNONYMOUS_CODON SO_TERM_NON_SYNONYMOUS_CODON SO_TERM_INFRAME_CODON_GAIN SO_TERM_INFRAME_CODON_LOSS SO_TERM_STOP_GAINED SO_TERM_STOP_LOST SO_TERM_STOP_RETAINED_VARIANT SO_TERM_INITIATOR_CODON_CHANGE SO_TERM_FRAMESHIFT_VARIANT SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT SO_TERM_NMD_TRANSCRIPT_VARIANT SO_TERM_NC_TRANSCRIPT_VARIANT SO_TERM_MATURE_MIRNA_VARIANT SO_TERM_CODING_SEQUENCE_VARIANT SO_TERM_REGULATORY_REGION_VARIANT SO_TERM_MIRNA_TARGET_SITE_VARIANT SO_TERM_BINDING_SITE_VARIANT SO_TERM_SNV SO_TERM_SUBSTITUTION SO_TERM_COPY_NUMBER_VARIATION SO_TERM_INSERTION SO_TERM_DELETION SO_TERM_INDEL SO_TERM_TANDEM_REPEAT SO_TERM_SEQUENCE_ALTERATION);

our %EXPORT_TAGS = ( attrib_types => [qw(ATTRIB_TYPE_SO_ACCESSION ATTRIB_TYPE_SO_TERM ATTRIB_TYPE_DISPLAY_TERM ATTRIB_TYPE_NCBI_TERM ATTRIB_TYPE_FEATURE_SO_TERM ATTRIB_TYPE_RANK)], SO_consequence_terms => [qw(SO_TERM_INTERGENIC_VARIANT SO_TERM_5KB_UPSTREAM_VARIANT SO_TERM_5KB_DOWNSTREAM_VARIANT SO_TERM_2KB_UPSTREAM_VARIANT SO_TERM_500B_DOWNSTREAM_VARIANT SO_TERM_SPLICE_DONOR_VARIANT SO_TERM_SPLICE_ACCEPTOR_VARIANT SO_TERM_SPLICE_REGION_VARIANT SO_TERM_INTRON_VARIANT SO_TERM_5_PRIME_UTR_VARIANT SO_TERM_3_PRIME_UTR_VARIANT SO_TERM_COMPLEX_CHANGE_IN_TRANSCRIPT SO_TERM_SYNONYMOUS_CODON SO_TERM_NON_SYNONYMOUS_CODON SO_TERM_INFRAME_CODON_GAIN SO_TERM_INFRAME_CODON_LOSS SO_TERM_STOP_GAINED SO_TERM_STOP_LOST SO_TERM_STOP_RETAINED_VARIANT SO_TERM_INITIATOR_CODON_CHANGE SO_TERM_FRAMESHIFT_VARIANT SO_TERM_INCOMPLETE_TERMINAL_CODON_VARIANT SO_TERM_NMD_TRANSCRIPT_VARIANT SO_TERM_NC_TRANSCRIPT_VARIANT SO_TERM_MATURE_MIRNA_VARIANT SO_TERM_CODING_SEQUENCE_VARIANT SO_TERM_REGULATORY_REGION_VARIANT SO_TERM_MIRNA_TARGET_SITE_VARIANT SO_TERM_BINDING_SITE_VARIANT)], SO_class_terms => [qw(SO_TERM_SNV SO_TERM_SUBSTITUTION SO_TERM_COPY_NUMBER_VARIATION SO_TERM_INSERTION SO_TERM_DELETION SO_TERM_INDEL SO_TERM_TANDEM_REPEAT SO_TERM_SEQUENCE_ALTERATION)],  );

use Bio::EnsEMBL::Variation::OverlapConsequence;

use constant ATTRIB_TYPE_SO_ACCESSION => 'SO_accession';
use constant ATTRIB_TYPE_SO_TERM => 'SO_term';
use constant ATTRIB_TYPE_DISPLAY_TERM => 'display_term';
use constant ATTRIB_TYPE_NCBI_TERM => 'NCBI_term';
use constant ATTRIB_TYPE_FEATURE_SO_TERM => 'feature_SO_term';
use constant ATTRIB_TYPE_RANK => 'rank';

use constant SO_TERM_SNV => 'SNV';
use constant SO_TERM_SUBSTITUTION => 'substitution';
use constant SO_TERM_COPY_NUMBER_VARIATION => 'copy_number_variation';
use constant SO_TERM_INSERTION => 'insertion';
use constant SO_TERM_DELETION => 'deletion';
use constant SO_TERM_INDEL => 'indel';
use constant SO_TERM_TANDEM_REPEAT => 'tandem_repeat';
use constant SO_TERM_SEQUENCE_ALTERATION => 'sequence_alteration';

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
use constant SO_TERM_MIRNA_TARGET_SITE_VARIANT => 'miRNA_target_site_variant';
use constant SO_TERM_BINDING_SITE_VARIANT => 'binding_site_variant';

our @OVERLAP_CONSEQUENCES = (
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'intergenic_variant',
          'rank' => '100',
          'SO_accession' => 'SO:0001628',
          'display_term' => 'INTERGENIC'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => '5KB_upstream_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream_5KB',
          'feature_SO_term' => 'transcript',
          'rank' => '20',
          'SO_accession' => 'SO:0001635',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'UPSTREAM'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => '5KB_downstream_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream_5KB',
          'feature_SO_term' => 'transcript',
          'rank' => '21',
          'SO_accession' => 'SO:0001633',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'DOWNSTREAM'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => '2KB_upstream_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream_2KB',
          'NCBI_term' => 'near-gene-5',
          'feature_SO_term' => 'transcript',
          'rank' => '18',
          'SO_accession' => 'SO:0001636',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'UPSTREAM'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => '500B_downstream_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream_500B',
          'NCBI_term' => 'near-gene-3',
          'feature_SO_term' => 'transcript',
          'rank' => '19',
          'SO_accession' => 'SO:0001634',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'DOWNSTREAM'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'splice_donor_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::donor_splice_site',
          'NCBI_term' => 'splice-5',
          'feature_SO_term' => 'primary_transcript',
          'rank' => '1',
          'SO_accession' => 'SO:0001575',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'ESSENTIAL_SPLICE_SITE'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'splice_acceptor_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::acceptor_splice_site',
          'NCBI_term' => 'splice-3',
          'feature_SO_term' => 'primary_transcript',
          'rank' => '1',
          'SO_accession' => 'SO:0001574',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'ESSENTIAL_SPLICE_SITE'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'splice_region_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::splice_region',
          'feature_SO_term' => 'primary_transcript',
          'rank' => '8',
          'SO_accession' => 'SO:0001630',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'SPLICE_SITE'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'intron_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_intron',
          'NCBI_term' => 'intron',
          'feature_SO_term' => 'primary_transcript',
          'rank' => '15',
          'SO_accession' => 'SO:0001627',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'INTRONIC'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => '5_prime_UTR_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_5_prime_utr',
          'NCBI_term' => 'untranslated_5',
          'feature_SO_term' => 'mRNA',
          'rank' => '13',
          'SO_accession' => 'SO:0001623',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => '5PRIME_UTR'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => '3_prime_UTR_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_3_prime_utr',
          'NCBI_term' => 'untranslated_3',
          'feature_SO_term' => 'mRNA',
          'rank' => '14',
          'SO_accession' => 'SO:0001624',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => '3PRIME_UTR'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'complex_change_in_transcript',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::complex_indel',
          'feature_SO_term' => 'primary_transcript',
          'rank' => '5',
          'SO_accession' => 'SO:0001577',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'COMPLEX_INDEL'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'synonymous_codon',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::synonymous_codon',
          'NCBI_term' => 'cds-synon',
          'feature_SO_term' => 'mRNA',
          'rank' => '10',
          'SO_accession' => 'SO:0001588',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'SYNONYMOUS_CODING'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'non_synonymous_codon',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::non_synonymous_codon',
          'NCBI_term' => 'missense',
          'feature_SO_term' => 'mRNA',
          'rank' => '7',
          'SO_accession' => 'SO:0001583',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'NON_SYNONYMOUS_CODING'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'inframe_codon_gain',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_codon_gain',
          'feature_SO_term' => 'mRNA',
          'rank' => '6',
          'SO_accession' => 'SO:0001651',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'NON_SYNONYMOUS_CODING'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'inframe_codon_loss',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::inframe_codon_loss',
          'feature_SO_term' => 'mRNA',
          'rank' => '6',
          'SO_accession' => 'SO:0001652',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'NON_SYNONYMOUS_CODING'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'stop_gained',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_gained',
          'NCBI_term' => 'nonsense',
          'feature_SO_term' => 'mRNA',
          'rank' => '3',
          'SO_accession' => 'SO:0001587',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'STOP_GAINED'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'stop_lost',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_lost',
          'feature_SO_term' => 'mRNA',
          'rank' => '4',
          'SO_accession' => 'SO:0001578',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'STOP_LOST'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'stop_retained_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_retained',
          'feature_SO_term' => 'mRNA',
          'rank' => '10',
          'SO_accession' => 'SO:0001567',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'SYNONYMOUS_CODING'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'initiator_codon_change',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::affects_start_codon',
          'feature_SO_term' => 'mRNA',
          'rank' => '7',
          'SO_accession' => 'SO:0001582',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'NON_SYNONYMOUS_CODING'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'frameshift_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::frameshift',
          'NCBI_term' => 'frameshift',
          'feature_SO_term' => 'mRNA',
          'rank' => '6',
          'SO_accession' => 'SO:0001589',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'FRAMESHIFT_CODING'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'incomplete_terminal_codon_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::partial_codon',
          'feature_SO_term' => 'mRNA',
          'rank' => '9',
          'SO_accession' => 'SO:0001626',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'PARTIAL_CODON'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'NMD_transcript_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_nmd_transcript',
          'feature_SO_term' => 'mRNA',
          'rank' => '16',
          'SO_accession' => 'SO:0001621',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'NMD_TRANSCRIPT'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'nc_transcript_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_non_coding_gene',
          'feature_SO_term' => 'ncRNA',
          'rank' => '17',
          'SO_accession' => 'SO:0001619',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'WITHIN_NON_CODING_GENE'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'mature_miRNA_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_mature_miRNA',
          'feature_SO_term' => 'miRNA',
          'rank' => '12',
          'SO_accession' => 'SO:0001620',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'WITHIN_MATURE_miRNA'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'coding_sequence_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::coding_unknown',
          'feature_SO_term' => 'mRNA',
          'rank' => '11',
          'SO_accession' => 'SO:0001580',
          'feature_class' => 'Bio::EnsEMBL::Transcript',
          'display_term' => 'CODING_UNKNOWN'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'regulatory_region_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_regulatory_feature',
          'feature_SO_term' => 'regulatory_region',
          'rank' => '50',
          'SO_accession' => 'SO:0001566',
          'feature_class' => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
          'display_term' => 'REGULATORY_REGION'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'miRNA_target_site_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_miRNA_target_site',
          'feature_SO_term' => 'binding_site',
          'rank' => '13',
          'SO_accession' => 'SO:X000004',
          'feature_class' => 'Bio::EnsEMBL::Funcgen::ExternalFeature',
          'display_term' => 'REGULATORY_REGION'
        }
),
Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
          'SO_term' => 'binding_site_variant',
          'predicate' => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::within_motif_feature',
          'feature_SO_term' => 'binding_site',
          'rank' => '49',
          'SO_accession' => 'SO:X000003',
          'feature_class' => 'Bio::EnsEMBL::Funcgen::MotifFeature',
          'display_term' => 'REGULATORY_REGION'
        }
),
);

1;

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

# EnsEMBL module for Bio::EnsEMBL::Variation::Utils::Sequence
#
#

=head1 NAME

Bio::EnsEMBL::Variation::Utils::VEP - Methods used by the Variant Effect Predictor

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::Utils::VEP qw(configure);

  my $config = configure();

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Utils::VEP;

# module list
use Getopt::Long;
use FileHandle;
use File::Path qw(mkpath);
use Storable qw(nstore_fd fd_retrieve freeze thaw);
use Scalar::Util qw(weaken looks_like_number);
use Digest::MD5 qw(md5_hex);
use IO::Socket;
use IO::Select;
use Exporter;

my ($CAN_USE_TABIX_PM, $CAN_USE_TABIX_CL);

BEGIN {
  if (eval { require Tabix; 1 }) {
    $CAN_USE_TABIX_PM = 1;
  }
  elsif (`which tabix 2>&1` =~ /tabix$/) {
    $CAN_USE_TABIX_CL = 1;
  }

  # use Sereal
  eval q{ use Sereal; 1; };
}

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code SO_variation_class);
use Bio::EnsEMBL::Variation::Utils::EnsEMBL2GFF3;
use Bio::EnsEMBL::Variation::StructuralVariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::TranscriptStructuralVariation;
use Bio::EnsEMBL::Variation::Source;

# we need to manually include all these modules for caching to work
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
use Bio::EnsEMBL::Funcgen::MotifFeature;
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::TranslationAdaptor;
use Bio::EnsEMBL::DBSQL::TranscriptAdaptor;
use Bio::EnsEMBL::DBSQL::MetaContainer;
use Bio::EnsEMBL::DBSQL::CoordSystemAdaptor;

use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Exporter);

@EXPORT_OK = qw(
    &detect_format
    &parse_line
    &vf_to_consequences
    &validate_vf
    &read_cache_info
    &get_version_data
    &load_dumped_variation_cache
    &load_dumped_transcript_cache
    &prefetch_transcript_data
    &get_all_consequences
    &get_slice
    &build_slice_cache
    &build_full_cache
    &regions_from_hash
    &get_time
    &debug
    &convert_to_vcf
    &progress
    &end_progress
    @REG_FEAT_TYPES
    @OUTPUT_COLS
    @VCF_COLS
    @EXTRA_HEADERS
    %COL_DESCS
    @VEP_WEB_CONFIG
    %FILTER_SHORTCUTS
    @PICK_ORDER
);

our @OUTPUT_COLS = qw(
    Uploaded_variation
    Location
    Allele
    Gene
    Feature
    Feature_type
    Consequence
    cDNA_position
    CDS_position
    Protein_position
    Amino_acids
    Codons
    Existing_variation
    Extra
);

our @VCF_COLS = (
  'Allele',
  'Consequence',
  'IMPACT',
  'SYMBOL',
  'Gene',
  'Feature_type',
  'Feature',
  'BIOTYPE',
  'EXON',
  'INTRON',
  'HGVSc',
  'HGVSp',
  'cDNA_position',
  'CDS_position',
  'Protein_position',
);

# define headers that would normally go in the extra field
# keyed on the config parameter used to turn it on
our @EXTRA_HEADERS = (

  # general
  { flag => 'individual',      cols => ['IND','ZYG'] },
  { flag => 'allele_number',   cols => ['ALLELE_NUM'] },
  { flag => 'user',            cols => ['IMPACT','DISTANCE','STRAND','FLAGS'] },
  { flag => 'flag_pick',       cols => ['PICK'] },
  { flag => 'flag_pick_allele',cols => ['PICK'] },
  { flag => 'variant_class',   cols => ['VARIANT_CLASS']},
  { flag => 'minimal',         cols => ['MINIMISED']},

  # gene-related
  { flag => 'symbol',          cols => ['SYMBOL','SYMBOL_SOURCE','HGNC_ID'] },
  { flag => 'biotype',         cols => ['BIOTYPE'] },
  { flag => 'canonical',       cols => ['CANONICAL'] },
  { flag => 'tsl',             cols => ['TSL']},
  { flag => 'appris',          cols => ['APPRIS']},
  { flag => 'ccds',            cols => ['CCDS'] },
  { flag => 'protein',         cols => ['ENSP'] },
  { flag => 'uniprot',         cols => ['SWISSPROT', 'TREMBL', 'UNIPARC'] },
  { flag => 'xref_refseq',     cols => ['RefSeq'] },
  { flag => 'refseq',          cols => ['REFSEQ_MATCH'] },
  { flag => 'merged',          cols => ['REFSEQ_MATCH'] },
  { flag => 'gene_phenotype',  cols => ['GENE_PHENO'] },

  # non-synonymous predictions
  { flag => 'sift',            cols => ['SIFT'] },
  { flag => 'polyphen',        cols => ['PolyPhen'] },

  # transcript/protein stuff
  { flag => 'numbers',         cols => ['EXON','INTRON'] },
  { flag => 'domains',         cols => ['DOMAINS'] },
  { flag => 'hgvs',            cols => ['HGVSc','HGVSp','HGVS_OFFSET'] },

  # frequency stuff
  { flag => 'gmaf',            cols => ['GMAF'] },
  { flag => 'maf_1kg',         cols => ['AFR_MAF','AMR_MAF','EAS_MAF','EUR_MAF','SAS_MAF'] },
  { flag => 'maf_esp',         cols => ['AA_MAF','EA_MAF'] },
  { flag => 'maf_exac',        cols => ['ExAC_MAF','ExAC_Adj_MAF','ExAC_AFR_MAF','ExAC_AMR_MAF','ExAC_EAS_MAF','ExAC_FIN_MAF','ExAC_NFE_MAF','ExAC_OTH_MAF','ExAC_SAS_MAF'] },
  { flag => 'check_frequency', cols => ['FREQS'] },

  # misc variation stuff
  { flag => 'check_existing',  cols => ['CLIN_SIG','SOMATIC','PHENO'] },
  { flag => 'pubmed',          cols => ['PUBMED'] },
  { flag => 'check_svs',       cols => ['SV'] },

  # regulatory
  { flag => 'regulatory',      cols => ['MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE'] },
  { flag => 'cell_type',       cols => ['CELL_TYPE'] },
);

our %COL_DESCS = (
    'Uploaded_variation' => 'Identifier of uploaded variant',
    'ID'                 => 'Identifier of uploaded variant',
    'Location'           => 'Location of variant in standard coordinate format (chr:start or chr:start-end)',
    'Allele'             => 'The variant allele used to calculate the consequence',
    'Gene'               => 'Stable ID of affected gene',
    'Feature'            => 'Stable ID of feature',
    'Feature_type'       => 'Type of feature - Transcript, RegulatoryFeature or MotifFeature',
    'Consequence'        => 'Consequence type',
    'cDNA_position'      => 'Relative position of base pair in cDNA sequence',
    'CDS_position'       => 'Relative position of base pair in coding sequence',
    'Protein_position'   => 'Relative position of amino acid in protein',
    'Amino_acids'        => 'Reference and variant amino acids',
    'Codons'             => 'Reference and variant codon sequence',
    'Existing_variation' => 'Identifier(s) of co-located known variants',
    'IMPACT'             => 'Subjective impact classification of consequence type',
    'CANONICAL'          => 'Indicates if transcript is canonical for this gene',
    'TSL'                => 'Transcript support level',
    'APPRIS'             => 'Annotates alternatively spliced transcripts as primary or alternate based on a range of computational methods',
    'CCDS'               => 'Indicates if transcript is a CCDS transcript',
    'SYMBOL'             => 'Gene symbol (e.g. HGNC)',
    'SYMBOL_SOURCE'      => 'Source of gene symbol',
    'SOURCE'             => 'Source of transcript in merged gene set',
    'HGNC_ID'            => 'Stable identifer of HGNC gene symbol',
    'ENSP'               => 'Protein identifer',
    'FLAGS'              => 'Transcript quality flags',
    'SWISSPROT'          => 'Best match UniProtKB/Swiss-Prot accession',
    'TREMBL'             => 'Best match UniProtKB/TrEMBL accession',
    'UNIPARC'            => 'Best match UniParc accession',
    'HGVSc'              => 'HGVS coding sequence name',
    'HGVSp'              => 'HGVS protein sequence name',
    'SIFT'               => 'SIFT prediction and/or score',
    'PolyPhen'           => 'PolyPhen prediction and/or score',
    'EXON'               => 'Exon number(s) / total',
    'INTRON'             => 'Intron number(s) / total',
    'DOMAINS'            => 'The source and identifer of any overlapping protein domains',
    'MOTIF_NAME'         => 'The source and identifier of a transcription factor binding profile (TFBP) aligned at this position',
    'MOTIF_POS'          => 'The relative position of the variation in the aligned TFBP',
    'HIGH_INF_POS'       => 'A flag indicating if the variant falls in a high information position of the TFBP',
    'MOTIF_SCORE_CHANGE' => 'The difference in motif score of the reference and variant sequences for the TFBP',
    'CELL_TYPE'          => 'List of cell types and classifications for regulatory feature',
    'IND'                => 'Individual name',
    'ZYG'                => 'Zygosity of individual genotype at this locus',
    'SV'                 => 'IDs of overlapping structural variants',
    'FREQS'              => 'Frequencies of overlapping variants used in filtering',
    'GMAF'               => 'Minor allele and frequency of existing variant in 1000 Genomes combined population',
    'AFR_MAF'            => 'Frequency of existing variant in 1000 Genomes combined African population',
    'AMR_MAF'            => 'Frequency of existing variant in 1000 Genomes combined American population',
    'ASN_MAF'            => 'Frequency of existing variant in 1000 Genomes combined Asian population',
    'EAS_MAF'            => 'Frequency of existing variant in 1000 Genomes combined East Asian population',
    'EUR_MAF'            => 'Frequency of existing variant in 1000 Genomes combined European population',
    'SAS_MAF'            => 'Frequency of existing variant in 1000 Genomes combined South Asian population',
    'AA_MAF'             => 'Frequency of existing variant in NHLBI-ESP African American population',
    'EA_MAF'             => 'Frequency of existing variant in NHLBI-ESP European American population',
    'ExAC_MAF',          => 'Frequency of existing variant in ExAC combined population',
    'ExAC_Adj_MAF',      => 'Adjusted frequency of existing variant in ExAC combined population',
    'ExAC_AFR_MAF',      => 'Frequency of existing variant in ExAC African/American population',
    'ExAC_AMR_MAF',      => 'Frequency of existing variant in ExAC American population',
    'ExAC_EAS_MAF',      => 'Frequency of existing variant in ExAC East Asian population',
    'ExAC_FIN_MAF',      => 'Frequency of existing variant in ExAC Finnish population',
    'ExAC_NFE_MAF',      => 'Frequency of existing variant in ExAC Non-Finnish European population',
    'ExAC_OTH_MAF',      => 'Frequency of existing variant in ExAC combined other combined populations',
    'ExAC_SAS_MAF',      => 'Frequency of existing variant in ExAC South Asian population',
    'DISTANCE'           => 'Shortest distance from variant to transcript',
    'CLIN_SIG'           => 'ClinVar clinical significance of the dbSNP variant',
    'BIOTYPE'            => 'Biotype of transcript or regulatory feature',
    'PUBMED'             => 'Pubmed ID(s) of publications that cite existing variant',
    'ALLELE_NUM'         => 'Allele number from input; 0 is reference, 1 is first alternate etc',
    'STRAND'             => 'Strand of the feature (1/-1)',
    'PICK'               => 'Indicates if this consequence has been picked as the most severe',
    'SOMATIC'            => 'Somatic status of existing variant',
    'REFSEQ_MATCH'       => 'RefSeq transcript match status',
    'VARIANT_CLASS'      => 'SO variant class',
    'PHENO'              => 'Indicates if existing variant(s) is associated with a phenotype, disease or trait; multiple values correspond to multiple variants',
    'GENE_PHENO'         => 'Indicates if gene is associated with a phenotype, disease or trait',
    'MINIMISED'          => 'Alleles in this variant have been converted to minimal representation before consequence calculation',
    'HGVS_OFFSET'        => 'Indicates by how many bases the HGVS notations for this variant have been shifted',
);

our @REG_FEAT_TYPES = qw(
    RegulatoryFeature
    MotifFeature
);

our @VEP_WEB_CONFIG = qw(
    format
    check_existing
    coding_only
    core_type
    symbol
    protein
    hgvs
    terms
    check_frequency
    freq_filter
    freq_gt_lt
    freq_freq
    freq_pop
    filter_common
    sift
    polyphen
    regulatory
);

our @VAR_CACHE_COLS = qw(
    variation_name
    failed
    somatic
    start
    end
    allele_string
    strand
    minor_allele
    minor_allele_freq
    clin_sig
    phenotype_or_disease
);

our @PICK_ORDER = qw(canonical appris tsl biotype ccds rank length ensembl refseq);

# parses a line of input, returns VF object(s)
sub parse_line {
    my $config = shift;
    my $line   = shift;

    # find out file format - will only do this on first line
    if(!defined($config->{format}) || (defined($config->{format}) && $config->{format} eq 'guess')) {
        $config->{format} = &detect_format($line);
        debug("Detected format of input file as ", $config->{format}) unless defined($config->{quiet});

        # HGVS and ID formats need DB
        die("ERROR: Can't use ".uc($config->{format})." format in offline mode") if $config->{format} =~ /id|hgvs/ && defined($config->{offline});
    }

    # check that format is vcf when using --individual
    die("ERROR: --individual only compatible with VCF input files\n") if defined($config->{individual}) && $config->{format} ne 'vcf';

    my $parse_method = 'parse_'.$config->{format};
    $parse_method =~ s/vep_//;
    my $method_ref   = \&$parse_method;

    my $vfs = &$method_ref($config, $line);

    $vfs = minimise_alleles($config, $vfs) if defined($config->{minimal});

    $vfs = add_lrg_mappings($config, $vfs) if defined($config->{lrg});

    $_->{_line} = $line for @$vfs;

    return $vfs;
}

# sub-routine to detect format of input
sub detect_format {
    my $line = shift;
    my @data = split /\s+/, $line;

    # HGVS: ENST00000285667.3:c.1047_1048insC
    if (
        scalar @data == 1 &&
        $data[0] =~ /^([^\:]+)\:.*?([cgmrp]?)\.?([\*\-0-9]+.*)$/i
    ) {
        return 'hgvs';
    }

    # variant identifier: rs123456
    elsif (
        scalar @data == 1
    ) {
        return 'id';
    }

    # VCF: 20  14370  rs6054257  G  A  29  0  NS=58;DP=258;AF=0.786;DB;H2  GT:GQ:DP:HQ
    elsif (
        $data[0] =~ /(chr)?\w+/ &&
        $data[1] =~ /^\d+$/ &&
        $data[3] =~ /^[ACGTN\-\.]+$/i &&
        $data[4] && $data[4] =~ /^([\.ACGTN\-\*]+\,?)+$|^(\<[\w]+\>)$/i
    ) {
        return 'vcf';
    }

    # pileup: chr1  60  T  A
    elsif (
        $data[0] =~ /(chr)?\w+/ &&
        $data[1] =~ /^\d+$/ &&
        $data[2] =~ /^[\*ACGTN-]+$/i &&
        $data[3] =~ /^[\*ACGTNRYSWKM\+\/-]+$/i
    ) {
        return 'pileup';
    }

    # ensembl: 20  14370  14370  A/G  +
    elsif (
        $data[0] =~ /\w+/ &&
        $data[1] =~ /^\d+$/ &&
        $data[2] =~ /^\d+$/ &&
        $data[3] =~ /(ins|dup|del)|([ACGTN-]+\/[ACGTN-]+)/i
    ) {
        return 'ensembl';
    }

    else {
        die("ERROR: Could not detect input file format\n");
    }
}

# parse a line of Ensembl format input into a variation feature object
sub parse_ensembl {
    my $config = shift;
    my $line = shift;

    my ($chr, $start, $end, $allele_string, $strand, $var_name) = split /\s+/, $line;

    # simple validity check
    unless($chr && $start && $end && $allele_string) {
      warning_msg($config, "Invalid input formatting on line ".$config->{line_number});
      return [];
    }

    $strand = 1 if !defined($strand);

    my $vf;

    # sv?
    if($allele_string !~ /\//) {
        my $so_term;

        # convert to SO term
        my %terms = (
            INS  => 'insertion',
            DEL  => 'deletion',
            TDUP => 'tandem_duplication',
            DUP  => 'duplication'
        );

        $so_term = defined $terms{$allele_string} ? $terms{$allele_string} : $allele_string;

        $vf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new_fast({
            start          => $start,
            end            => $end,
            strand         => $strand =~ /\-/ ? -1 : 1,
            adaptor        => $config->{svfa},
            variation_name => $var_name,
            chr            => $chr,
            class_SO_term  => $so_term,
        });
    }

    # normal vf
    else {
        $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
            start          => $start,
            end            => $end,
            allele_string  => $allele_string,
            strand         => $strand =~ /\-/ ? -1 : 1,
            map_weight     => 1,
            adaptor        => $config->{vfa},
            variation_name => $var_name,
            chr            => $chr,
        });
    }

    return [$vf];
}

# parse a line of VCF input into a variation feature object
sub parse_vcf {
    my $config = shift;
    my $line = shift;

    my @data = split /\s+/, $line;

    # get relevant data
    my ($chr, $start, $end, $ref, $alt) = ($data[0], $data[1], $data[1], $data[3], $data[4]);

    # simple validity check
    unless($chr && $start && $end && $ref && $alt) {
      warning_msg($config, "Invalid input formatting on line ".$config->{line_number});
      return [];
    }

    # non-variant
    my $non_variant = 0;

    if($alt eq '.') {
        if(defined($config->{allow_non_variant})) {
            $non_variant = 1;
        }
        else {
            return [];
        }
    }

    # some VCF files have a GRCh37 pos defined in GP flag in INFO column
    # if user has requested, we can use that as the position instead
    if(defined $config->{gp}) {
        $chr = undef;
        $start = undef;

        foreach my $pair(split ';', $data[7]) {
            my ($key, $value) = split '=', $pair;
            if($key eq 'GP') {
                ($chr, $start) = split ':', $value;
                $end = $start;
            }
        }

        unless(defined($chr) and defined($start)) {
            warning_msg($config, "No GP flag found in INFO column on line ".$config->{line_number});
            return [];
        }
    }

    # adjust end coord
    $end += (length($ref) - 1);

    # structural variation
    if(
      $ref.$alt !~ /^[ACGT\,]+$/ &&
      (
        (defined($data[7]) && $data[7] =~ /SVTYPE/) ||
        $alt =~ /[\<|\[]^\*[\]\>]/
      )
    ) {

        # parse INFO field
        my %info = ();

        foreach my $bit(split ';', ($data[7] || '')) {
            my ($key, $value) = split '=', $bit;
            $info{$key} = $value;
        }

        # like indels, SVs have the base before included for reference
        $start++;

        # work out the end coord
        if(defined($info{END})) {
            $end = $info{END};
        }
        elsif(defined($info{SVLEN})) {
            $end = $start + abs($info{SVLEN}) - 1;
        }

        # check for imprecise breakpoints
        my ($min_start, $max_start, $min_end, $max_end);

        if(defined($info{CIPOS})) {
            my ($low, $high) = split ',', $info{CIPOS};
            $min_start = $start + $low;
            $max_start = $start + $high;
        }

        if(defined($info{CIEND})) {
            my ($low, $high) = split ',', $info{CIEND};
            $min_end = $end + $low;
            $max_end = $end + $high;
        }

        # get type
        my $type;

        if($alt =~ /\<|\[|\]|\>/) {
            $type = $alt;
            $type =~ s/\<|\>//g;
            $type =~ s/\:.+//g;

            if($start >= $end && $type =~ /del/i) {
                warning_msg($config, "WARNING: VCF line on line ".$config->{line_number}." looks incomplete, skipping:\n$line\n");
                return [];
            }

        }
        else {
            $type = $info{SVTYPE};
        }

        my $so_term;

        if(defined($type)) {
            # convert to SO term
            my %terms = (
                INS  => 'insertion',
                DEL  => 'deletion',
                TDUP => 'tandem_duplication',
                DUP  => 'duplication'
            );

            $so_term = defined $terms{$type} ? $terms{$type} : $type;
        }

        my $svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new_fast({
            start          => $start,
            inner_start    => $max_start,
            outer_start    => $min_start,
            end            => $end,
            inner_end      => $min_end,
            outer_end      => $max_end,
            strand         => 1,
            adaptor        => $config->{svfa},
            variation_name => $data[2] eq '.' ? undef : $data[2],
            chr            => $chr,
            class_SO_term  => $so_term,
        });

        return [$svf];
    }

    # normal variation
    else {

        # find out if any of the alt alleles make this an insertion or a deletion
        my ($is_indel, $is_sub, $ins_count, $total_count);
        foreach my $alt_allele(split ',', $alt) {
            $is_indel = 1 if $alt_allele =~ /^[DI]/;
            $is_indel = 1 if length($alt_allele) != length($ref);
            $is_sub = 1 if length($alt_allele) == length($ref);
            $ins_count++ if length($alt_allele) > length($ref);
            $total_count++;
        }

        # multiple alt alleles?
        if($alt =~ /\,/) {
            if($is_indel) {
                my @alts;

                # find out if all the alts start with the same base
                # ignore "*"-types
                my %first_bases = map {substr($_, 0, 1) => 1} grep {!/\*/} ($ref, split(',', $alt));

                if(scalar keys %first_bases == 1) {
                    $ref = substr($ref, 1) || '-';
                    $start++;

                    foreach my $alt_allele(split ',', $alt) {
                        $alt_allele = substr($alt_allele, 1) unless $alt_allele =~ /\*/;
                        $alt_allele = '-' if $alt_allele eq '';
                        push @alts, $alt_allele;
                    }
                }
                else {
                    push @alts, split(',', $alt);
                }

                $alt = join "/", @alts;
            }

            else {
                # for substitutions we just need to replace ',' with '/' in $alt
                $alt =~ s/\,/\//g;
            }
        }

        elsif($is_indel) {

            # insertion or deletion (VCF 4+)
            if(substr($ref, 0, 1) eq substr($alt, 0, 1)) {

                # chop off first base
                $ref = substr($ref, 1) || '-';
                $alt = substr($alt, 1) || '-';

                $start++;
            }
        }

        my $original_alt = $alt;

        # create VF object
        my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
            start          => $start,
            end            => $end,
            allele_string  => $non_variant ? $ref : $ref.'/'.$alt,
            strand         => 1,
            map_weight     => 1,
            adaptor        => $config->{vfa},
            variation_name => $data[2] eq '.' ? undef : $data[2],
            chr            => $chr,
        });

        # flag as non-variant
        $vf->{non_variant} = 1 if $non_variant;

        # individuals?
        if(defined($config->{individual})) {
            my @alleles = split '\/', $ref.'/'.$original_alt;

            my @return;

            foreach my $ind(keys %{$config->{ind_cols}}) {

                # get alleles present in this individual
                my @bits;
                my $gt = (split ':', $data[$config->{ind_cols}->{$ind}])[0];

                my $phased = ($gt =~ /\|/ ? 1 : 0);

                foreach my $bit(split /\||\/|\\/, $gt) {
                    push @bits, $alleles[$bit] unless $bit eq '.';
                }

                # shallow copy VF
                my $vf_copy = { %$vf };
                bless $vf_copy, ref($vf);

                # get non-refs, remembering to exclude "*"-types
                my %non_ref = map {$_ => 1} grep {$_ ne $ref && $_ !~ /\*/} @bits;

                # construct allele_string
                if(scalar keys %non_ref) {
                    $vf_copy->{allele_string} = $ref."/".(join "/", keys %non_ref);
                }
                else {
                    $vf_copy->{allele_string} = $ref;
                    $vf_copy->{hom_ref} = 1;

                    if(defined($config->{process_ref_homs})) {
                        $vf_copy->{allele_string} .= "/".$ref;
                    }
                    else {
                        $vf_copy->{non_variant} = 1;
                    }
                }

                # store phasing info
                $vf_copy->{phased} = defined($config->{phased} ? 1 : $phased);

                # store GT
                $vf_copy->{genotype} = \@bits;

                # store individual name
                $vf_copy->{individual} = $ind;

                push @return, $vf_copy;
            }

            return \@return;
        }
        else {
            return [$vf];
        }
    }
}

# parse a line of pileup input into variation feature objects
sub parse_pileup {
    my $config = shift;
    my $line = shift;

    my @data = split /\s+/, $line;

    # simple validity check
    unless($data[0] && $data[1] && $data[2]) {
      warning_msg($config, "Invalid input formatting on line ".$config->{line_number});
      return [];
    }

    # pileup can produce more than one VF per line
    my @return;

    # normal variant
    if($data[2] ne "*"){
        my $var;

        if($data[3] =~ /^[A|C|G|T]$/) {
            $var = $data[3];
        }
        else {
            ($var = (unambiguity_code($data[3]) || $data[3])) =~ s/$data[2]//ig;
        }

        for my $alt(split //, $var){
            push @return, Bio::EnsEMBL::Variation::VariationFeature->new_fast({
                start          => $data[1],
                end            => $data[1],
                allele_string  => $data[2].'/'.$alt,
                strand         => 1,
                map_weight     => 1,
                adaptor        => $config->{vfa},
                chr            => $data[0],
            });
        }
    }

    # in/del
    else {
        my %tmp_hash = map {$_ => 1} split '\/', $data[3];
        my @genotype = keys %tmp_hash;

        foreach my $allele(@genotype){
            if(substr($allele,0,1) eq "+") { #ins
                push @return, Bio::EnsEMBL::Variation::VariationFeature->new_fast({
                    start          => $data[1] + 1,
                    end            => $data[1],
                    allele_string  => '-/'.substr($allele, 1),
                    strand         => 1,
                    map_weight     => 1,
                    adaptor        => $config->{vfa},
                    chr            => $data[0],
                });
            }
            elsif(substr($allele,0,1) eq "-"){ #del
                push @return, Bio::EnsEMBL::Variation::VariationFeature->new_fast({
                    start          => $data[1] + 1,
                    end            => $data[1] + length(substr($allele, 1)),
                    allele_string  => substr($allele, 1).'/-',
                    strand         => 1,
                    map_weight     => 1,
                    adaptor        => $config->{vfa},
                    chr            => $data[0],
                });
            }
            elsif($allele ne "*"){
                warning_msg($config, "WARNING: invalid pileup indel genotype: $line\n");
            }
        }
    }

    return \@return;
}

# parse a line of HGVS input into a variation feature object
sub parse_hgvs {
    my $config = shift;
    my $line = shift;

    # simple validity check
    unless($line) {
      warning_msg($config, "Invalid input formatting on line ".$config->{line_number});
      return [];
    }

    my $vf;

    # not all hgvs notations are supported yet, so we have to wrap it in an eval
    eval { $vf = $config->{vfa}->fetch_by_hgvs_notation($line, $config->{sa}, $config->{ta}) };

    if((!defined($vf) || (defined $@ && length($@) > 1)) && defined($config->{coordinator})) {
        eval { $vf = $config->{vfa}->fetch_by_hgvs_notation($line, $config->{ofsa}, $config->{ofta}) };
    }

    if(!defined($vf) || (defined $@ && length($@) > 1)) {
        warning_msg($config, "WARNING: Unable to parse HGVS notation \'$line\'\n$@");
        return [];
    }

    # get whole chromosome slice
    my $slice = $vf->slice->adaptor->fetch_by_region($vf->slice->coord_system->name, $vf->slice->seq_region_name);

    $vf = $vf->transfer($slice);

    # name it after the HGVS
    $vf->{variation_name} = $line;

    # add chr attrib
    $vf->{chr} = $vf->slice->seq_region_name;

    return [$vf];
}

# parse a variation identifier e.g. a dbSNP rsID
sub parse_id {
    my $config = shift;
    my $line = shift;

    # simple validity check
    unless($line) {
      warning_msg($config, "Invalid input formatting on line ".$config->{line_number});
      return [];
    }

    # tell adaptor to fetch failed variants
    # but store state to restore afterwards
    my $prev = $config->{va}->db->include_failed_variations;
    $config->{va}->db->include_failed_variations(1);

    my $v_obj = $config->{va}->fetch_by_name($line);

    return [] unless defined $v_obj;

    my @vfs = @{$v_obj->get_all_VariationFeatures};
    for(@vfs) {
      delete $_->{dbID};
      delete $_->{overlap_consequences};
      $_->{chr} = $_->seq_region_name;
      $config->{slice_cache}->{$_->{chr}} = $_->slice;
      $_->{variation_name} = $line;
    }

    # restore state
    $config->{va}->db->include_failed_variations($prev);

    return \@vfs;
}



# converts to VCF format
sub convert_to_vcf {
    my $config = shift;
    my $vf = shift;

    # look for imbalance in the allele string
    if($vf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
        my %allele_lengths;
        my @alleles = split '\/', $vf->allele_string;

        map {reverse_comp(\$_)} @alleles if $vf->strand < 0;

        foreach my $allele(@alleles) {
            $allele =~ s/\-//g;
            $allele_lengths{length($allele)} = 1;
        }

        # in/del/unbalanced
        if(scalar keys %allele_lengths > 1) {

            # we need the ref base before the variation
            # default to N in case we can't get it
            my $prev_base = 'N';

            if(defined($vf->slice) && ref($vf->slice) eq 'Bio::EnsEMBL::Slice') {
                my $slice = $vf->slice->sub_Slice($vf->start - 1, $vf->start - 1);
                $prev_base = $slice->seq if defined($slice);
            }

            for my $i(0..$#alleles) {
                $alleles[$i] =~ s/\-//g;
                $alleles[$i] = $prev_base.$alleles[$i];
            }

            return [
                $vf->{chr} || $vf->seq_region_name,
                $vf->start - 1,
                $vf->variation_name,
                shift @alleles,
                (join ",", @alleles),
                '.', '.', '.'
            ];

        }

        # balanced sub
        else {
            return [
                $vf->{chr} || $vf->seq_region_name,
                $vf->start,
                $vf->variation_name,
                shift @alleles,
                (join ",", @alleles),
                '.', '.', '.'
            ];
        }
    }

    # SV
    else {

        # convert to SO term
        my %terms = (
            'insertion' => 'INS',
            'deletion' => 'DEL',
            'tandem_duplication' => 'TDUP',
            'duplication' => 'DUP'
        );

        my $alt = '<'.($terms{$vf->class_SO_term} || $vf->class_SO_term).'>';

        return [
            $vf->{chr} || $vf->seq_region_name,
            $vf->start,
            $vf->variation_name,
            '.',
            $alt,
            '.', '.', '.'
        ];
    }
}


# tries to map a VF to the LRG coordinate system
sub add_lrg_mappings {
    my $config = shift;
    my $vfs = shift;

    my @new_vfs;

    foreach my $vf(@$vfs) {

        # add the unmapped VF to the array
        push @new_vfs, $vf;

        # make sure the VF has an attached slice
        $vf->{slice} ||= get_slice($config, $vf->{chr}, undef, 1);
        next unless defined($vf->{slice});

        # transform LRG <-> chromosome
        my $new_vf = $vf->transform($vf->{slice}->coord_system->name eq 'lrg' ? 'chromosome' : 'lrg');

        # add it to the array if transformation worked
        if(defined($new_vf)) {

            # update new VF's chr entry
            $new_vf->{chr} = $new_vf->seq_region_name;
            push @new_vfs, $new_vf;
        }
    }

    return \@new_vfs;
}

sub minimise_alleles {
  my $config = shift;
  my $vfs = shift;

  my @new_vfs;

  foreach my $vf(@$vfs) {

    # skip VFs with more than one alt
    # they get taken care of later by split_variants/rejoin_variants
    if(!$vf->{allele_string} || $vf->{allele_string} =~ /.+\/.+\/.+/ || $vf->{allele_string} !~ /.+\/.+/) {
      push @new_vfs, $vf;
    }

    else {
      my @alleles = split('/', $vf->{allele_string});
      my $ref = shift @alleles;
      my $changed = 0;

      foreach my $alt(@alleles) {

        my $start = $vf->{start};
        my $end   = $vf->{end};

        # trim from left
        while($ref && $alt && substr($ref, 0, 1) eq substr($alt, 0, 1)) {
          $ref = substr($ref, 1);
          $alt = substr($alt, 1);
          $start++;
          $changed = 1;
        }

        # trim from right
        while($ref && $alt && substr($ref, -1, 1) eq substr($alt, -1, 1)) {
          $ref = substr($ref, 0, length($ref) - 1);
          $alt = substr($alt, 0, length($alt) - 1);
          $end--;
          $changed = 1;
        }

        $ref ||= '-';
        $alt ||= '-';

        # create a copy
        my $new_vf;
        %$new_vf = %{$vf};
        bless $new_vf, ref($vf);

        # give it a new allele string and coords
        $new_vf->{allele_string}          = $ref.'/'.$alt;
        $new_vf->{start}                  = $start;
        $new_vf->{end}                    = $end;
        $new_vf->{original_allele_string} = $vf->{allele_string};
        $new_vf->{original_start}         = $vf->{start};
        $new_vf->{original_end}           = $vf->{end};
        $new_vf->{minimised}              = 1;

        push @new_vfs, $new_vf;
      }
    }
  }

  return \@new_vfs;
}


# wrapper for whole_genome_fetch and vf_to_consequences
# takes config and a listref of VFs, returns listref of line hashes for printing
sub get_all_consequences {
    my $config     = shift;
    my $listref    = shift;

    # set some package variables to optimal values for speed
    # don't assert refs
    my $assertions_bak = $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;
    $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0;

    # don't use rearrange
    my $no_rearrange_bak = $Bio::EnsEMBL::Utils::Argument::NO_REARRANGE;
    $Bio::EnsEMBL::Utils::Argument::NO_REARRANGE = 1;

    # avoid using transfer
    my $no_transfer_bak = $Bio::EnsEMBL::Variation::TranscriptVariationAllele::NO_TRANSFER;
    $Bio::EnsEMBL::Variation::TranscriptVariationAllele::NO_TRANSFER = 1;

    if ($config->{extra}) {
        eval "use Plugin qw($config);"
    }

    # check we can use MIME::Base64
    if(defined($config->{fork})) {
        eval q{ use MIME::Base64; };

        if($@) {
            warning_msg($config, "WARNING: Unable to load MIME::Base64, forking disabled");
            delete $config->{fork};
        }
    }

    # log sorted order for VCF input
    if(($config->{format} && $config->{format} eq 'vcf') || defined($config->{rest})) {
      my $i = 0;
      $_->{_order} = sprintf("%09d", ++$i) for @$listref;
    }

    # cache link to config on VFs
    # this allows some stuff to be cached on here
    # which would otherwise be cached on adaptors that don't exist in offline mode
    for(@$listref) {
      $_->{config} = $config;
      weaken($_->{config});
    }

    my (@temp_array, @return, %by_pid, @pids);
    my $active_forks = 0;

    if(defined($config->{fork})) {
      my $total_size = scalar @$listref;
      my $done_vars = 0;

      # this variable stores child process references
      my $sel = IO::Select->new;

      debug("Calculating consequences") unless defined($config->{quiet});
      progress($config, 0, 1);

      my $delta = 0.5;
      my $minForkSize = 50;
      my $maxForkSize = int(($config->{buffer_size} || 5000) / (2 * $config->{fork}));

      # loop while variants in $listref or forks running
      while (scalar @$listref or $active_forks ) {

        # only spawn new forks if we have space
        if ($active_forks <= $config->{fork} ) {
          my $numLines = scalar @$listref;
          my $forkSize = int($numLines / ($config->{fork} + $delta*$config->{fork}) + $minForkSize ) + 1;

          $forkSize =  $maxForkSize if $forkSize > $maxForkSize;

          while($active_forks <= $config->{fork} && scalar @$listref) {
              my $tmp_vf = shift @$listref;

              push @temp_array, $tmp_vf;

              # fork
              if(scalar @temp_array >= $forkSize || scalar @$listref == 0) {
                my ($child,$parent);

                socketpair($child, $parent, AF_UNIX, SOCK_STREAM, PF_UNSPEC) or die "ERROR: Failed to open socketpair: $!";
                $child->autoflush(1);
                $parent->autoflush(1);
                $sel->add($child);

                my $pid = fork;

                if(!defined($pid)) {
                  die("ERROR: Failed to fork -") unless defined($config->{quiet});
                  last;
                }
                elsif($pid) {
                  push @pids, $pid;
                  $active_forks++;
                  @temp_array = ();
                }
                elsif($pid == 0) {
                  $config->{forked} = $$;
                  $config->{quiet} = 1;
                  $config->{stats} = {};

                  *PARENT = $parent;
                  *STDERR = *PARENT;

                  # reinitialise FASTA if using FAIDX
                  # the XS code doesn't seem to like being forked
                  if($config->{fasta_db} && $config->{fasta_db}->isa('Bio::DB::HTS::Faidx')) {
                    delete($config->{fasta_db});
                    $config->{fasta_db} = Bio::DB::HTS::Faidx->new($config->{fasta});

                    no warnings 'once';
                    $Bio::EnsEMBL::Slice::fasta_db = $config->{fasta_db};
                  }

                  #die("TEST DEATH\n") if rand() < 0.1;

                  my $cons = vf_list_to_cons($config, \@temp_array);

                  # what we're doing here is sending a serialised hash of the
                  # results through to the parent process through the socket.
                  # This is then thawed by the parent process.
                  # $$, or the PID, is added so that the input can be sorted
                  # back into the correct order for output

                  print  PARENT $$." ".encode_base64(freeze($_), "\t")."\n" for grep {defined($_)} @$cons;

                  # some plugins may cache stuff, check for this and try and
                  # reconstitute it into parent's plugin cache
                  foreach my $plugin(@{$config->{plugins}}) {
                    next unless defined($plugin->{has_cache});

                    # delete unnecessary stuff and stuff that can't be serialised
                    delete $plugin->{$_} for qw(config feature_types variant_feature_types version feature_types_wanted variant_feature_types_wanted params);
                    print PARENT $$." PLUGIN ".ref($plugin)." ".encode_base64(freeze($plugin), "\t")."\n";
                  }

                  # tell parent about stats
                  print PARENT $$." STATS ".encode_base64(freeze($config->{stats}), "\t")."\n" if defined($config->{stats});

                  # we need to tell the parent this child is finished
                  # otherwise it keeps listening
                  print PARENT "DONE $$\n";

                  exit(0);
                }
              }
          }
        }

        # read child input
        while(my @ready = $sel->can_read()) {
            my $no_read = 1;

            foreach my $fh(@ready) {
                $no_read++;

                my $line = $fh->getline();
                next unless defined($line) && $line;

                $no_read = 0;

                # child finished
                if($line =~ /^DONE/) {
                    $sel->remove($fh);
                    $fh->close;
                    $active_forks--;
                    last;
                }

                # variant finished / progress indicator
                elsif($line =~ /^BUMP/) {
                    $line =~ m/BUMP ?(\d*)/;
                    $done_vars += $1 || 1;
                    progress($config, $done_vars, $total_size);
                }

                # output
                elsif($line =~ /^\-?\d+ /) {

                    # plugin
                    if($line =~ /^\-?\d+ PLUGIN/) {

                        $line =~ m/^(\-?\d+) PLUGIN (\w+) /;
                        my ($pid, $plugin) = ($1, $2);

                        # remove the PID
                        $line =~ s/^\-?\d+ PLUGIN \w+ //;
                        chomp $line;

                        my $tmp = thaw(decode_base64($line));

                        next unless defined($plugin);

                        # copy data to parent plugin
                        my ($parent_plugin) = grep {ref($line) eq $plugin} @{$config->{plugins}};

                        next unless defined($parent_plugin);

                        merge_hashes($parent_plugin, $tmp);
                    }

                    # filtered count
                    elsif($line =~ /^\-?\d+ STATS/) {
                        $line =~ s/^\-?\d+\sSTATS\s//;
                        my $tmp = thaw(decode_base64($line));
                        $config->{stats} ||= {};

                        # special case chr lengths
                        my %chr_lengths;

                        if(defined($config->{stats}->{chr_lengths})) {
                          merge_hashes($config->{stats}->{chr_lengths}, $tmp->{chr_lengths});
                          %chr_lengths = %{$config->{stats}->{chr_lengths}};
                        }

                        merge_hashes($config->{stats}, $tmp, 1);

                        $config->{stats}->{chr_lengths} = \%chr_lengths;
                    }

                    else {
                        # grab the PID
                        $line =~ m/^(\-?\d+)\s/;
                        my $pid = $1;
                        die "ERROR: Could not parse forked PID from line $line" unless defined($pid);

                        # remove the PID
                        $line =~ s/^\-?\d+\s//;
                        chomp $line;

                        # decode and thaw "output" from forked process
                        push @{$by_pid{$pid}}, thaw(decode_base64($line));
                    }
                }

                # something's wrong
                else {
                    print STDERR "\n$line\n";
                }
            }

            # read-through detected, DIE
            die("\nERROR: Forked process(es) died\n") if $no_read;

            last if $active_forks < $config->{fork};
        }
      }

      end_progress($config);

      debug("Writing output") unless defined($config->{quiet});

      waitpid($_, 0) for @pids;

      # add the sorted data to the return array
      push @return, @{$by_pid{$_} || []} for @pids;
    }

    # no forking
    else {
        push @return, @{vf_list_to_cons($config, $listref)};
    }

    if(defined($config->{debug})) {
        eval q{use Devel::Size qw(total_size)};
        my $mem = memory();
        my $tot;
        $tot += $_ for @$mem;

        if($tot > 1000000) {
            $tot = sprintf("%.2fGB", $tot / (1024 * 1024));
        }

        elsif($tot > 1000) {
            $tot = sprintf("%.2fMB", $tot / 1024);
        }

        my $mem_diff = mem_diff($config);
        debug(
            "LINES ", $config->{line_number},
            "\tMEMORY $tot ", (join " ", @$mem),
            "\tDIFF ", (join " ", @$mem_diff),
            "\tCONFIG ", total_size($config)
        );
        #exit(0) if grep {$_ < 0} @$mem_diff;
    }

    # check and order
    my $test = $return[0];
    if(defined($test) && ref($test) ne 'HASH' && $$test =~ /^\#\#\#ORDER\#\#\#/) {
      @return = sort {$$a cmp $$b} @return;
      $$_ =~ s/\#\#\#ORDER\#\#\# \d+ // for @return;
    }
    elsif(defined($test) && ref($test) eq 'HASH' && defined($test->{_order})) {
      @return = map {delete $_->{_order}; $_} sort {$a->{_order} cmp $b->{_order}} @return;
    }

    # restore changed package variable values
    $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = $assertions_bak;
    $Bio::EnsEMBL::Utils::Argument::NO_REARRANGE = $no_rearrange_bak;
    $Bio::EnsEMBL::Variation::TranscriptVariationAllele::NO_TRANSFER = $no_transfer_bak;

    return \@return;
}

sub vf_list_to_cons {
    my $config = shift;
    my $listref = shift;

    # initialize caches
    $config->{$_.'_cache'} ||= {} for qw(tr rf slice);

    # build hash
    my %vf_hash;
    push @{$vf_hash{$_->{chr}}{int($_->{start} / $config->{chunk_size})}{$_->{start}}}, $_ for @$listref;

    # get chr list
    my @chrs = sort {natural_sort($a,$b)} keys %{{map {$_->{chr} => 1} @$listref}};

    # get non-variants
    my @non_variants = grep {$_->{non_variant}} @$listref;

    # check existing VFs
    if(defined($config->{check_existing})) {
      if(defined($config->{'cache_var_type'}) && $config->{'cache_var_type'} eq 'tabix' && !defined($config->{database})) {
        if($CAN_USE_TABIX_PM) {
          check_existing_tabix_pm($config, $listref);
        }
        elsif($CAN_USE_TABIX_CL) {
          check_existing_tabix($config, $listref);
        }
        else {
          die("ERROR: tabix does not seem to be in your path\n");
        }
      }
      else {
        check_existing_hash($config, \%vf_hash);
      }
    }

    my $new_listref = [];

    # skip any based on frequency checks?
    if(defined($config->{check_frequency})) {
        foreach my $vf(@$listref) {
            if(defined($vf->{existing}) && scalar @{$vf->{existing}}) {
                my @passed = grep {$_} map {check_frequencies($config, $_)} reverse @{$vf->{existing}};
                push @$new_listref, $vf if scalar @passed == scalar @{$vf->{existing}};
                $vf->{freqs} = $config->{filtered_freqs};
            }
            else {
                push @$new_listref, $vf;
            }
        }
    }
    else {
        $new_listref = $listref;
    }

    # if using consequence filter, we're not interested in how many remain yet
    $config->{stats}->{filter_count} += scalar @$new_listref unless defined($config->{filter}) || defined($config->{no_stats});

    # get overlapping SVs
    &check_svs_hash($config, \%vf_hash) if defined($config->{check_svs});

    # split variants with complex allele strings to get minimal reps
    my $rejoin_required = 0;

    if(defined($config->{minimal})) {
      my $before_split = scalar @$new_listref;
      $new_listref = split_variants($config, $new_listref);
      $rejoin_required = scalar @$new_listref == $before_split ? 0 : 1;
    }

    # remake hash without non-variants
    %vf_hash = ();
    push @{$vf_hash{$_->{chr}}{int($_->{start} / $config->{chunk_size})}{$_->{start}}}, $_ for grep {!defined($_->{non_variant})} @$new_listref;

    # get regions
    my $regions = &regions_from_hash($config, \%vf_hash);
    my $trim_regions = $regions;

    # prune caches
    if(!defined($config->{forked})) {
      prune_cache($config, $config->{tr_cache}, $regions, $config->{loaded_tr});
      prune_cache($config, $config->{rf_cache}, $regions, $config->{loaded_rf});
    }

    my $fetched_tr_count = 0;
    $fetched_tr_count = fetch_transcripts($config, $regions, $trim_regions)
        unless defined($config->{no_consequences});

    my $fetched_rf_count = 0;
    $fetched_rf_count = fetch_regfeats($config, $regions, $trim_regions)
        if defined($config->{regulatory})
        && !defined($config->{no_consequences});

    my @return;

    foreach my $chr(@chrs) {
        my $finished_vfs = whole_genome_fetch($config, $chr, \%vf_hash);

        # non-variants?
        if(scalar @non_variants) {
            push @$finished_vfs, grep {$_->{chr} eq $chr} @non_variants;

            # need to re-sort
            @$finished_vfs = sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @$finished_vfs;
        }

        # rejoin required?
        $finished_vfs = rejoin_variants($config, $finished_vfs) if $rejoin_required;

        debug("Calculating consequences") unless defined($config->{quiet});

        my $vf_count = scalar @$finished_vfs;
        my $vf_counter = 0;

        while(my $vf = shift @$finished_vfs) {
            progress($config, $vf_counter++, $vf_count) unless $vf_count == 1;

            # GVF output
            if(defined($config->{gvf})) {
                $vf->source(Bio::EnsEMBL::Variation::Source->new_fast({name => "User"}));

                $config->{gvf_id} ||= 1;

                # get custom annotation
                my $custom_annotation = defined($config->{custom}) ? get_custom_annotation($config, $vf) : {};
                $custom_annotation->{ID} = $config->{gvf_id}++;


                my $tmp = $vf->to_gvf(
                    include_consequences => defined($config->{no_consequences}) ? 0 : 1,
                    extra_attrs          => $custom_annotation,
                );
                push @return, \$tmp;
            }

            # VCF output
            elsif(defined($config->{vcf})) {

                # convert to VCF, otherwise get line
                my $line = $config->{format} eq 'vcf' ? [split /\s+/, $vf->{_line}] : convert_to_vcf($config, $vf);

                if(!defined($line->[7]) || $line->[7] eq '.') {
                    $line->[7] = '';
                }

                my $fieldname = $config->{vcf_info_field} || 'CSQ';

                # nuke existing CSQ field
                if($line->[7] =~ /(^|\;)$fieldname\=/ && !defined($config->{keep_csq})) {
                  $line->[7] =~ s/(^|\;)$fieldname\=\S+?(\;|$)(\S|$)/$1$3/;
                }

                # get all the lines the normal way
                # and process them into VCF-compatible string
                my $string = $fieldname.'=';

                foreach my $line_hash(grep {defined($_)} @{vf_to_consequences($config, $vf)}) {

                    # use the field list (can be user-defined by setting --fields)
                    for my $col(@{$config->{fields}}) {

                        # skip fields already represented in the VCF
                        next if $col eq 'Uploaded_variation' or $col eq 'Location' or $col eq 'Extra';

                        # search for data in main line hash as well as extra field
                        my $data = defined $line_hash->{$col} ? $line_hash->{$col} : $line_hash->{Extra}->{$col};
                        reverse_comp(\$data) if $vf->strand < 0 and $col eq 'Allele' and $config->{format} ne 'vcf';

                        # "-" means null for everything except the Allele field (confusing...)
                        $data = undef if defined($data) and $data eq '-' and $col ne 'Allele';
                        $data =~ s/\,/\&/g if defined $data;
                        $data =~ s/\;/\%3B/g if defined $data;
                        $string .= defined($data) ? $data : '';
                        $string .= '|';
                    }

                    $string =~ s/\|$//;
                    $string .= ',';
                }

                $string =~ s/\,$//;

                if(!defined($config->{no_consequences}) && $string ne $fieldname.'=') {
                    $line->[7] .= ($line->[7] ? ';' : '').$string;
                }

                # get custom annotation
                if(defined($config->{custom}) && scalar @{$config->{custom}}) {
                    my $custom_annotation = get_custom_annotation($config, $vf);
                    foreach my $key(keys %{$custom_annotation}) {
                        $line->[7] .= ($line->[7] ? ';' : '').$key.'='.$custom_annotation->{$key};
                    }
                }

                for(@$line) {
                    $_ = '.' unless defined($_);
                }

                my $tmp = join "\t", @$line;

                # add order
                $tmp = "###ORDER### ".$vf->{_order}." ".$tmp if defined($vf->{_order});

                push @return, \$tmp;
            }

            # XML output for Solr
            elsif(defined($config->{solr})) {
                eval q{
                  use CGI qw(escape);
                };

                foreach my $con(grep {defined($_)} @{vf_to_consequences($config, $vf)}) {
                    my $line = "<doc>\n";

                    # create unique ID
                    $line .= sprintf(qq{  <field name="id">%s_%i_%i_%s_%s</field>\n}, $vf->{chr}, $vf->{start}, $vf->{end}, $con->{Allele} || '-', $con->{Feature} || '-');

                    # add proper location fields that can be indexed
                    $line .= sprintf(qq{  <field name="chr">%s</field>\n}, $vf->{chr});
                    $line .= sprintf(qq{  <field name="start">%s</field>\n}, $vf->{start});
                    $line .= sprintf(qq{  <field name="end">%s</field>\n}, $vf->{end});

                    foreach my $col(@{$config->{fields}}) {

                        # search for data in main line hash as well as extra field
                        my $val = defined $con->{$col} ? $con->{$col} : $con->{Extra}->{$col};
                        next unless defined($val) && $val ne '-';

                        # some have multiple values
                        foreach my $data(split(',', $val)) {

                          # split SIFT and PolyPhen
                          if($col eq 'SIFT' || $col eq 'PolyPhen') {
                            if($data =~ m/([a-z\_]+)?\(?([\d\.]+)?\)?/i) {
                              my ($pred, $score) = ($1, $2);
                              $line .= sprintf(qq{  <field name="%s">%s</field>\n}, $col.'_pred', $pred) if $pred;
                              $line .= sprintf(qq{  <field name="%s">%s</field>\n}, $col.'_score', $score) if defined($score);
                            }
                          }

                          # GMAF
                          elsif($col eq 'GMAF') {
                            if($data =~ m/([\d\.]+)/) {
                              $line .= sprintf(qq{  <field name="%s">%s</field>\n}, $col, $1) if defined($1);
                            }
                          }

                          else {
                            $line .= sprintf(qq{  <field name="%s">%s</field>\n}, $col, escape($data)) if defined($data);
                          }
                        }
                    }

                    $line .= "</doc>\n";

                    push @return, \$line;
                }
            }

            # structured hash output for REST API
            elsif(defined($config->{rest})) {
              push @return, format_rest_output($config, $vf);
            }

            # normal output
            else {
                push @return, @{vf_to_consequences($config, $vf)};
            }

            print PARENT "BUMP\n" if defined($config->{forked}) && !defined($config->{no_progress});
        }

        end_progress($config) unless scalar @$listref == 1;
    }

    return \@return;
}

sub split_variants {
  my ($config, $listref) = @_;

  # split and link multi-allele VFs
  my @split_list;

  foreach my $original_vf(@$listref)  {
    if($original_vf->{allele_string} =~ /.+\/.+\/.+/) {

      my @alleles = split('/', $original_vf->{allele_string});
      my $original_ref = shift @alleles;
      my $first;

      my @tmp;
      my $changed = 0;
      my $base_allele_number = 1;

      foreach my $alt(@alleles) {

        my $ref   = $original_ref;
        my $start = $original_vf->{start};
        my $end   = $original_vf->{end};

        # trim from left
        while($ref && $alt && substr($ref, 0, 1) eq substr($alt, 0, 1)) {
          $ref = substr($ref, 1);
          $alt = substr($alt, 1);
          $start++;
          $changed = 1;
        }

        # trim from right
        while($ref && $alt && substr($ref, -1, 1) eq substr($alt, -1, 1)) {
          $ref = substr($ref, 0, length($ref) - 1);
          $alt = substr($alt, 0, length($alt) - 1);
          $end--;
          $changed = 1;
        }

        $ref ||= '-';
        $alt ||= '-';

        # create a copy
        my $new_vf;
        %$new_vf = %{$original_vf};
        bless $new_vf, ref($original_vf);

        # give it a new allele string and coords
        $new_vf->{allele_string} = $ref.'/'.$alt;
        $new_vf->{start} = $start;
        $new_vf->{end} = $end;
        $new_vf->{alt_allele} = $alt;

        # $new_vf->{variation_name} = 'merge_'.$alt;

        # not the first one ($first already exists)
        if($first) {
          $new_vf->{merge_with} = $first;
          $new_vf->{_base_allele_number} = $base_allele_number++;
        }

        # this is the first one
        else {
          $first = $new_vf;

          # $new_vf->{variation_name} = 'first_'.$alt;

          # store the original allele string and coords
          $first->{original_allele_string} = $original_vf->{allele_string};
          $first->{original_start}         = $original_vf->{start};
          $first->{original_end}           = $original_vf->{end};
          $first->{minimised}              = 1
        }

        push @tmp, $new_vf;
      }

      if($changed) {
        push @split_list, @tmp;
      }
      else {
        push @split_list, $original_vf;
      }
    }
    else {
      push @split_list, $original_vf;
    }
  }

  return \@split_list;
}

sub rejoin_variants {
  my ($config, $listref) = @_;

  my @joined_list = ();

  # backup stats here as methods below will increment stats counts
  my %stats_backup = %{$config->{stats} || {}};

  foreach my $vf(@$listref) {

    # reset original one
    if(defined($vf->{original_allele_string})) {

      # do consequence stuff
      vf_to_consequences($config, $vf);

      $vf->{allele_string} = $vf->{original_allele_string};
      $vf->{start}         = $vf->{original_start};
      $vf->{end}           = $vf->{original_end};

      push @joined_list, $vf;
    }

    # this one needs to be merged in
    elsif(defined($vf->{merge_with})) {
      my $original = $vf->{merge_with};

      # do consequence stuff
      vf_to_consequences($config, $vf);

      # now we have to copy the [Feature]Variation objects
      # we can't simply copy the alleles as the coords will be different
      # better to make new keys
      # we also have to set the VF pointer to the original

      # copy transcript variations etc
      foreach my $type(map {$_.'_variations'} ('transcript', map {lc_rf_type($_)} @REG_FEAT_TYPES)) {
        foreach my $key(keys %{$vf->{$type} || {}}) {
          my $val = $vf->{$type}->{$key};
          $val->base_variation_feature($original);

          # rename the key they're stored under
          $original->{$type}->{$vf->{allele_string}.'_'.$key} = $val;
        }
      }

      # intergenic variation is a bit different
      # there is only one, and no reference feature to key on
      # means we have to copy over alleles manually
      if(my $iv = $vf->{intergenic_variation}) {

        $iv->base_variation_feature($original);

        if(my $oiv = $original->{intergenic_variation}) {
            push @{$oiv->{alt_alleles}}, @{$iv->{alt_alleles}};
            $oiv->{_alleles_by_seq}->{$_->variation_feature_seq} = $_ for @{$oiv->{alt_alleles}};
        }

        # this probably won't happen, but can't hurt to cover all bases
        else {
            $original->{intergenic_variation} = $iv;
        }
      }

      # reset these keys, they can be recalculated
      delete $original->{$_} for qw(overlap_consequences _most_severe_consequence);
    }

    # normal
    else {
      push @joined_list, $vf;
    }
  }

  $config->{stats} = \%stats_backup if $config->{stats};
  return \@joined_list;
}

sub natural_sort {
  my ($a, $b) = @_;
  if($a =~ /^[0-9]+$/ && $b =~ /^[0-9]+$/) {
    return $a <=> $b;
  }
  else {
    return $a cmp $b;
  }
}

sub format_rest_output {
  my ($config, $vf) = @_;

  # define the set of fields that are lists
  my @list_fields = qw(domains cell_type consequence refseq pubmed clin_sig);

  # define some to delete
  my @delete_fields = qw(
    Location
    Uploaded_variation
    Existing_variation
    GMAF AFR_MAF AMR_MAF ASN_MAF EAS_MAF SAS_MAF EUR_MAF
    AA_MAF EA_MAF
    ExAC_MAF ExAC_Adj_MAF ExAC_AFR_MAF ExAC_AMR_MAF ExAC_EAS_MAF ExAC_FIN_MAF ExAC_NFE_MAF ExAC_OTH_MAF ExAC_SAS_MAF
    PUBMED CLIN_SIG SOMATIC VARIANT_CLASS PHENO
  );

  # define some fields to rename
  my %rename = (
    'consequence' => 'consequence_terms',
    'gene' => 'gene_id',
    'allele' => 'variant_allele',
    'symbol' => 'gene_symbol',
    'symbol_source' => 'gene_symbol_source',
    'overlapbp' => 'bp_overlap',
    'overlappc' => 'percentage_overlap',
    'refseq' => 'refseq_transcript_ids',
    'ensp' => 'protein_id',
    'chr' => 'seq_region_name',
    'variation_name' => 'id',
  );

  my $hash = {
    id => $vf->{variation_name},
    seq_region_name => $vf->{chr},
    start => $vf->{start},
    end => $vf->{end},
    strand => $vf->{strand},
    _order => $vf->{_order},
  };

  # add original input for use by POST endpoints
  $hash->{input} = $vf->{_line} if defined($vf->{_line});

  # add variant class
  $hash->{variant_class} = $vf->class_SO_term if defined($config->{variant_class});

  if(defined($vf->{allele_string})) {
    $hash->{allele_string} = $vf->{allele_string};
  }
  else {
    $hash->{variant_class} = $vf->{class_SO_term};
  }

  # record all cons terms so we can get the most severe
  my @con_terms;

  # add consequence stuff
  foreach my $con(grep {defined($_)} @{vf_to_consequences($config, $vf)}) {

    # flatten
    $con->{$_} = $con->{Extra}->{$_} for keys %{$con->{Extra}};
    delete $con->{Extra};

    # remove unwanted keys
    delete $con->{$_} for @delete_fields;

    # lc and remove empty
    foreach my $key(keys %$con) {
      my $tmp = $con->{$key};
      delete $con->{$key};

      next if !defined($tmp) || $tmp eq '-';

      # convert YES to 1
      $tmp = 1 if $tmp eq 'YES';

      # fix position fields into start and end
      if($key =~ /(\w+?)\_position$/i) {
        my $coord_type = lc($1);
        my ($s, $e) = split('-', $tmp);
        $con->{$coord_type.'_start'} = $s;
        $con->{$coord_type.'_end'} = defined($e) && $e =~ /^\d+$/ ? $e : $s;

        # on rare occasions coord can be "?"; for now just don't print anything
        delete $con->{$coord_type.'_start'} unless $con->{$coord_type.'_start'} =~ /^\d+$/;
        delete $con->{$coord_type.'_end'} unless $con->{$coord_type.'_end'} =~ /^\d+$/;
        next;
      }

      $con->{lc($key)} = $tmp;
    }

    my $ftype = lc($con->{feature_type} || 'intergenic');
    $ftype =~ s/feature/\_feature/;
    delete $con->{feature_type};

    # fix SIFT and PolyPhen
    foreach my $tool(qw(sift polyphen)) {
      if(defined($con->{$tool}) && $con->{$tool} =~ m/([a-z\_]+)?\(?([\d\.]+)?\)?/i) {
        my ($pred, $score) = ($1, $2);
        $con->{$tool.'_prediction'} = $pred if $pred;
        $con->{$tool.'_score'} = $score if defined($score);
        delete $con->{$tool};
      }
    }

    # fix comma-separated lists into arrays
    foreach my $key(grep {defined($con->{$_})} @list_fields) {
      $con->{$key} = [split(',', $con->{$key})];
    }

    # fix domains
    if(defined($con->{domains})) {
      my @dom;

      foreach(@{$con->{domains}}) {
        m/(\w+)\:(\w+)/;
        push @dom, {"db" => $1, "name" => $2} if $1 && $2;
      }
      $con->{domains} = \@dom;
    }

    # log observed consequence terms
    push @con_terms, @{$con->{consequence}};

    # rename
    $rename{feature} = lc($ftype).'_id';
    foreach my $key(grep {defined($con->{$_})} keys %rename) {
      $con->{$rename{$key}} = $con->{$key};
      delete $con->{$key};
    }

    push @{$hash->{$ftype.'_consequences'}}, $con;
  }

  # get most severe consequence from those logged in @con_terms
  my %all_cons = %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;
  $hash->{most_severe_consequence} = (sort {$all_cons{$a}->rank <=> $all_cons{$b}->rank} grep {$_ ne '?'} @con_terms)[0] || '?';

  # add assembly
  $hash->{assembly_name} = $config->{assembly} || $config->{cache_assembly};

  # add existing variants
  if(defined($vf->{existing}) && scalar @{$vf->{existing}}) {

    foreach my $ex_orig(@{$vf->{existing}}) {

      # work on a copy as we're going to modify/delete things
      my $ex;
      %$ex = %$ex_orig;

      delete $ex->{$_} for qw(failed);

      # frequencies
      foreach my $pop(grep {defined($ex->{$_})} qw(
        AFR AMR ASN EAS SAS EUR
        AA EA
        ExAC ExAC_Adj ExAC_AFR ExAC_AMR ExAC_EAS ExAC_FIN ExAC_NFE ExAC_OTH ExAC_SAS
      )) {
        my $tmp = $ex->{$pop};

        if($tmp =~ /(\w)\:([\d\.\-e]+)/) {
          $ex->{lc($pop).'_maf'} = $2;
          $ex->{lc($pop).'_allele'} = $1;
        }
        else {
          $ex->{lc($pop).'_maf'} = $tmp;
        }

        delete $ex->{$pop};
      }

      # remove empty
      foreach my $key(keys %$ex) {
        delete $ex->{$key} if !defined($ex->{$key}) || $ex->{$key} eq '' || ($key !~ /maf/ && $ex->{$key} eq 0);
      }

      # fix comma-separated lists into arrays
      foreach my $key(grep {defined($ex->{$_})} @list_fields) {
        $ex->{$key} = [split(',', $ex->{$key})];
      }

      # rename
      foreach my $key(grep {defined($ex->{$_})} keys %rename) {
        $ex->{$rename{$key}} = $ex->{$key};
        delete $ex->{$key};
      }

      push @{$hash->{colocated_variants}}, $ex;
    }
  }

  numberify($hash);

  return $hash;
}

sub numberify {
  my $ref = shift;

  if(ref($ref) eq 'HASH') {
    foreach my $k(keys %$ref) {
      if(ref($ref->{$k}) =~ /HASH|ARRAY/) {
        numberify($ref->{$k});
      }
      else {
        $ref->{$k} = $ref->{$k} + 0 if defined($ref->{$k}) && $k ne 'seq_region_name' && $k ne 'id' && looks_like_number($ref->{$k});
      }
    }
  }
  elsif(ref($ref) eq 'ARRAY') {
    foreach my $i(0..((scalar @$ref) - 1)) {
      if(ref($ref->[$i]) =~ /HASH|ARRAY/) {
        numberify($ref->[$i]);
      }
      else {
        $ref->[$i] = $ref->[$i] + 0 if defined($ref->[$i]) && looks_like_number($ref->[$i]);
      }
    }
  }
}

# takes a variation feature and returns ready to print consequence information
sub vf_to_consequences {
  my $config = shift;
  my $vf = shift;

  # force empty hash into object's transcript_variations if undefined from whole_genome_fetch
  # this will stop the API trying to go off and fill it again
  if(defined $config->{whole_genome}) {
    $vf->{$_} ||= {} for map {$_.'_variations'} ('transcript', map {lc_rf_type($_)} @REG_FEAT_TYPES);
  }

  my $vf_ref = ref($vf);

  # prefetch intergenic variation
  # pass a true argument to get_IntergenicVariation to stop it doing a reference allele check
  # (to stay consistent with the rest of the VEP)
  $vf->get_IntergenicVariation(1) if $vf_ref eq 'Bio::EnsEMBL::Variation::VariationFeature';

  # stats
  unless(defined($config->{no_stats})) {

    # position
    $config->{stats}->{chr}->{$vf->{chr}}->{1e6 * int($vf->start / 1e6)}++;

    # most severe consequence
    $config->{stats}->{var_cons}->{$vf->display_consequence}++;

    # known variants
    $config->{stats}->{existing}++ if defined($vf->{existing}) && scalar @{$vf->{existing}};

    # get stats
    my $so_term = $vf->class_SO_term;
    if(defined($so_term)) {
      $config->{stats}->{classes}->{$so_term}++;
      $config->{stats}->{allele_changes}->{$vf->allele_string}++ if $so_term eq 'SNV';
    }
  }

  # use a different method for SVs
  return svf_to_consequences($config, $vf) if $vf_ref eq 'Bio::EnsEMBL::Variation::StructuralVariationFeature';

  my @return = ();
  return \@return if defined($config->{no_intergenic}) && defined($vf->{intergenic_variation});

  # method name stub for getting *VariationAlleles
  my $allele_method = defined($config->{process_ref_homs}) ? 'get_all_' : 'get_all_alternate_';

  # get all VFOAs
  # need to be sensitive to whether --regulatory or --coding_only is switched on
  my $vfos;
  my $method = $allele_method.'VariationFeatureOverlapAlleles';

  # include regulatory stuff?
  if(!defined $config->{coding_only} && defined $config->{regulatory}) {
    $vfos = $vf->get_all_VariationFeatureOverlaps;
  }
  # otherwise just get transcript & intergenic ones
  else {
    @$vfos = grep {defined($_)} (
      @{$vf->get_all_TranscriptVariations},
      $vf->get_IntergenicVariation
    );
  }

  # grep out non-coding?
  @$vfos = grep {$_->can('affects_cds') && $_->affects_cds} @$vfos if defined($config->{coding_only});

  # get alleles
  my @vfoas = map {@{$_->$method}} @{$vfos};

  # only most severe or summary?
  # these options don't need any feature-specific columns
  if(defined($config->{most_severe}) || defined($config->{summary})) {

    my $line = init_line($config, $vf);

    my $term_method = $config->{terms}.'_term';

    my @ocs = sort {$a->rank <=> $b->rank} map {@{$_->get_all_OverlapConsequences}} @vfoas;

    if(@ocs) {

      # summary is just all unique consequence terms
      if(defined($config->{summary})) {
        $line->{Consequence} = join ",", keys %{{map {$_ => 1} map {$_->$term_method || $_->SO_term} @ocs}};
      }

      # most severe is the consequence term with the lowest rank
      else {
        $line->{Consequence} = $ocs[0]->$term_method || $ocs[0]->SO_term;
      }

      unless(defined($config->{no_stats})) {
        $config->{stats}->{consequences}->{$_}++ for split(',', $line->{Consequence});
      }
    }
    else {
      warning_msg(
        $config,
        "Unable to assign consequence type on line ".$vf->{_line_number}
      );
    }

    push @return, $line;
  }

  # otherwise do normal consequence processing
  else {

    # pick worst?
    if(defined($config->{pick})) {
      @vfoas = (pick_worst_vfoa($config, \@vfoas)) if defined($config->{pick});
    }

    # pick worst per allele?
    elsif(defined($config->{pick_allele})) {
      my %by_allele;
      push @{$by_allele{$_->variation_feature_seq}}, $_ for @vfoas;
      @vfoas = ();
      push @vfoas, pick_worst_vfoa($config, $by_allele{$_}) for keys %by_allele;
    }

    # flag picked?
    elsif(defined($config->{flag_pick})) {
      if(my $worst = pick_worst_vfoa($config, \@vfoas)) {
        $worst->{PICK} = 1;
      }
    }

    # flag worst per allele?
    elsif(defined($config->{flag_pick_allele})) {
      my %by_allele;
      push @{$by_allele{$_->variation_feature_seq}}, $_ for @vfoas;
      pick_worst_vfoa($config, $by_allele{$_})->{PICK} = 1 for keys %by_allele;
    }

    # pick per gene?
    elsif(defined($config->{per_gene})) {
      @vfoas = @{pick_vfoa_per_gene($config, \@vfoas)};
    }

    # process remaining
    push @return, map {vfoa_to_line($config, $_)} grep {defined($_)} @vfoas;
  }

  return \@return;
}


# picks the worst of a list of VariationFeatureOverlapAlleles
# VFOAs are ordered by a heirarchy:
# 1: canonical
# 2: transcript support level
# 3: biotype (protein coding favoured)
# 4: consequence rank
# 5: transcript length
# 6: transcript from Ensembl?
# 7: transcript from RefSeq?
sub pick_worst_vfoa {
  my $config = shift;
  my $vfoas = shift;

  my @vfoa_info;
  my %ranks = map {$_->SO_term => $_->rank} values %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;

  return $vfoas->[0] if scalar @$vfoas == 1;

  foreach my $vfoa(@$vfoas) {

    # create a hash of info for this VFOA that will be used to rank it
    my $info = {
      vfoa => $vfoa,
      rank => undef,

      # these will only be used by transcript types, default to 1 for others
      # to avoid writing an else clause below
      canonical => 1,
      ccds => 1,
      length => 0,
      biotype => 1,
      tsl => 100,
      appris => 100,
      ensembl => 1,
      refseq => 1,
    };

    if($vfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {
      my $tr = $vfoa->feature;

      # 0 is "best"
      $info->{canonical} = $tr->is_canonical ? 0 : 1;
      $info->{biotype} = $tr->biotype eq 'protein_coding' ? 0 : 1;
      $info->{ccds} = $tr->{_ccds} && $tr->{_ccds} ne '-' ? 0 : 1;
      $info->{lc($tr->{_source_cache})} = 0 if exists($tr->{_source_cache});

      # "invert" length so longer is best
      $info->{length} = 0 - (
        $tr->translation ?
        length($tr->{_variation_effect_feature_cache}->{translateable_seq} || $tr->translateable_seq) :
        $tr->length()
      );

      # lower TSL is best
      if(my ($tsl) = @{$tr->get_all_Attributes('TSL')}) {
        if($tsl->value =~ m/tsl(\d+)/) {
          $info->{tsl} = $1 if $1;
        }
      }

      # lower APPRIS is best
      if(my ($appris) = @{$tr->get_all_Attributes('appris')}) {
        if($appris->value =~ m/([A-Za-z]).+(\d+)/) {
          my ($type, $grade) = ($1, $2);

          # values are principal1, principal2, ..., alternative1, alternative2
          # so add 10 to grade if alternate
          $grade += 10 if substr($type, 0, 1) eq 'a';

          $info->{appris} = $grade if $grade;
        }
      }
    }

    push @vfoa_info, $info;
  }

  if(scalar @vfoa_info) {
    my @order = defined($config->{pick_order}) ? @{$config->{pick_order}} : @PICK_ORDER;
    my $picked;

    # go through each category in order
    foreach my $cat(@order) {

      # get ranks here as it saves time
      if($cat eq 'rank') {
        foreach my $info(@vfoa_info) {
          my @ocs = sort {$a->rank <=> $b->rank} @{$info->{vfoa}->get_all_OverlapConsequences};
          $info->{rank} = scalar @ocs ? $ranks{$ocs[0]->SO_term} : 1000;
        }
      }

      # sort on that category
      @vfoa_info = sort {$a->{$cat} <=> $b->{$cat}} @vfoa_info;

      # take the first (will have the lowest value of $cat)
      $picked = shift @vfoa_info;
      my @tmp = ($picked);

      # now add to @tmp those vfoas that have the same value of $cat as $picked
      push @tmp, shift @vfoa_info while @vfoa_info && $vfoa_info[0]->{$cat} eq $picked->{$cat};

      # if there was only one, return
      return $picked->{vfoa} if scalar @tmp == 1;

      # otherwise shrink the array to just those that had the lowest
      # this gives fewer to sort on the next round
      @vfoa_info = @tmp;
    }

    # probably shouldn't get here, but if we do, return the first
    return $vfoa_info[0]->{vfoa};
  }

  return undef;
}

# pick one vfoa per gene
# allow non-transcript types to pass through
sub pick_vfoa_per_gene {
  my $config = shift;
  my $vfoas = shift;

  my @return;
  my @tvas;

  # pick out TVAs
  foreach my $vfoa(@$vfoas) {
    if($vfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {
      push @tvas, $vfoa;
    }
    else {
      push @return, $vfoa;
    }
  }

  # sort the TVA objects into a hash by gene
  my %by_gene;

  foreach my $tva(@tvas) {
    my $gene = $tva->transcript->{_gene_stable_id} || $config->{ga}->fetch_by_transcript_stable_id($tva->transcript->stable_id)->stable_id;
    push @{$by_gene{$gene}}, $tva;
  }

  foreach my $gene(keys %by_gene) {
    push @return, grep {defined($_)} pick_worst_vfoa($config, $by_gene{$gene});
  }

  return \@return;
}

# get consequences for a structural variation feature
sub svf_to_consequences {
    my $config = shift;
    my $svf    = shift;

    my @return = ();

    # stats
    $config->{stats}->{classes}->{$svf->{class_SO_term}}++ unless defined($config->{no_stats});

    my $term_method = $config->{terms}.'_term';

    if(defined $config->{whole_genome}) {
        $svf->{transcript_structural_variations} ||= [];
        $svf->{regulation_structural_variations}->{$_} ||= [] for @REG_FEAT_TYPES;
    }

    if ((my $iv = $svf->get_IntergenicStructuralVariation(1)) && !defined($config->{no_intergenic})) {

        for my $iva (@{ $iv->get_all_alternate_IntergenicStructuralVariationAlleles }) {

            my $line = init_line($config, $svf);

            $line->{Allele} = '-';

            my $cons = $iva->get_all_OverlapConsequences->[0];

            $line->{Consequence} = $cons->$term_method || $cons->SO_term;

            $config->{stats}->{consequences}->{$cons->$term_method || $cons->SO_term}++ unless defined($config->{no_stats});

            $line = run_plugins($iva, $line, $config);

            push @return, $line;
        }
    }

    foreach my $svo(@{$svf->get_all_StructuralVariationOverlaps}) {

        next if $svo->isa('Bio::EnsEMBL::Variation::IntergenicStructuralVariation');

        my $feature = $svo->feature;

        # get feature type
        my $feature_type = (split '::', ref($feature))[-1];

        my $base_line = {
            Feature_type     => $feature_type,
            Feature          => $feature->stable_id,
            Allele           => $svf->class_SO_term,
        };

        if($svo->isa('Bio::EnsEMBL::Variation::BaseTranscriptVariation')) {
            $base_line->{cDNA_position}    = format_coords($svo->cdna_start, $svo->cdna_end).
              (defined($config->{total_length}) ? '/'.$feature->length : '');
            $base_line->{CDS_position}     = format_coords($svo->cds_start, $svo->cds_end).
              (defined($config->{total_length}) && $feature->{_variation_effect_feature_cache}->{translateable_seq} ?
                '/'.length($feature->{_variation_effect_feature_cache}->{translateable_seq}) : ''
              );
            $base_line->{Protein_position} = format_coords($svo->translation_start, $svo->translation_end).
              (defined($config->{total_length}) && $feature->{_variation_effect_feature_cache}->{peptide} ?
                '/'.length($feature->{_variation_effect_feature_cache}->{peptide}) : ''
              );
        }

        foreach my $svoa(@{$svo->get_all_StructuralVariationOverlapAlleles}) {
            my $line = init_line($config, $svf, $base_line);

            $line->{Consequence} = join ",",
                #map {s/feature/$feature_type/e; $_}
                map {$_->$term_method}
                sort {$a->rank <=> $b->rank}
                @{$svoa->get_all_OverlapConsequences};

            unless(defined($config->{no_stats})) {
              map {$config->{stats}->{consequences}->{$_->$term_method}++} @{$svoa->get_all_OverlapConsequences};
            }

            # work out overlap amounts
            my $overlap_start  = (sort {$a <=> $b} ($svf->start, $feature->start))[-1];
            my $overlap_end    = (sort {$a <=> $b} ($svf->end, $feature->end))[0];
            my $overlap_length = ($overlap_end - $overlap_start) + 1;
            my $overlap_pc     = 100 * ($overlap_length / (($feature->end - $feature->start) + 1));

            $line->{Extra}->{OverlapBP} = $overlap_length if $overlap_length > 0;
            $line->{Extra}->{OverlapPC} = sprintf("%.2f", $overlap_pc) if $overlap_pc > 0;

            add_extra_fields($config, $line, $svoa);

            $line = run_plugins($svoa, $line, $config);

            push @return, $line;
        }
    }

    return \@return;
}

# run all of the configured plugins on a VariationFeatureOverlapAllele instance
# and store any results in the provided line hash
sub run_plugins {

    my ($bvfoa, $line_hash, $config) = @_;

    my $skip_line = 0;

    for my $plugin (@{ $config->{plugins} || [] }) {

        # check that this plugin is interested in this type of variation feature

        if ($plugin->check_variant_feature_type(ref $bvfoa->base_variation_feature)) {

            # check that this plugin is interested in this type of feature

            if ($plugin->check_feature_type(ref $bvfoa->feature || 'Intergenic')) {

                eval {
                    my $plugin_results = $plugin->run($bvfoa, $line_hash);

                    if (defined $plugin_results) {
                        if (ref $plugin_results eq 'HASH') {
                            for my $key (keys %$plugin_results) {
                                $line_hash->{Extra}->{$key} = $plugin_results->{$key};
                            }
                        }
                        else {
                            warning_msg($config, "Plugin '".(ref $plugin)."' did not return a hashref, output ignored!\n");
                        }
                    }
                    else {
                        # if a plugin returns undef, that means it want to filter out this line
                        $skip_line = 1;
                    }
                };
                if ($@) {
                    warning_msg($config, "Plugin '".(ref $plugin)."' went wrong: $@");
                }

                # there's no point running any other plugins if we're filtering this line,
                # because the first plugin to skip the line wins, so we might as well last
                # out of the loop now and avoid any unnecessary computation

                last if $skip_line;
            }
        }
    }

    return $skip_line ? undef : $line_hash;
}

# turn a generic VariationFeatureOverlapAllele into a line hash
sub vfoa_to_line {
  my $config = shift;
  my $vfoa = shift;

  # stats
  unless(defined($config->{no_stats})) {
    my $term_method = $config->{terms}.'_term';
    map {$config->{stats}->{consequences}->{$_->$term_method || $_->SO_term}++} @{$vfoa->get_all_OverlapConsequences};
  }

  my $line;

  if($vfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {
    $line = tva_to_line($config, $vfoa);
  }
  elsif($vfoa->isa('Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele')) {
    $line = rfva_to_line($config, $vfoa);
  }
  elsif($vfoa->isa('Bio::EnsEMBL::Variation::MotifFeatureVariationAllele')) {
    $line = mfva_to_line($config, $vfoa);
  }
  elsif($vfoa->isa('Bio::EnsEMBL::Variation::IntergenicVariationAllele')) {
    $line = iva_to_line($config, $vfoa);
  }
  else {
    return undef;
  }

  # warn user if no consequence type was added
  if(!$line->{Consequence}) {
    $line->{Consequence} = '?';

    warning_msg(
        $config,
        "Unable to assign consequence type on line ".$vfoa->variation_feature->{_line_number}.
        ($vfoa->can('feature') && $vfoa->feature ? " for ".$vfoa->feature->stable_id : "")
    );
  }

  # add extra fields
  $line = add_extra_fields($config, $line, $vfoa);

  # run plugins
  $line = run_plugins($vfoa, $line, $config);

  return $line;
}

# process IntergenicVariationAllele
sub iva_to_line {
  my $config = shift;
  my $iva = shift;

  my $line = init_line($config, $iva->variation_feature);

  $line->{Allele} = $iva->variation_feature_seq;

  my $cons = $iva->get_all_OverlapConsequences->[0];

  # method name for consequence terms
  my $term_method = $config->{terms}.'_term';
  $line->{Consequence} = $cons->$term_method || $cons->SO_term;

  return $line;
}

# process TranscriptVariationAllele
sub tva_to_line {
  my $config = shift;
  my $tva = shift;

  my $tv = $tva->base_variation_feature_overlap;
  my $t  = $tv->transcript;

  # method name for consequence terms
  my $term_method = $config->{terms}.'_term';

  my $csq = join ",", map {$_->$term_method || $_->SO_term} sort {$a->rank <=> $b->rank} @{$tva->get_all_OverlapConsequences};

  my $pre = $tva->_pre_consequence_predicates();

  my $base_line = {
    Feature_type     => 'Transcript',
    Feature          => (defined $t ? $t->stable_id : undef),

    # use pre_consequence_predicates to avoid calling coord methods
    cDNA_position    =>
      ($pre->{exon} ? format_coords($tv->cdna_start, $tv->cdna_end) : '-').
      (defined($config->{total_length}) ? '/'.$t->length : ''),

    CDS_position     =>
      ($pre->{coding} ? format_coords($tv->cds_start, $tv->cds_end) : '-').
      (
        defined($config->{total_length}) && $t->{_variation_effect_feature_cache}->{translateable_seq} ?
        '/'.length($t->{_variation_effect_feature_cache}->{translateable_seq}) : ''
      ),

    Protein_position =>
      ($pre->{coding} ? format_coords($tv->translation_start, $tv->translation_end) : '-').
      (
        defined($config->{total_length}) && $t->{_variation_effect_feature_cache}->{peptide} ?
        '/'.length($t->{_variation_effect_feature_cache}->{peptide}) : ''
      ),

    Allele           => $tva->variation_feature_seq,
    Amino_acids      => ($pre->{coding} ? $tva->pep_allele_string : undef),
    Codons           => ($pre->{coding} ? $tva->display_codon_allele_string : undef),
    Consequence      => $csq,
  };

  if(!defined($config->{no_stats}) && $pre->{coding} && defined($tv->translation_start)) {
    $config->{stats}->{protein_pos}->{int(10 * ($tv->translation_start / ($t->{_variation_effect_feature_cache}->{peptide} ? length($t->{_variation_effect_feature_cache}->{peptide}) : $t->translation->length)))}++;
  }

  my $line = init_line($config, $tv->base_variation_feature, $base_line);

  # HGVS
  if(defined $config->{hgvs} && $pre->{within_feature}) {
    my $hgvs_t = $tva->hgvs_transcript;
    my $hgvs_p = $tva->hgvs_protein;
    my $offset = $tva->hgvs_offset;

    # URI encode "="
    $hgvs_p =~ s/\=/\%3D/g if $hgvs_p && !(defined($config->{no_escape}) || defined($config->{json}) || defined($config->{rest}));

    $line->{Extra}->{HGVSc} = $hgvs_t if $hgvs_t;
    $line->{Extra}->{HGVSp} = $hgvs_p if $hgvs_p;

    $line->{Extra}->{HGVS_OFFSET} = $offset if $offset;
  }

  if($pre->{coding}) {
    foreach my $tool (qw(SIFT PolyPhen)) {
      my $lc_tool = lc($tool);

      if (my $opt = $config->{$lc_tool}) {
        my $want_pred  = $opt =~ /^p/i;
        my $want_score = $opt =~ /^s/i;
        my $want_both  = $opt =~ /^b/i;

        if ($want_both) {
          $want_pred  = 1;
          $want_score = 1;
        }

        next unless $want_pred || $want_score;

        my $pred_meth  = $lc_tool.'_prediction';
        my $score_meth = $lc_tool.'_score';
        my $analysis   = $config->{polyphen_analysis} if $lc_tool eq 'polyphen';

        my $pred = $tva->$pred_meth($analysis);

        if($pred) {

          if ($want_pred) {
            $pred =~ s/\s+/\_/g;
            $pred =~ s/\_\-\_/\_/g;
            $line->{Extra}->{$tool} = $pred;
          }

          if ($want_score) {
            my $score = $tva->$score_meth($analysis);

            if(defined $score) {
              if($want_pred) {
                $line->{Extra}->{$tool} .= "($score)";
              }
              else {
                $line->{Extra}->{$tool} = $score;
              }
            }
          }
        }

        # update stats
        $config->{stats}->{$tool}->{$tva->$pred_meth}++ if $tva->$pred_meth && !defined($config->{no_stats});
      }
    }
  }

  return $line;
}

# process RegulatoryFeatureVariationAllele
sub rfva_to_line {
  my $config = shift;
  my $rfva = shift;

  # method name for consequence terms
  my $term_method = $config->{terms}.'_term';

  # method name stub for getting *VariationAlleles
  my $allele_method = defined($config->{process_ref_homs}) ? 'get_all_' : 'get_all_alternate_';

  my $rf = $rfva->regulatory_feature;

  my $base_line = {
    Feature_type => 'RegulatoryFeature',
    Feature      => $rf->stable_id,
  };

  if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {
    $base_line->{Extra}->{CELL_TYPE} = join ",",
      map {$_.':'.$rf->{cell_types}->{$_}}
      grep {$rf->{cell_types}->{$_}}
      @{$config->{cell_type}};

    $base_line->{Extra}->{CELL_TYPE} =~ s/\s+/\_/g;
  }

  $base_line->{Extra}->{BIOTYPE} = ref($rf->{feature_type}) ? $rf->{feature_type}->{so_name} : $rf->{feature_type} if defined($rf->{feature_type});

  my $method = $allele_method.'RegulatoryFeatureVariationAlleles';

  my $line = init_line($config, $rfva->variation_feature, $base_line);

  $line->{Allele}         = $rfva->variation_feature_seq;
  $line->{Consequence}    = join ',',
    map { $_->$term_method || $_->SO_term }
    @{ $rfva->get_all_OverlapConsequences };

  $line = run_plugins($rfva, $line, $config);

  return $line;
}

# process MotifFeatureVariationAllele
sub mfva_to_line {
  my $config = shift;
  my $mfva = shift;

  # method name for consequence terms
  my $term_method = $config->{terms}.'_term';

  # method name stub for getting *VariationAlleles
  my $allele_method = defined($config->{process_ref_homs}) ? 'get_all_' : 'get_all_alternate_';

  my $mf = $mfva->motif_feature;

  # check that the motif has a binding matrix, if not there's not
  # much we can do so don't return anything
  return undef unless defined $mf->binding_matrix;

  my $matrix = ($mf->binding_matrix->description ? $mf->binding_matrix->description.' ' : '').$mf->display_label;
  $matrix =~ s/\s+/\_/g;

  my $base_line = {
    Feature_type => 'MotifFeature',
    Feature      => $mf->binding_matrix->name,
    Extra        => {
      MOTIF_NAME  => $matrix,
      STRAND      => $mf->strand + 0
    }
  };

  if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {
    $base_line->{Extra}->{CELL_TYPE} = join ",",
      map {$_.':'.$mf->{cell_types}->{$_}}
      grep {$mf->{cell_types}->{$_}}
      @{$config->{cell_type}};

    $base_line->{Extra}->{CELL_TYPE} =~ s/\s+/\_/g;
  }

  my $method = $allele_method.'MotifFeatureVariationAlleles';

  my $line = init_line($config, $mfva->variation_feature, $base_line);

  $line->{Extra}->{MOTIF_POS}          = $mfva->motif_start if defined $mfva->motif_start;
  $line->{Extra}->{HIGH_INF_POS}       = ($mfva->in_informative_position ? 'Y' : 'N');

  my $delta = $mfva->motif_score_delta if $mfva->variation_feature_seq =~ /^[ACGT]+$/;

  $line->{Extra}->{MOTIF_SCORE_CHANGE} = sprintf("%.3f", $delta) if defined $delta;

  $line->{Allele}         = $mfva->variation_feature_seq;
  $line->{Consequence}    = join ',',
    map { $_->$term_method || $_->SO_term }
    @{ $mfva->get_all_OverlapConsequences };

  return $line;
}

sub add_extra_fields {
    my $config = shift;
    my $line   = shift;
    my $bvfoa  = shift;

    my $ocs = $bvfoa->get_all_OverlapConsequences;

    # impact
    $line->{Extra}->{IMPACT} = (sort {$a->rank <=> $b->rank} @$ocs)[0]->impact() if scalar @$ocs;

    # overlapping SVs
    if(defined $config->{check_svs} && defined $bvfoa->base_variation_feature->{overlapping_svs}) {
        $line->{Extra}->{SV} = $bvfoa->base_variation_feature->{overlapping_svs};
    }

    # allele number
    if(defined($config->{allele_number})) {
      $line->{Extra}->{ALLELE_NUM} = $bvfoa->allele_number if $bvfoa->can('allele_number');
    }

    # add transcript-specific fields
    $line = add_extra_fields_transcript($config, $line, $bvfoa) if $bvfoa->isa('Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele');

    # picked?
    $line->{Extra}->{PICK} = 1 if defined($bvfoa->{PICK});

    # stats
    unless(defined($config->{no_stats})) {
      $config->{stats}->{gene}->{$line->{Gene}}++ if defined($line->{Gene});
      $config->{stats}->{lc($line->{Feature_type})}->{$line->{Feature}}++ if defined($line->{Feature_type}) && defined($line->{Feature});
    }

    return $line;
}

sub add_extra_fields_transcript {
    my $config = shift;
    my $line = shift;
    my $tva = shift;

    my $tv = $tva->base_variation_feature_overlap;
    my $tr = $tv->transcript;
    my $pre = $tva->_pre_consequence_predicates;

    # get gene
    $line->{Gene} = $tr->{_gene_stable_id};

    # strand
    $line->{Extra}->{STRAND} = $tr->strand + 0;

    my @attribs = @{$tr->get_all_Attributes()};

    # flags
    my @flags = grep {$_->code =~ /^cds_/} @attribs;
    $line->{Extra}->{FLAGS} = join(",", map {$_->code} @flags) if scalar @flags;

    # exon/intron numbers
    if ($config->{numbers}) {
      if($pre->{exon}) {
        if(my $num = $tv->exon_number) {
          $line->{Extra}->{EXON} = $num;
        }
      }
      if($pre->{intron}) {
        if(my $num = $tv->intron_number) {
          $line->{Extra}->{INTRON} = $num;
        }
      }
    }

    if ($config->{domains} && $pre->{coding}) {
        my $feats = $tv->get_overlapping_ProteinFeatures;

        my @strings;

        for my $feat (@$feats) {

            # do a join/grep in case of missing data
            my $label = join(':', grep {$_} ($feat->analysis->display_label, $feat->hseqname));

            # replace any special characters
            $label =~ s/[\s;=]/_/g;

            push @strings, $label;
        }

        $line->{Extra}->{DOMAINS} = join ',', @strings if @strings;
    }

    # distance to transcript
    if($line->{Consequence} =~ /(up|down)stream/i) {
        $line->{Extra}->{DISTANCE} = $tv->distance_to_transcript;
    }

    # gene symbol
    if(defined $config->{symbol}) {
        my $symbol  = $tr->{_gene_symbol} || $tr->{_gene_hgnc};
        my $source  = $tr->{_gene_symbol_source};
        my $hgnc_id = $tr->{_gene_hgnc_id} if defined($tr->{_gene_hgnc_id});

        $line->{Extra}->{SYMBOL} = $symbol if defined($symbol) && $symbol ne '-';
        $line->{Extra}->{SYMBOL_SOURCE} = $source if defined($source) && $source ne '-';
        $line->{Extra}->{HGNC_ID} = $hgnc_id if defined($hgnc_id) && $hgnc_id ne '-';
    }

    # CCDS
    $line->{Extra}->{CCDS} = $tr->{_ccds} if
      defined($config->{ccds}) &&
      defined($tr->{_ccds}) &&
      $tr->{_ccds} ne '-';

    # refseq xref
    $line->{Extra}->{RefSeq} = $tr->{_refseq} if
      defined($config->{xref_refseq}) &&
      defined($tr->{_refseq}) &&
      $tr->{_refseq} ne '-';

    # refseq match info
    if(defined($config->{refseq}) || defined($config->{merged})) {
      my @rseq_attrs = grep {$_->code =~ /^rseq/} @attribs;
      $line->{Extra}->{REFSEQ_MATCH} = join(",", map {$_->code} @rseq_attrs) if scalar @rseq_attrs;
    }

    # protein ID
    $line->{Extra}->{ENSP} = $tr->{_protein} if
      defined($config->{protein}) &&
      defined($tr->{_protein}) &&
      $tr->{_protein} ne '-';

    # uniprot
    if(defined $config->{uniprot}) {
        for my $db(qw(swissprot trembl uniparc)) {
            my $id = $tr->{'_'.$db};
            $id = undef if defined($id) && $id eq '-';
            $line->{Extra}->{uc($db)} = $id if defined($id);
        }
    }

    # canonical transcript
    if(defined $config->{canonical}) {
        $line->{Extra}->{CANONICAL} = 'YES' if $tr->is_canonical;
    }

    # biotype
    if(defined $config->{biotype}) {
        $line->{Extra}->{BIOTYPE} = $tr->biotype;
    }

    # source cache of transcript if using --merged
    if(defined $config->{merged} && defined $tr->{_source_cache}) {
        $line->{Extra}->{SOURCE} = $tr->{_source_cache};
    }

    # transcript support level
    if(defined($config->{tsl}) && (my ($tsl) = grep {$_->code eq 'TSL'} @attribs)) {
        if($tsl->value =~ m/tsl(\d+)/) {
            $line->{Extra}->{TSL} = $1 if $1;
        }
    }

    # APPRIS
    if(defined($config->{appris}) && (my ($appris) = grep {$_->code eq 'appris'} @attribs)) {
        if(my $value = $appris->value) {
            $value =~ s/principal/P/;
            $value =~ s/alternative/A/;
            $line->{Extra}->{APPRIS} = $value;
        }
    }

    # gene phenotype
    if(defined($config->{gene_phenotype}) && $tr->{_gene_phenotype}) {
        $line->{Extra}->{GENE_PHENO} = 1;
    }

    return $line;
}

# initialize a line hash
sub init_line {
    my $config = shift;
    my $vf = shift;
    my $base_line = shift;

    my $line = {
        Uploaded_variation  => $vf->variation_name,
        Location            => ($vf->{chr} || $vf->seq_region_name).':'.format_coords($vf->{start}, $vf->{end}),
        Existing_variation  => defined $vf->{existing} && scalar @{$vf->{existing}} ? join ",", map {$_->{variation_name} || ''} @{$vf->{existing}} : '-',
        Extra               => {},
    };

    # add custom info
    if(defined($config->{custom}) && scalar @{$config->{custom}}) {
        # merge the custom hash with the extra hash
        my $custom = get_custom_annotation($config, $vf);

        for my $key (keys %$custom) {
            $line->{Extra}->{$key} = $custom->{$key};
        }
    }

    # individual?
    if(defined($vf->{individual})) {
      $line->{Extra}->{IND} = $vf->{individual};

      # zygosity
      if(defined($vf->{genotype})) {
        my %unique = map {$_ => 1} @{$vf->{genotype}};
        $line->{Extra}->{ZYG} = (scalar keys %unique > 1 ? 'HET' : 'HOM').(defined($vf->{hom_ref}) ? 'REF' : '');
      }
    }

    # variant class
    $line->{Extra}->{VARIANT_CLASS} = $vf->class_SO_term() if defined($config->{variant_class});

    # frequencies?
    $line->{Extra}->{FREQS} = join ",", @{$vf->{freqs}} if defined($vf->{freqs});

    # minimised?
    $line->{Extra}->{MINIMISED} = 1 if $vf->{minimised};

    # gmaf?
    if(defined($config->{gmaf}) && defined($vf->{existing}) && scalar @{$vf->{existing}}) {
        my @gmafs =
          map {$_->{minor_allele}.':'.$_->{minor_allele_freq}}
          grep {defined($_->{minor_allele}) && $_->{minor_allele_freq} =~ /\d/}
          @{$vf->{existing}};

        $line->{Extra}->{GMAF} = join ",", @gmafs if scalar @gmafs;
    }

    # existing var stuff
    if(defined($vf->{existing}) && scalar @{$vf->{existing}}) {

        # 1KG MAFs?
        if(defined($config->{maf_1kg})) {
            my @pops = qw(AFR AMR ASN EAS EUR SAS);

            foreach my $var(@{$vf->{existing}}) {
                foreach my $pop(grep {defined($var->{$_})} @pops) {
                    my $freq = $var->{$pop};
                    $freq = '-' unless defined($freq);
                    $line->{Extra}->{$pop.'_MAF'} =
                        exists($line->{Extra}->{$pop.'_MAF'}) ?
                        $line->{Extra}->{$pop.'_MAF'}.','.$freq :
                        $freq;
                }
            }
        }

        # ESP MAFs?
        if(defined($config->{maf_esp})) {
            my @pops = qw(AA EA);

            foreach my $var(@{$vf->{existing}}) {
                foreach my $pop(grep {defined($var->{$_})} @pops) {
                    my $freq = $var->{$pop};
                    $freq = '-' unless defined($freq);
                    $line->{Extra}->{$pop.'_MAF'} =
                        exists($line->{Extra}->{$pop.'_MAF'}) ?
                        $line->{Extra}->{$pop.'_MAF'}.','.$freq :
                        $freq;
                }
            }
        }

        # ExAC MAFs?
        if(defined($config->{maf_exac})) {
            my @pops = ('ExAC', map {'ExAC_'.$_} qw(Adj AFR AMR EAS FIN NFE OTH SAS));

            foreach my $var(@{$vf->{existing}}) {
                foreach my $pop(grep {defined($var->{$_})} @pops) {
                    my $freq = $var->{$pop};
                    $freq = '-' unless defined($freq);
                    $line->{Extra}->{$pop.'_MAF'} =
                        exists($line->{Extra}->{$pop.'_MAF'}) ?
                        $line->{Extra}->{$pop.'_MAF'}.','.$freq :
                        $freq;
                }
            }
        }

        # clin sig and pubmed?
        foreach my $var(@{$vf->{existing}}) {
            if(defined($var->{clin_sig}) && $var->{clin_sig}) {
                $line->{Extra}->{CLIN_SIG} =
                    exists($line->{Extra}->{CLIN_SIG}) ?
                    $line->{Extra}->{CLIN_SIG}.','.$var->{clin_sig} :
                    $var->{clin_sig};
            }

            if(defined($config->{pubmed}) && defined($var->{pubmed}) && $var->{pubmed}) {
                $line->{Extra}->{PUBMED} =
                    exists($line->{Extra}->{PUBMED}) ?
                    $line->{Extra}->{PUBMED}.','.$var->{pubmed} :
                    $var->{pubmed};
            }
        }

        # somatic?
        my @somatic = map {$_->{somatic}} @{$vf->{existing}};
        $line->{Extra}->{SOMATIC} = join(",", map {$_ ||= 0; $_} @somatic) if grep {defined($_) && $_ > 0} @somatic;

        # phenotype or disease
        my @p_or_d = map {$_->{phenotype_or_disease}} @{$vf->{existing}};
        $line->{Extra}->{PHENO} = join(",", @p_or_d) if grep {defined($_) && $_ > 0} @p_or_d;
    }

    # copy entries from base_line
    merge_hashes($line, $base_line) if defined($base_line);

    return $line;
}


# get custom annotation for a single VF
sub get_custom_annotation {
    my $config = shift;
    my $vf = shift;
    my $cache = shift;

    return $vf->{custom} if defined($vf->{custom});

    my $annotation = {};

    my $chr = $vf->{chr};

    if(!defined($cache)) {
        # spoof regions
        my $regions;
        $regions->{$chr} = [$vf->{start}.'-'.$vf->{end}];
        $cache = cache_custom_annotation($config, $regions, $chr);
    }

    foreach my $custom(@{$config->{custom}}) {

        next unless defined($cache->{$chr}->{$custom->{name}});

        my ($s, $e) = ($vf->{start}, $vf->{end});

        # adjust start for BED as it is 0-based
        $s-- if $custom->{format} eq 'bed';

        # exact type must match coords of variant exactly
        if($custom->{type} eq 'exact') {

            foreach my $feature(values %{$cache->{$chr}->{$custom->{name}}->{$s}}) {

                next unless
                    $feature->{chr}   eq $chr &&
                    $feature->{start} == $s &&
                    $feature->{end}   == $e;

                $annotation->{$custom->{name}} .= $feature->{name}.',';

                foreach my $field(@{$custom->{fields}}) {
                  $annotation->{$custom->{name}."_".$field} .= $feature->{$field}.',' if defined($feature->{$field});
                }
            }
        }

        # overlap type only needs to overlap, but we need to search the whole range
        elsif($custom->{type} eq 'overlap') {
            ($s, $e) = ($e, $s) if $e < $s;

            foreach my $pos(keys %{$cache->{$chr}->{$custom->{name}}}) {
                foreach my $feature(values %{$cache->{$chr}->{$custom->{name}}->{$pos}}) {
                    next unless
                        $feature->{chr}   eq $chr &&
                        $feature->{end}   >= $s &&
                        $feature->{start} <= $e;

                    $annotation->{$custom->{name}} .= $feature->{name}.',';
                    foreach my $field(@{$custom->{fields}}) {
                      $annotation->{$custom->{name}."_".$field} = $feature->{$field} if defined($feature->{$field});
                    }
                }
            }
        }

        # trim off trailing commas
        $annotation->{$custom->{name}} =~ s/\,$//g if defined($annotation->{$custom->{name}});
        foreach my $field(@{$custom->{fields}}) {
          $annotation->{$custom->{name}."_".$field} =~ s/\,$//g if defined($annotation->{$custom->{name}."_".$field});
        }
    }

    return $annotation;
}

# takes VFs created from input, fixes and checks various things
sub validate_vf {
    my $config = shift;
    my $vf = shift;

    # user specified chr skip list
    return 0 if defined($config->{chr}) && !$config->{chr}->{$vf->{chr}};

    # check valid chromosomes
    my $valid_chrs = {};

    if($config->{cache}) {
      my $synonyms = $config->{chromosome_synonyms} || {};

      foreach my $cache_chr(keys %{get_cache_chromosomes($config)}) {
        $valid_chrs->{$cache_chr} = 1;
        $valid_chrs->{$_} = 1 for keys %{$synonyms->{$cache_chr} || {}};
      }
    }

    # only try to modify chromosome name if it's not currently valid
    unless($valid_chrs->{$vf->{chr}}) {      
      $vf->{chr} =~ s/^chr//ig unless $vf->{chr} =~ /^chromosome$/i || $vf->{chr} =~ /^CHR\_/;
      $vf->{chr} = 'MT' if $vf->{chr} eq 'M';
    }

    # fix inputs
    $vf->{strand} ||= 1;
    $vf->{strand} = ($vf->{strand} =~ /\-/ ? "-1" : "1");

    # sanity checks
    unless($vf->{start} =~ /^\d+$/ && $vf->{end} =~ /^\d+$/) {
      warning_msg($config, "WARNING: Start ".$vf->{start}." or end ".$vf->{end}." coordinate invalid on line ".$config->{line_number});
      return 0;
    }

    # check chromosome exists
    # transform if necessary
    if(defined($config->{cache})) {

      # map to top level?  
      unless($valid_chrs->{$vf->{chr}}) {

        # slice adaptor required
        if(defined($config->{sa})) {
          $vf->{slice} ||= get_slice($config, $vf->{chr}, undef, 1);

          if($vf->{slice}) {
            my $transformed = $vf->transform('toplevel');

            # copy to VF
            if($transformed) {
              $vf->{$_} = $transformed->{$_} for keys %$transformed;
              $vf->{original_chr} = $vf->{chr};
              $vf->{chr} = $vf->{slice}->seq_region_name;
            }

            # could not transform
            else {
              warning_msg($config, "WARNING: Chromosome ".$vf->{chr}." not found in cache and could not transform to toplevel on line ".$config->{line_number});
              return 0;
            }
          }

          # no slice
          else {
            warning_msg($config, "WARNING: Could not fetch slice for chromosome ".$vf->{chr}." on line ".$config->{line_number});
            return 0;
          }
        }

        # offline, can't transform
        else {
          warning_msg($config, "WARNING: Chromosome ".$vf->{chr}." not found in cache on line ".$config->{line_number});
          return 0;
        }
      }
    }

    # structural variation?
    return validate_svf($config, $vf) if $vf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature');

    # uppercase allele string
    $vf->{allele_string} =~ tr/[a-z]/[A-Z]/;

    unless($vf->{allele_string} =~ /([ACGT-]+\/*)+/) {
      warning_msg($config, "WARNING: Invalid allele string ".$vf->{allele_string}." on line ".$config->{line_number}." or possible parsing error\n");
      return 0;
    }

    # insertion should have start = end + 1
    if($vf->{allele_string} =~ /^\-\// && $vf->{start} != $vf->{end} + 1) {
        warning_msg(
            $config,
            "WARNING: Alleles look like an insertion (".
            $vf->{allele_string}.
            ") but coordinates are not start = end + 1 (START=".
            $vf->{start}.", END=".$vf->{end}.
            ") on line ".$config->{line_number}."\n"
        );
        return 0;
    }

    # check start <= end + 1
    if($vf->{start} > $vf->{end} + 1) {
        warning_msg(
            $config,
            "WARNING: start > end+1 : (START=".$vf->{start}.
            ", END=".$vf->{end}.
            ") on line ".$config->{line_number}."\n"
        );
        return 0;
    }

    # check length of reference matches seq length spanned
    my @alleles = split '\/', $vf->{allele_string};
    my $ref_allele = shift @alleles;
    my $tmp_ref_allele = $ref_allele;
    $tmp_ref_allele =~ s/\-//g;

    if($tmp_ref_allele =~ /^[ACGT]*$/ && ($vf->{end} - $vf->{start}) + 1 != length($tmp_ref_allele)) {
       warning_msg(
           $config,
           "WARNING: Length of reference allele (".$ref_allele.
           " length ".length($tmp_ref_allele).") does not match co-ordinates ".$vf->{start}."-".$vf->{end}.
           " on line ".$config->{line_number}
       );
       return 0;
    }

    # flag as unbalanced
    foreach my $allele(@alleles) {
        $allele =~ s/\-//g;
        $vf->{indel} = 1 unless length($allele) == length($tmp_ref_allele);
    }

    # check reference allele if requested
    if(defined $config->{check_ref}) {
        my $ok = 0;
        my $slice_ref_allele;

        # insertion, therefore no ref allele to check
        if($ref_allele eq '-') {
            $ok = 1;
        }
        else {
            my $slice_ref = $vf->{slice}->sub_Slice($vf->{start}, $vf->{end}, $vf->{strand});

            if(!defined($slice_ref)) {
                warning_msg($config, "WARNING: Could not fetch sub-slice from ".$vf->{chr}.":".$vf->{start}."\-".$vf->{end}."\(".$vf->{strand}."\) on line ".$config->{line_number});
            }

            else {
                $slice_ref_allele = $slice_ref->seq;
                $ok = (uc($slice_ref_allele) eq uc($ref_allele) ? 1 : 0);
            }
        }

        if(!$ok) {
            warning_msg(
                $config,
                "WARNING: Specified reference allele $ref_allele ".
                "does not match Ensembl reference allele".
                ($slice_ref_allele ? " $slice_ref_allele" : "").
                " on line ".$config->{line_number}
            );
            return 0;
        }
    }

    return 1;
}


# validate a structural variation
sub validate_svf {
    my $config = shift;
    my $svf = shift;

    if($svf->{start} > $svf->{end} + 1) {
        warning_msg(
            $config,
            "WARNING: start > end+1 : (START=".$svf->{start}.
            ", END=".$svf->{end}.
            ") on line ".$config->{line_number}."\n"
        );
        return 0;
    }

    return 1;
}


# takes a hash of VFs and fetches consequences by pre-fetching overlapping transcripts
# from database and/or cache
sub whole_genome_fetch {
    my $config = shift;
    my $chr = shift;
    my $vf_hash = shift;

    my (%vf_done, @finished_vfs, %seen_rfs);

    my $cache_chr = get_cache_chr_name($config, $chr);

    if(defined($config->{offline}) && !-e $config->{dir}.'/'.$cache_chr) {
        debug("No cache found for chromsome $cache_chr") unless defined($config->{quiet});

        foreach my $chunk(keys %{$vf_hash->{$chr}}) {
            foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
                push @finished_vfs, @{$vf_hash->{$chr}{$chunk}{$pos}};
            }
        }

        return \@finished_vfs;
    }

    my $slice_cache = $config->{slice_cache};
    build_slice_cache($config, $config->{tr_cache}) unless defined($slice_cache->{$chr});
    build_slice_cache($config, $config->{rf_cache}) unless defined($slice_cache->{$chr});

    debug("Analyzing chromosome $chr".($cache_chr ne $chr ? " ($cache_chr)" : "")) unless defined($config->{quiet});

    # custom annotations
    whole_genome_fetch_custom($config, $vf_hash, $chr) if defined($config->{custom});

    # split up normal variations from SVs
    my ($tmp_vf_hash, @svfs);

    foreach my $chunk(keys %{$vf_hash->{$chr}}) {
        foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
            foreach my $vf(@{$vf_hash->{$chr}{$chunk}{$pos}}) {

                # copy slice while we're here
                $vf->{slice} ||= $slice_cache->{$chr};
                $vf->{slice} = $slice_cache->{$chr} if defined($vf->{slice}->{is_fake}) && defined($slice_cache->{$chr});

                if($vf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {
                    push @svfs, $vf;
                }
                else {
                    push @{$tmp_vf_hash->{$chr}{$chunk}{$pos}}, $vf;
                }
            }
        }
    }

    $vf_hash = $tmp_vf_hash;

    # transcript annotations
    whole_genome_fetch_transcript($config, $vf_hash, $chr)
        unless defined($config->{no_consequences});

    # regulatory annotations
    whole_genome_fetch_reg($config, $vf_hash, $chr)
        if defined($config->{regulatory})
        && !defined($config->{no_consequences});

    # structural variations
    @finished_vfs = @{whole_genome_fetch_sv($config, \@svfs, $chr)}
        if scalar @svfs;

    # sort results into @finished_vfs array
    foreach my $chunk(keys %{$vf_hash->{$chr}}) {
        foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {

            if(defined($config->{regulatory})) {
                foreach my $type(@REG_FEAT_TYPES) {
                    $_->{lc_rf_type($type).'_variations'} ||= {} for @{$vf_hash->{$chr}{$chunk}{$pos}};
                }
            }

            if(defined($config->{custom})) {
                $_->{custom} ||= {} for @{$vf_hash->{$chr}{$chunk}{$pos}};
            }

            $_->{transcript_variations} ||= {} for @{$vf_hash->{$chr}{$chunk}{$pos}};

            # add to final array
            push @finished_vfs, @{$vf_hash->{$chr}{$chunk}{$pos}};
        }
    }

    # sort
    @finished_vfs = sort {
      ($a->{_line_number} || 1) <=> ($b->{_line_number} || 1) ||
      $a->{start} <=> $b->{start} ||
      $a->{end} <=> $b->{end}
    } @finished_vfs;

    # clean hash
    delete $vf_hash->{$chr};

    return \@finished_vfs;
}

sub whole_genome_fetch_custom {
    my $config = shift;
    my $vf_hash = shift;
    my $chr = shift;

    return unless scalar @{$config->{custom}};

    # create regions based on VFs instead of chunks
    my $tmp_regions;

    foreach my $chunk(keys %{$vf_hash->{$chr}}) {
        foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
            foreach my $vf(@{$vf_hash->{$chr}{$chunk}{$pos}}) {
                push @{$tmp_regions->{$chr}}, ($vf->{start}-1).'-'.($vf->{end}+1);
            }
        }
    }

    return unless defined($tmp_regions->{$chr});

    # cache annotations
    my $annotation_cache = cache_custom_annotation($config, $tmp_regions, $chr);

    # count and report
    my $total_annotations = 0;
    $total_annotations += scalar keys %{$annotation_cache->{$chr}->{$_}} for keys %{$annotation_cache->{$chr}};
    debug("Retrieved $total_annotations custom annotations (", (join ", ", map {(scalar keys %{$annotation_cache->{$chr}->{$_}}).' '.$_} keys %{$annotation_cache->{$chr}}), ")") unless defined($config->{quiet});

    # compare annotations to variations in hash
    debug("Analyzing custom annotations") unless defined($config->{quiet});
    my $total = scalar keys %{$vf_hash->{$chr}};
    my $i = 0;

    foreach my $chunk(keys %{$vf_hash->{$chr}}) {
        progress($config, $i++, $total);

        foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
            foreach my $vf(@{$vf_hash->{$chr}{$chunk}{$pos}}) {

                $vf->{custom} = get_custom_annotation($config, $vf, $annotation_cache);
            }
        }
    }

    end_progress($config);
}

sub whole_genome_fetch_transcript {
    my $config = shift;
    my $vf_hash = shift;
    my $chr = shift;

    my $tr_cache = $config->{tr_cache};

    my $up_size   = $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE;
    my $down_size = $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE;

    # check we have defined regions
    return unless defined($vf_hash->{$chr}) && defined($tr_cache->{$chr});

    debug("Analyzing variants") unless defined($config->{quiet});

    my $tr_counter = 0;
    my $tr_count   = scalar @{$tr_cache->{$chr}};

    while($tr_counter < $tr_count) {

        progress($config, $tr_counter, $tr_count);

        my $tr = $tr_cache->{$chr}->[$tr_counter++];

        # do each overlapping VF
        my $s = $tr->{start} - ($tr->strand == 1 ? $up_size : $down_size);
        my $e = $tr->{end} + ($tr->strand == 1 ? $down_size : $up_size);

        # get the chunks this transcript overlaps
        my %chunks;
        $chunks{$_} = 1 for (int($s/$config->{chunk_size})..int($e/$config->{chunk_size}));
        map {delete $chunks{$_} unless defined($vf_hash->{$chr}{$_})} keys %chunks;

        # pointer to previous VF
        # used to tell plugins this is the last variant analysed in this transcript
        my $previous_vf;

        foreach my $chunk(keys %chunks) {
            foreach my $vf(
                grep {$_->{start} <= $e && $_->{end} >= $s}
                map {@{$vf_hash->{$chr}{$chunk}{$_}}}
                keys %{$vf_hash->{$chr}{$chunk}}
            ) {
                my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
                    -transcript        => $tr,
                    -variation_feature => $vf,
                    -adaptor           => $config->{tva},
                    -no_ref_check      => 1,
                    -no_transfer       => 1
                );

                # prefetching stuff here prevents doing loads at the
                # end and makes progress reporting more useful
                # $tv->_prefetch_for_vep;

                $vf->add_TranscriptVariation($tv);

                if(defined($config->{individual})) {

                    # store VF on transcript, weaken reference to avoid circularity
                    push @{$tr->{vfs}->{$vf->{individual}}}, $vf;
                    weaken($tr->{vfs}->{$vf->{individual}}->[-1]);

                    delete $previous_vf->{last_in_transcript}->{$tr->stable_id};
                    $vf->{last_in_transcript}->{$tr->stable_id} = 1;
                }

                $previous_vf = $vf;
            }
        }
    }

    end_progress($config);
}

sub whole_genome_fetch_reg {
  my $config = shift;
  my $vf_hash = shift;
  my $chr = shift;

  my $rf_cache = $config->{rf_cache};

  foreach my $type(keys %{$rf_cache->{$chr}}) {
    debug("Analyzing ".$type."s") unless defined($config->{quiet});

    my $constructor = 'Bio::EnsEMBL::Variation::'.$type.'Variation';
    my $add_method  = 'add_'.$type.'Variation';

    my $rf_counter = 0;
    my $rf_count = scalar @{$rf_cache->{$chr}->{$type}};

    while($rf_counter < $rf_count) {

      progress($config, $rf_counter, $rf_count);

      my $rf = $rf_cache->{$chr}->{$type}->[$rf_counter++];

      # do each overlapping VF
      my $s = $rf->{start};
      my $e = $rf->{end};

      # get the chunks this transcript overlaps
      my %chunks;
      $chunks{$_} = 1 for (int($s/$config->{chunk_size})..int($e/$config->{chunk_size}));
      map {delete $chunks{$_} unless defined($vf_hash->{$chr}{$_})} keys %chunks;

      foreach my $chunk(keys %chunks) {
        foreach my $vf(
          grep {$_->{start} <= $e && $_->{end} >= $s}
          map {@{$vf_hash->{$chr}{$chunk}{$_}}}
          keys %{$vf_hash->{$chr}{$chunk}}
        ) {
          $vf->$add_method(
            $constructor->new(
              -variation_feature  => $vf,
              -feature            => $rf,
              -no_ref_check       => 1,
              -no_transfer        => 1
            )
          );
        }
      }
    }

    end_progress($config);
  }
}

sub whole_genome_fetch_sv {
    my $config = shift;
    my $svfs = shift;
    my $chr = shift;

    my $tr_cache = $config->{tr_cache};
    my $rf_cache = $config->{rf_cache};

    my $up_size   = $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE;
    my $down_size = $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE;

    debug("Analyzing structural variations") unless defined($config->{quiet});

    my($i, $total) = (0, scalar @$svfs);

    my @finished_vfs;

    foreach my $svf(@$svfs) {
        progress($config, $i++, $total);

        my %done_genes = ();

        if(defined($tr_cache->{$chr})) {
            foreach my $tr(
              grep {
                overlap(
                  $_->{start} - ($_->strand == 1 ? $up_size : $down_size),
                  $_->{end} + ($_->strand == 1 ? $down_size : $up_size),
                  $svf->{start},
                  $svf->{end}
                )
              } @{$tr_cache->{$chr}}
            ) {
                my $svo = Bio::EnsEMBL::Variation::TranscriptStructuralVariation->new(
                    -transcript                   => $tr,
                    -structural_variation_feature => $svf,
                    -no_transfer                  => 1
                );

                $svf->add_TranscriptStructuralVariation($svo);
            }
        }

        $svf->{transcript_structural_variations} ||= {};

        # do regulatory features
        if(defined($config->{regulatory}) && defined($rf_cache->{$chr})) {
            foreach my $rf_type(qw/RegulatoryFeature/) {#keys %{$rf_cache->{$chr}}) {
                foreach my $rf(
                    grep {
                        overlap(
                            $svf->{start}, $svf->{end},
                            $_->{start}, $_->{end}
                        )
                    }
                    @{$rf_cache->{$chr}->{$rf_type}}
                ) {
                    my $svo = Bio::EnsEMBL::Variation::StructuralVariationOverlap->new(
                        -feature                      => $rf,
                        -structural_variation_feature => $svf,
                        -no_transfer                  => 1
                    );

                    push @{$svf->{regulation_structural_variations}->{$rf_type}}, $svo;
                }

                $svf->{regulation_structural_variations}->{$rf_type} ||= [];
            }
        }

        # sort them
        #$svf->_sort_svos;
        push @finished_vfs, $svf;
    }

    end_progress($config);

    return \@finished_vfs;
}

# retrieves transcripts given region list
sub fetch_transcripts {
    my $config = shift;
    my $regions = shift;
    my $trim_regions = shift;

    my $tr_cache = $config->{tr_cache};
    my $slice_cache = $config->{slice_cache};

    my ($count_from_mem, $count_from_db, $count_from_cache, $count_duplicates, $count_trimmed) = (0, 0, 0, 0, 0);

    my %seen_trs;

    $count_from_mem = 0;
    my $region_count = 0;
    foreach my $chr(keys %{$regions}) {
        $count_from_mem += scalar @{$tr_cache->{$chr}} if defined($tr_cache->{$chr}) && ref($tr_cache->{$chr}) eq 'ARRAY';
        $region_count += scalar @{$regions->{$chr}};
    }

    my ($counter, $gencode_skip_count, $refseq_skip_count);

    debug("Reading transcript data from cache and/or database") unless defined($config->{quiet});

    foreach my $chr(keys %{$regions}) {

        my $cache_chr = get_cache_chr_name($config, $chr);
      
        ## hack to copy HGNC IDs
        my %hgnc_ids = ();
        my %refseq_stuff = ();

        foreach my $region(sort {(split '-', $a)[0] <=> (split '-', $b)[1]} @{$regions->{$chr}}) {
            progress($config, $counter++, $region_count);

            # skip regions beyond the end of the chr
            next if defined($slice_cache->{$chr}) && (split '-', $region)[0] > $slice_cache->{$chr}->length;

            next if defined($config->{loaded_tr}->{$chr}->{$region});

            # force quiet so other methods don't mess up the progress bar
            my $quiet = $config->{quiet};
            $config->{quiet} = 1;

            # try and load cache from disk if using cache
            my $tmp_cache;
            if(defined($config->{cache})) {
                #$tmp_cache = (
                #    defined($config->{cache_tr_type}) && $config->{cache_tr_type} eq 'tabix' ?
                #    load_dumped_transcript_cache_tabix($config, $chr, $region) :
                #    load_dumped_transcript_cache($config, $chr, $region)
                #);
                $tmp_cache = load_dumped_transcript_cache($config, $cache_chr, $region);
                $tmp_cache->{$chr} = delete($tmp_cache->{$cache_chr}) if $cache_chr ne $chr && $tmp_cache->{$cache_chr};
                $count_from_cache += scalar @{$tmp_cache->{$chr}} if defined($tmp_cache->{$chr});
                $config->{loaded_tr}->{$chr}->{$region} = 1;
            }

            # no cache found on disk or not using cache
            if(!defined($tmp_cache->{$chr})) {

                unless(defined($config->{write_cache}) || defined($config->{database}) || $chr =~ /LRG/) {
                    # restore quiet status
                    $config->{quiet} = $quiet;

                    warning_msg($config, "WARNING: Could not find cache for $chr\:$region; correct assembly used?");
                    next;
                }

                # spoof temporary region hash
                my $tmp_hash;
                push @{$tmp_hash->{$chr}}, $region;

                $tmp_cache = cache_transcripts($config, $tmp_hash);

                # make it an empty arrayref that gets cached
                # so we don't get confused and reload next time round
                $tmp_cache->{$chr} ||= [];

                $count_from_db += scalar @{$tmp_cache->{$chr}};

                # dump to disk if writing to cache
                (defined($config->{tabix}) ? dump_transcript_cache_tabix($config, $tmp_cache, $chr, $region) : dump_transcript_cache($config, $tmp_cache, $chr, $region)) if defined($config->{write_cache});

                $config->{loaded_tr}->{$chr}->{$region} = 1;
            }

            # add loaded transcripts to main cache
            if(defined($tmp_cache->{$chr})) {
                TRANSCRIPT: while(my $tr = shift @{$tmp_cache->{$chr}}) {

                    # there are some transcripts in the otherfeatures DB with no stable ID!
                    next unless $tr->stable_id;

                    # track already added transcripts by dbID
                    my $dbID = $tr->dbID;
                    if($seen_trs{$dbID}) {
                        $count_duplicates++;
                        next;
                    }

                    # trim out?
                    #if(defined($trim_regions) && defined($trim_regions->{$chr})) {
                    #    my $tmp_count = scalar grep {
                    #        overlap(
                    #            (split /\-/, $_)[0], (split /\-/, $_)[1],
                    #            $tr->{start}, $tr->{end}
                    #        )
                    #    } @{$trim_regions->{$chr}};
                    #
                    #    if(!$tmp_count) {
                    #        $count_trimmed++;
                    #        next;
                    #    }
                    #}

                    # using gencode basic?
                    if(defined($config->{gencode_basic}) && !(grep {$_->{code} eq 'gencode_basic'} @{$tr->get_all_Attributes})) {
                      $gencode_skip_count++;
                      next;
                    }

                    # using all_refseq?
                    if(
                      !defined($config->{all_refseq}) &&
                      (
                        (
                          defined($config->{refseq}) &&
                          ($tr->stable_id || '') !~ /^[A-Z]{2}\_\d+/
                        ) ||
                        (
                          defined($config->{merged}) &&
                          ($tr->{_source_cache} || '') eq 'RefSeq' &&
                          ($tr->stable_id || '') !~ /^[A-Z]{2}\_\d+/
                        )
                      )
                    ) {
                      $refseq_skip_count++;
                      next;
                    }

                    ## hack to copy HGNC IDs
                    $hgnc_ids{$tr->{_gene_symbol}} = $tr->{_gene_hgnc_id} if defined($tr->{_gene_hgnc_id});

                    ## hack to copy RefSeq gene stuff
                    if(defined($config->{refseq}) || defined($config->{merged})) {
                      $refseq_stuff{$tr->{_gene}->stable_id}->{$_} ||= $tr->{$_} for qw(_gene_symbol _gene_symbol_source _gene_hgnc_id);

                      # check for transcripts with overlapping exons
                      # Argh, too slow!
                      # my @exons = @{$tr->get_all_Exons};

                      # for my $i(0..($#exons - 1)) {
                      #   my $e1 = $exons[$i];

                      #   for my $j($i+1..$#exons) {
                      #       my $e2 = $exons[$j];

                      #       if(overlap($e1->start, $e1->end, $e2->start, $e2->end)) {
                      #           warning_msg($config, "Overlapping exons found in transcript ".$tr->stable_id.", skipping");
                      #           next TRANSCRIPT;
                      #       }
                      #   }
                      # }
                    }

                    $seen_trs{$dbID} = 1;

                    push @{$tr_cache->{$chr}}, $tr;
                }

                ## hack to copy HGNC IDs and RefSeq stuff
                my %counts;
                foreach my $tr(@{$tr_cache->{$chr}}) {
                  $tr->{_gene_hgnc_id} = $hgnc_ids{$tr->{_gene_symbol}} if defined($tr->{_gene_symbol}) && defined($hgnc_ids{$tr->{_gene_symbol}});

                  if(defined($config->{refseq}) || defined($config->{merged})) {
                    $tr->{$_} ||= $refseq_stuff{$tr->{_gene}->stable_id}->{$_} for qw(_gene_symbol _gene_symbol_source _gene_hgnc_id);
                    $counts{$tr->{stable_id}}++
                  }
                }

                ## now a further hack to remove duplicates...
                if(defined($config->{refseq}) || defined($config->{merged})) {
                  my @new;
                  foreach my $tr(@{$tr_cache->{$chr}}) {
                    if($counts{$tr->{stable_id}} > 1) {
                      push @new, $tr unless $tr->{source} eq 'ensembl';
                    }
                    else {
                      push @new, $tr;
                    }
                  }

                  $tr_cache->{$chr} = \@new;
                }
            }

            $tr_cache->{$chr} ||= [];

            undef $tmp_cache;

            # restore quiet status
            $config->{quiet} = $quiet;
        }
    }

    end_progress($config);

    my $tr_count = 0;
    $tr_count += scalar @{$tr_cache->{$_}} for keys %$tr_cache;

    debug("Skipped $gencode_skip_count transcripts not in Gencode basic set") if $gencode_skip_count && !defined($config->{quiet});
    debug("Retrieved $tr_count transcripts ($count_from_mem mem, $count_from_cache cached, $count_from_db DB, $count_duplicates duplicates)") unless defined($config->{quiet});

    return $tr_count;
}

sub fetch_regfeats {
    my $config = shift;
    my $regions = shift;
    my $trim_regions = shift;

    my $rf_cache = $config->{rf_cache};
    my $slice_cache = $config->{slice_cache};

    my ($count_from_mem, $count_from_db, $count_from_cache, $count_duplicates, $count_trimmed) = (0, 0, 0, 0, 0);

    my $seen_rfs = {};

    $count_from_mem = 0;
    my $region_count = 0;

    foreach my $chr(keys %$regions) {
        if(defined($rf_cache->{$chr}) && ref($rf_cache->{$chr}) eq 'HASH') {
            $count_from_mem += scalar @{$rf_cache->{$chr}->{$_}} for keys %{$rf_cache->{$chr}};
        }
        $region_count += scalar @{$regions->{$chr}};
    }

    my $counter = 0;

    debug("Reading regulatory data from cache and/or database") unless defined($config->{quiet});

    foreach my $chr(keys %$regions) {

        my $cache_chr = get_cache_chr_name($config, $chr);

        foreach my $region(sort {(split '-', $a)[0] cmp (split '-', $b)[1]} @{$regions->{$chr}}) {
            progress($config, $counter++, $region_count);

            next if defined($config->{loaded_rf}->{$chr}->{$region});

            # skip regions beyond the end of the chr
            next if defined($slice_cache->{$chr}) && (split '-', $region)[0] > $slice_cache->{$chr}->length;

            # force quiet so other methods don't mess up the progress bar
            my $quiet = $config->{quiet};
            $config->{quiet} = 1;

            # try and load cache from disk if using cache
            my $tmp_cache;
            if(defined($config->{cache})) {
                $tmp_cache = load_dumped_reg_feat_cache($config, $chr, $region);

                $tmp_cache->{$chr} = delete($tmp_cache->{$cache_chr}) if $cache_chr ne $chr && $tmp_cache->{$cache_chr};

                #$tmp_cache =
                #    defined($config->{tabix}) ?
                #    load_dumped_reg_feat_cache_tabix($config, $chr, $region, $trim_regions) :
                #    load_dumped_reg_feat_cache($config, $chr, $region);


                if(defined($tmp_cache->{$chr})) {
                    $count_from_cache += scalar @{$tmp_cache->{$chr}->{$_}} for keys %{$tmp_cache->{$chr}};
                }

                # flag as loaded
                $config->{loaded_rf}->{$chr}->{$region} = 1;
            }

            # no cache found on disk or not using cache
            if(!defined($tmp_cache->{$chr})) {

                unless(defined($config->{write_cache}) || defined($config->{database}) || $chr =~ /LRG/) {

                    # restore quiet status
                    $config->{quiet} = $quiet;

                    warning_msg($config, "WARNING: Could not find cache for $chr\:$region");
                    next;
                }

                # spoof temporary region hash
                my $tmp_hash;
                push @{$tmp_hash->{$chr}}, $region;

                $tmp_cache = cache_reg_feats($config, $tmp_hash);

                # make it an empty arrayref that gets cached
                # so we don't get confused and reload next time round
                $tmp_cache->{$chr} ||= {};

                $count_from_db += scalar @{$tmp_cache->{$chr}->{$_}} for keys %{$tmp_cache->{$chr}};

                # dump to disk if writing to cache
                #dump_reg_feat_cache($config, $tmp_cache, $chr, $region) if defined($config->{write_cache});
                (defined($config->{tabix}) ? dump_reg_feat_cache_tabix($config, $tmp_cache, $chr, $region) : dump_reg_feat_cache($config, $tmp_cache, $chr, $region)) if defined($config->{write_cache});

                # restore deleted coord_system adaptor
                foreach my $type(keys %{$tmp_cache->{$chr}}) {
                    $_->{slice}->{coord_system}->{adaptor} = $config->{csa} for @{$tmp_cache->{$chr}->{$type}};
                }

                # flag as loaded
                $config->{loaded_rf}->{$chr}->{$region} = 1;
            }

            # add loaded reg_feats to main cache
            if(defined($tmp_cache->{$chr})) {
                foreach my $type(keys %{$tmp_cache->{$chr}}) {
                    while(my $rf = shift @{$tmp_cache->{$chr}->{$type}}) {

                        # filter on cell type
                        if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {
                            next unless grep {$rf->{cell_types}->{$_}} @{$config->{cell_type}};
                        }

                        # trim out?
                        #if(defined($trim_regions) && defined($trim_regions->{$chr})) {
                        #    my $tmp_count = scalar grep {
                        #        overlap(
                        #            (split /\-/, $_)[0], (split /\-/, $_)[1],
                        #            $rf->{start}, $rf->{end}
                        #        )
                        #    } @{$trim_regions->{$chr}};
                        #
                        #    if(!$tmp_count) {
                        #        $count_trimmed++;
                        #        next;
                        #    }
                        #}

                        # track already added reg_feats by dbID
                        my $dbID = $rf->{dbID};

                        if($seen_rfs->{$type}->{$dbID}) {
                            $count_duplicates++;
                            next;
                        }
                        $seen_rfs->{$type}->{$dbID} = 1;

                        push @{$rf_cache->{$chr}->{$type}}, $rf;
                    }
                }
            }

            undef $tmp_cache;

            # restore quiet status
            $config->{quiet} = $quiet;
        }
    }

    end_progress($config);

    my $rf_count = 0;

    foreach my $chr(keys %$rf_cache) {
        foreach my $type(keys %{$rf_cache->{$chr}}) {
            $rf_count += scalar @{$rf_cache->{$chr}->{$type}};
        }
    }

    debug("Retrieved $rf_count regulatory features ($count_from_mem mem, $count_from_cache cached, $count_from_db DB, $count_duplicates duplicates)") unless defined($config->{quiet});

    return $rf_count;
}

# gets existing VFs for a vf_hash
sub check_existing_hash {
    my $config = shift;
    my $vf_hash = shift;
    my $variation_cache;

    # we only care about non-SVs here
    my %new_hash;

    foreach my $chr(keys %{$vf_hash}) {
        foreach my $chunk(keys %{$vf_hash->{$chr}}) {
            foreach my $pos(keys %{$vf_hash->{$chr}->{$chunk}}) {
                foreach my $var(grep {$_->isa('Bio::EnsEMBL::Variation::VariationFeature')} @{$vf_hash->{$chr}->{$chunk}->{$pos}}) {
                    push @{$new_hash{$chr}->{$chunk}->{$pos}}, $var;
                }
            }
        }
    }

    $vf_hash = \%new_hash;

    debug("Checking for existing variations") unless defined($config->{quiet});

    my ($chunk_count, $counter);
    $chunk_count += scalar keys %{$vf_hash->{$_}} for keys %{$vf_hash};

    foreach my $chr(keys %{$vf_hash}) {

        my $cache_chr = get_cache_chr_name($config, $chr);

        my %loaded_regions;

        foreach my $chunk(keys %{$vf_hash->{$chr}}) {
            progress($config, $counter++, $chunk_count);

            # get the VFs for this chunk
            my ($start, $end);

            # work out start and end using chunk_size
            $start = $config->{chunk_size} * $chunk;
            $end = $config->{chunk_size} * ($chunk + 1);

            # using cache?
            if(defined($config->{cache})) {
                my $tmp_regions;
                push @{$tmp_regions->{$chr}}, $start.'-'.$end;

                my $converted_regions = convert_regions($config, $tmp_regions);

                foreach my $region(@{$converted_regions->{$chr}}) {

                    unless($loaded_regions{$region}) {
                        my $tmp_cache = load_dumped_variation_cache($config, $cache_chr, $region);
                        $tmp_cache->{$chr} = delete($tmp_cache->{$cache_chr}) if $cache_chr ne $chr && $tmp_cache->{$cache_chr};

                        # load from DB if not found in cache
                        if(!defined($tmp_cache->{$chr})) {
                            unless(defined($config->{write_cache}) || defined($config->{database})) {
                                warning_msg($config, "WARNING: Could not find variation cache for $chr\:$region");
                                next;
                            }

                            $tmp_cache->{$chr} = get_variations_in_region($config, $chr, $region);
                            dump_variation_cache($config, $tmp_cache, $chr, $region) if defined($config->{write_cache});
                        }

                        # merge tmp_cache with the main cache
                        foreach my $key(keys %{$tmp_cache->{$chr}}) {
                            $variation_cache->{$chr}->{$key} = $tmp_cache->{$chr}->{$key};
                            delete $tmp_cache->{$chr}->{$key};
                        }

                        # clear memory
                        undef $tmp_cache;

                        # record this region as fetched
                        $loaded_regions{$region} = 1;
                    }
                }
            }

            # no cache, get all variations in region from DB
            else {

                my ($min, $max);

                # we can fetch smaller region when using DB
                foreach my $pos(keys %{$vf_hash->{$chr}->{$chunk}}) {
                    foreach my $var(@{$vf_hash->{$chr}->{$chunk}->{$pos}}) {
                        foreach my $coord(qw(start end)) {
                            $min = $var->{$coord} if !defined($min) || $var->{$coord} < $min;
                            $max = $var->{$coord} if !defined($max) || $var->{$coord} > $max;
                        }
                    }
                }

                $variation_cache->{$chr} = get_variations_in_region($config, $chr, $min.'-'.$max);
            }

            # now compare retrieved vars with vf_hash
            foreach my $pos(keys %{$vf_hash->{$chr}->{$chunk}}) {
                foreach my $var(@{$vf_hash->{$chr}->{$chunk}->{$pos}}) {
                    my @found;

                    if(defined($variation_cache->{$chr})) {
                        if(my $existing_vars = $variation_cache->{$chr}->{$pos}) {
                            foreach my $existing_var(grep {$_->{failed} <= $config->{failed}} @$existing_vars) {
                                unless(is_var_novel($config, $existing_var, $var)) {
                                    push @found, $existing_var;
                                }
                            }
                        }
                    }

                    $var->{existing}   = \@found;
                    $var->{existing} ||= [];
                }
            }
        }

        delete $variation_cache->{$chr};
    }

    end_progress($config);
}

# gets existing VFs for a list of VFs
sub check_existing_tabix {
  my $config = shift;
  my $listref = shift;

  debug("Checking for existing variations") unless defined($config->{quiet});

  # we only care about non-SVs here
  my %by_chr;
  push @{$by_chr{$_->{chr}}}, $_ for grep {$_->isa('Bio::EnsEMBL::Variation::VariationFeature')} @$listref;

  my $max = 200;
  my $total = scalar @$listref;
  my $p = 0;

  foreach my $chr(keys %by_chr) {
    my $list = $by_chr{$chr};

    my $cache_chr = get_cache_chr_name($config, $chr);

    while(scalar @$list) {
      my @tmp_list = sort {$a->{start} <=> $b->{start}} splice @$list, 0, $max;
      progress($config, $p, $total);
      $p += scalar @tmp_list;

      my $region_string = join " ", map {$cache_chr.':'.($_->{start} > $_->{end} ? $_->{end}.'-'.$_->{start} : $_->{start}.'-'.$_->{end})} @tmp_list;

      my $file = get_dump_file_name($config, $cache_chr, "all", "vars");
      next unless -e $file;
      #die("ERROR: Could not read from file $file\n") unless -e $file;

      open VARS, "tabix -f $file $region_string 2>&1 |"
        or die "\nERROR: Could not open tabix pipe for $file\n";

      # convert list to hash so we can look up quickly by position
      my %hash;
      push @{$hash{$_->{start}}}, $_ for @tmp_list;

      VAR: while(<VARS>) {
        chomp;
        my $existing = parse_variation($config, $_);

        foreach my $input(@{$hash{$existing->{start}} || []}) {
          if(
            $existing->{start} == $input->{start} &&
            $existing->{failed} <= $config->{failed} &&
            !is_var_novel($config, $existing, $input)
          ) {
            push @{$input->{existing}}, $existing unless
              grep {$_->{variation_name} eq $existing->{variation_name}}
              @{$input->{existing} || []};
          }
        }
      }

      close VARS;

      $_->{existing} ||= [] for @tmp_list;
    }
  }

  end_progress($config);
}

sub check_existing_tabix_pm {
  my $config = shift;
  my $listref = shift;

  debug("Checking for existing variations") unless defined($config->{quiet});

  # we only care about non-SVs here
  my %by_chr;
  push @{$by_chr{$_->{chr}}}, $_ for grep {$_->isa('Bio::EnsEMBL::Variation::VariationFeature')} @$listref;

  my $total = scalar @$listref;
  my $p = 0;

  foreach my $chr(keys %by_chr) {

    my $cache_chr = get_cache_chr_name($config, $chr);

    my $file = get_dump_file_name($config, $cache_chr, "all", "vars");
    next unless -e $file;
    my $tabix_obj = $config->{_vf_tabix}->{$chr} ||= Tabix->new(-data => $file);
    next unless $tabix_obj;

    foreach my $vf(@{$by_chr{$chr}}) {
      progress($config, $p++, $total);

      my $iter = $tabix_obj->query($cache_chr, $vf->{start} - 1, $vf->{end} + 1);
      next unless $iter && $iter->{_};

      while(my $line = $tabix_obj->read($iter)) {
        chomp $line;
        my $existing = parse_variation($config, $line);

        if(
          $existing->{start} == $vf->{start} &&
          $existing->{failed} <= $config->{failed} &&
          !is_var_novel($config, $existing, $vf)
        ) {
          push @{$vf->{existing}}, $existing unless
            grep {$_->{variation_name} eq $existing->{variation_name}}
            @{$vf->{existing} || []};
        }
      }

      $vf->{existing} ||= [];
    }
  }

  end_progress($config);
}


# gets overlapping SVs for a vf_hash
sub check_svs_hash {
    my $config = shift;
    my $vf_hash = shift;

    debug("Checking for overlapping structural variations") unless defined($config->{quiet});

    my ($chunk_count, $counter);
    $chunk_count += scalar keys %{$vf_hash->{$_}} for keys %{$vf_hash};

    foreach my $chr(keys %{$vf_hash}) {
        foreach my $chunk(keys %{$vf_hash->{$chr}}) {

            progress($config, $counter++, $chunk_count);

            # work out start and end using chunk_size
            my ($start, $end);
            $start = $config->{chunk_size} * $chunk;
            $end = $config->{chunk_size} * ($chunk + 1);

            # check for structural variations
            if(defined($config->{sa})) {
                my $slice = $config->{sa}->fetch_by_region(undef, $chr, $start, $end);

                if(defined($slice)) {
                    my $svs = $config->{svfa}->fetch_all_by_Slice($slice);

                    foreach my $pos(keys %{$vf_hash->{$chr}->{$chunk}}) {
                        foreach my $var(@{$vf_hash->{$chr}->{$chunk}->{$pos}}) {
                            my $string = join ",",
                                map {$_->variation_name}
                                grep {$_->seq_region_start <= $var->{end} && $_->seq_region_end >= $var->{start}}
                                @$svs;

                            $var->{overlapping_svs} = $string if $string;
                        }
                    }
                }
            }
        }
    }

    end_progress($config);
}

# gets a slice from the slice adaptor
sub get_slice {
    my $config = shift;
    my $chr = shift;
    my $otherfeatures = shift;
    my $use_db = shift;
    $otherfeatures ||= '';

    $chr = get_cache_chr_name($config, $chr);

    return $config->{slice_cache}->{$chr} if defined($config->{slice_cache}) && defined($config->{slice_cache}->{$chr});

    my $slice;

    # with a FASTA DB we can just spoof slices
    if(defined($config->{fasta_db}) && !defined($use_db)) {

        my $fa_length = $config->{fasta_db}->length($chr);
        my $length = $fa_length && $fa_length > 0 ? $fa_length : 1;

        $slice = Bio::EnsEMBL::Slice->new(
          -COORD_SYSTEM      => $config->{coord_system},
          -START             => 1,
          -END               => $length,
          -SEQ_REGION_NAME   => $chr,
          -SEQ_REGION_LENGTH => $length
        );

        $slice->{is_fake} = 1;

        return $slice;
    }

    return undef unless defined($config->{sa}) && defined($chr);

    # first try to get a chromosome
    eval { $slice = $config->{$otherfeatures.'sa'}->fetch_by_region(undef, $chr); };

    $config->{slice_cache}->{$chr} ||= $slice;

    return $slice;
}




# METHODS THAT DEAL WITH "REGIONS"
##################################

# gets up/down size - this can be changed by UpDownDistance plugin
sub up_down_size {
  return (sort {$a <=> $b} ($Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE, $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE))[-1];
}

# gets regions from VF hash
sub regions_from_hash {
    my $config = shift;
    my $vf_hash = shift;

    my %include_regions;

    my $up_down_size = up_down_size();

    # if using cache we just want the regions of cache_region_size
    # since that's what we'll get from the cache (or DB if no cache found)
    if(defined($config->{cache})) {

        my $region_size = $config->{cache_region_size};

        foreach my $chr(keys %$vf_hash) {
            $include_regions{$chr} = [];
            my %temp_regions;

            foreach my $chunk(keys %{$vf_hash->{$chr}}) {
                foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
                    my @tmp = sort {$a <=> $b} map {($_->{start}, $_->{end})} @{$vf_hash->{$chr}{$chunk}{$pos}};
                    my ($s, $e) = ($tmp[0] - $up_down_size, $tmp[-1] + $up_down_size);

                    my $low = int ($s / $region_size);
                    my $high = int ($e / $region_size) + 1;

                    for my $i($low..($high - 1)) {
                        $temp_regions{(($i * $region_size) + 1).'-'.(($i + 1) * $region_size)} = 1;
                    }
                }
            }

            @{$include_regions{$chr}} = keys %temp_regions;
        }
    }

    # if no cache we don't want to fetch more than is necessary, so find the
    # minimum covered region of the variations in the hash
    else {
        foreach my $chr(keys %$vf_hash) {
            $include_regions{$chr} = [];

            foreach my $chunk(keys %{$vf_hash->{$chr}}) {
                foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
                    add_region($_->{start}, $_->{end}, $include_regions{$chr}) for @{$vf_hash->{$chr}{$chunk}{$pos}};
                }
            }
        }

        # merge regions
        merge_regions(\%include_regions, $config);
    }

    return \%include_regions;
}

# adds a region to region list, expanding existing one if overlaps
sub add_region {
    my $start = shift;
    my $end = shift;
    my $region_list = shift;

    my $up_down_size = up_down_size();

    # fix end for insertions
    $end = $start if $end < $start;

    my $added = 0;
    my $i = 0;

    while ($i < scalar @$region_list) {
        my ($region_start, $region_end) = split '-', $region_list->[$i];

        if($start <= $region_end && $end >= $region_start) {
            my $new_region_start = ($start < $end ? $start : $end) - $up_down_size;
            my $new_region_end = ($start > $end ? $start : $end) + $up_down_size;

            $new_region_start = 1 if $new_region_start < 1;

            $region_start = $new_region_start if $new_region_start < $region_start;
            $region_end = $new_region_end if $new_region_end > $region_end;

            $region_list->[$i] = $region_start.'-'.$region_end;
            $added = 1;
        }

        $i++;
    }

    unless($added) {
        my $s = $start - $up_down_size;
        $s = 1 if $s < 1;

        push @{$region_list}, $s.'-'.($end + $up_down_size);
    }
}

# merges overlapping regions from scans
sub merge_regions {
    my $include_regions = shift;
    my $config = shift;
    my $consecutive = shift;
    $consecutive ||= 0;

    # now merge overlapping regions
    foreach my $chr(keys %$include_regions) {
        my $max_index = $#{$include_regions->{$chr}};
        my (@new_regions, %skip);

        for my $i(0..$max_index) {
            next if $skip{$i};
            my ($s, $e) = split '-', $include_regions->{$chr}[$i];

            for my $j(($i+1)..$max_index) {
                next if $skip{$j};
                my ($ns, $ne) = split '-', $include_regions->{$chr}[$j];

                if($s <= ($ne + $consecutive) && $e >= ($ns - $consecutive)) {
                    $s = $ns if $ns < $s;
                    $e = $ne if $ne > $e;

                    $skip{$j} = 1;
                }
            }

            push @new_regions, $s.'-'.$e;
        }

        # replace original
        $include_regions->{$chr} = \@new_regions;

        $config->{region_count} += scalar @new_regions;
    }

    return $include_regions;
}

# converts regions as determined by scan_file to regions loadable from cache
sub convert_regions {
    my $config = shift;
    my $regions = shift;

    return undef unless defined $regions;

    my $region_size = $config->{cache_region_size};

    my %new_regions;

    foreach my $chr(keys %$regions) {
        my %temp_regions;

        foreach my $region(@{$regions->{$chr}}) {
            my ($s, $e) = split '-', $region;

            my $low = int ($s / $region_size);
            my $high = int ($e / $region_size) + 1;

            for my $i($low..($high - 1)) {
                $temp_regions{(($i * $region_size) + 1).'-'.(($i + 1) * $region_size)} = 1;
            }
        }

        @{$new_regions{$chr}} = keys %temp_regions;
    }

    return \%new_regions;
}





# CACHE METHODS
###############

# prunes a cache to get rid of features not in regions in use
sub prune_cache {
    my $config  = shift;
    my $cache   = shift;
    my $regions = shift;
    my $loaded  = shift;

    # delete no longer in use chroms
    foreach my $chr(keys %$cache) {
        unless(defined $regions->{$chr} && scalar @{$regions->{$chr}}) {
            delete $cache->{$chr};
            delete $loaded->{$chr};
        }
    }

    my $new_count = 0;

    foreach my $chr(keys %$cache) {

        # get total area spanned by regions
        my ($min, $max);
        foreach my $region(@{$regions->{$chr}}) {
            my ($s, $e) = split '-', $region;
            $min = $s if !defined($min) or $s < $min;
            $max = $e if !defined($max) or $e > $max;
        }

        # transcript cache
        if(ref($cache->{$chr}) eq 'ARRAY') {
            $cache->{$chr} = prune_min_max($cache->{$chr}, $min, $max);
            $new_count += scalar @{$cache->{$chr}};
        }
        # regfeat cache
        elsif(ref($cache->{$chr}) eq 'HASH') {
            for(keys %{$cache->{$chr}}) {
                $cache->{$chr}->{$_} = prune_min_max($cache->{$chr}->{$_}, $min, $max);
                $new_count += scalar @{$cache->{$chr}->{$_}};
            }
        }

        # update loaded regions
        my %have_regions = map {$_ => 1} @{$regions->{$chr}};

        foreach my $region(keys %{$loaded->{$chr}}) {
            delete $loaded->{$chr}->{$region} unless defined $have_regions{$region};
        }
    }

    return $new_count;
}

# does the actual pruning
sub prune_min_max {
    my $array = shift;
    my $min   = shift;
    my $max   = shift;

    # splice out features not in area spanned by min/max
    my $i = 0;
    my $f_count = scalar @{$array};
    my @new_cache;

    while($i < $f_count) {
        my $f = $array->[$i];

        $i++;

        if($max - $f->{start} > 0 && $f->{end} - $min > 0) {
            push @new_cache, $f;
        }

        # do some cleaning for transcripts
        elsif(defined $f->{translation}) {
            delete $f->{translation}->{transcript};
            delete $f->{translation};
        }
    }

    undef $array;
    return \@new_cache;
}

sub get_cache_chr_name {
  my $config = shift;
  my $chr = shift;

  my $chr_name_map = $config->{_chr_name_map} ||= {};

  if(!exists($chr_name_map->{$chr})) {
    my $mapped_name = $chr;

    my $valid = get_cache_chromosomes($config);

    unless($valid->{$chr}) {

      # try synonyms first
      my $synonyms = $config->{chromosome_synonyms} || {};

      foreach my $syn(keys %{$synonyms->{$chr} || {}}) {
        if($valid->{$syn}) {
          $mapped_name = $syn;
          last;
        }
      }

      # still haven't got it
      if($mapped_name eq $chr) {

        # try adding/removing "chr"
        if($chr =~ /^chr/i) {
          my $tmp = $chr;
          $tmp =~ s/^chr//i;

          $mapped_name = $tmp if $valid->{$tmp};
        }
        elsif($valid->{'chr'.$chr}) {
          $mapped_name = 'chr'.$chr;
        }
      }
    }

    $chr_name_map->{$chr} = $mapped_name;
  }

  return $chr_name_map->{$chr};
}

# get transcripts for slices
sub cache_transcripts {
    my $config = shift;
    my $include_regions = shift;

    my $tr_cache;
    my $i;

    debug("Caching transcripts") unless defined($config->{quiet});

    foreach my $chr(keys %$include_regions) {

        my $slice = get_slice($config, $chr, undef, 1);

        next unless defined $slice;

        # get a seq_region_Slice as for patch regions $slice won't cover the whole seq_region
        my $sr_slice = $slice->seq_region_Slice();

        # prefetch some things
        $slice->is_circular;

        # trim bumf off the slice
        delete $slice->{coord_system}->{adaptor} if defined($config->{write_cache});

        # no regions?
        if(!scalar @{$include_regions->{$chr}}) {
            my $start = 1;
            my $end = $config->{cache_region_size};

            while($start < $slice->end) {
                push @{$include_regions->{$chr}}, $start.'-'.$end;
                $start += $config->{cache_region_size};
                $end += $config->{cache_region_size};
            }
        }

        my $region_count;

        if(scalar keys %$include_regions == 1) {
            my ($chr) = keys %$include_regions;
            $region_count = scalar @{$include_regions->{$chr}};
            debug("Caching transcripts for chromosome $chr") unless defined($config->{quiet});
        }

        foreach my $region(@{$include_regions->{$chr}}) {
            progress($config, $i++, $region_count || $config->{region_count});

            my ($s, $e) = split '-', $region;

            # adjust relative to seq_region
            $s = ($s - $slice->start) + 1;
            $e = ($e - $slice->start) + 1;

            # sanity check start and end
            $s = 1 if $s < 1;
            $e = $slice->length if $e > $slice->length;

            # get sub-slice
            my $sub_slice = $slice->sub_Slice($s, $e);

            # add transcripts to the cache, via a transfer to the chrom's slice
            if(defined($sub_slice)) {

                # for some reason unless seq is called here the sequence becomes Ns later
                $sub_slice->seq;

                foreach my $gene(map {$_->transfer($sr_slice)} @{$sub_slice->get_all_Genes(undef, undef, 1)}) {
                    my $gene_stable_id = $gene->stable_id;
                    my $canonical_tr_id = $gene->{canonical_transcript_id};

                    # any phenotypes?
                    my $gene_has_phenotype = 0;

                    if($config->{pfa}) {
                        my $pfs = $config->{pfa}->fetch_all_by_Gene($gene);
                        $gene_has_phenotype = $pfs && scalar @$pfs;
                    }

                    my @trs;

                    foreach my $tr(@{$gene->get_all_Transcripts}) {
                        # there are some transcripts in the otherfeatures DB with no stable ID!!!
                        next unless $tr->stable_id;

                        # in human and mouse otherfeatures DB, there may be duplicate genes
                        # skip those from analysis refseq_human_import and refseq_mouse_import
                        # $DB::single = 1;
                        next if defined($config->{refseq}) && $config->{assembly} !~ /GRCh37/i && $tr->analysis && $tr->analysis->logic_name =~ /^refseq_[a-z]+_import$/;

                        $tr->{_gene_stable_id} = $gene_stable_id;
                        $tr->{_gene} = $gene;

                        # indicate if canonical
                        $tr->{is_canonical} = 1 if defined $canonical_tr_id and $tr->dbID eq $canonical_tr_id;

                        # indicate phenotype
                        $tr->{_gene_phenotype} = $gene_has_phenotype;

                        prefetch_transcript_data($config, $tr);

                        # strip some unnecessary data from the transcript object
                        clean_transcript($tr) if defined($config->{write_cache});

                        push @trs, $tr;
                    }

                    # sort the transcripts by translation so we can share sift/polyphen stuff
                    # between transcripts and save cache space
                    if(defined($config->{write_cache}) && (defined($config->{sift}) || defined($config->{polyphen}))) {

                        my $prev_tr;

                        # sort them by peptide seqeuence as transcripts with identical peptides
                        # will have identical SIFT/PolyPhen prediction strings
                        foreach my $tr(sort {$a->{_variation_effect_feature_cache}->{peptide} cmp $b->{_variation_effect_feature_cache}->{peptide}} grep {$_->{_variation_effect_feature_cache}->{peptide}} @trs)  {

                            if(
                                defined($prev_tr) &&
                                $prev_tr->{_variation_effect_feature_cache}->{peptide}
                                    eq $tr->{_variation_effect_feature_cache}->{peptide}
                            ) {

                                foreach my $analysis(qw(sift polyphen)) {
                                    next unless defined($config->{$analysis});
                                    $tr->{_variation_effect_feature_cache}->{protein_function_predictions}->{$analysis} = $prev_tr->{_variation_effect_feature_cache}->{protein_function_predictions}->{$analysis};
                                }
                            }

                            $prev_tr = $tr;
                        }
                    }

                    # clean the gene
                    clean_gene($gene);

                    push @{$tr_cache->{$chr}}, @trs;
                }
            }
        }
    }

    end_progress($config);

    return $tr_cache;
}

# gets rid of extra bits of info attached to the transcript that we don't need
sub clean_transcript {
    my $tr = shift;

    foreach my $key(qw(display_xref external_db external_display_name external_name external_status created_date status description edits_enabled modified_date dbentries is_current analysis transcript_mapper)) {
        delete $tr->{$key} if defined($tr->{$key});
    }

    # clean attributes
    if(defined($tr->{attributes})) {
        my @new_atts;
        my %keep = map {$_ => 1} qw(gencode_basic miRNA ncRNA cds_start_NF cds_end_NF TSL appris rseq_mrna_match rseq_mrna_nonmatch rseq_5p_mismatch rseq_cds_mismatch rseq_3p_mismatch rseq_nctran_mismatch rseq_no_comparison rseq_ens_match_wt rseq_ens_match_cds rseq_ens_no_match enst_refseq_compare);
        foreach my $att(@{$tr->{attributes}}) {
            delete $att->{description};
            push @new_atts, $att if defined($keep{$att->{code}});
        }
        $tr->{attributes} = \@new_atts;
    }

    # clean the translation
    if(defined($tr->translation)) {

        # sometimes the translation points to a different transcript?
        $tr->{translation}->{transcript} = $tr;
        weaken($tr->{translation}->{transcript});

        for my $key(qw(attributes protein_features created_date modified_date dbentries)) {
            delete $tr->translation->{$key};
        }
    }
}

# gets rid of extra bits of info attached to genes. At the moment this is almost
# everything as genes are only used for their locations when looking at
# structural variations
sub clean_gene {
    my $gene = shift;

    # delete almost everything in the gene
    map {delete $gene->{$_}}
        grep {
            $_ ne 'start' &&
            $_ ne 'end' &&
            $_ ne 'strand' &&
            $_ ne 'stable_id'
        }
    keys %{$gene};
}

# build slice cache from transcript cache
sub build_slice_cache {
    my $config = shift;
    my $tr_cache = shift;

    $config->{slice_cache} ||= {};

    foreach my $chr(keys %$tr_cache) {

        my $tmp = $tr_cache->{$chr};

        if(ref($tmp) eq 'HASH') {
          foreach my $type(keys %$tmp) {
            $config->{slice_cache}->{$chr} ||= scalar @{$tmp->{$type}} ? $tmp->{$type}->[0]->slice : &get_slice($config, $chr);
          }
        }
        else {
          $config->{slice_cache}->{$chr} ||= scalar @$tmp ? $tmp->[0]->slice : &get_slice($config, $chr);
        }

        if(!defined($config->{slice_cache}->{$chr})) {
            delete $config->{slice_cache}->{$chr}
        }

        else {
            # reattach adaptor to the coord system
            $config->{slice_cache}->{$chr}->{coord_system}->{adaptor} ||= $config->{csa};

            # log length for stats
            $config->{stats}->{chr_lengths}->{$chr} ||= $config->{slice_cache}->{$chr}->end;
        }
    }

    return $config->{slice_cache};
}

# pre-fetches per-transcript data
sub prefetch_transcript_data {
    my $config = shift;
    my $tr = shift;

    # introns
    my $introns = $tr->get_all_Introns;

    if(defined($introns)) {
        foreach my $intron(@$introns) {
            foreach my $key(qw(adaptor analysis dbID next prev seqname)) {
                delete $intron->{$key};
            }
        }
    }

    $tr->{_variation_effect_feature_cache}->{introns} ||= $introns;
    $tr->{_variation_effect_feature_cache}->{sorted_exons} ||= [sort {$a->start <=> $b->start} @{$tr->get_all_Exons}];

    # translateable_seq, mapper
    $tr->{_variation_effect_feature_cache}->{translateable_seq} ||= $tr->translateable_seq;
    $tr->{_variation_effect_feature_cache}->{mapper} ||= $tr->get_TranscriptMapper;

    # three prime UTR
    my $transferred = $tr->transfer($tr->feature_Slice());

    eval {
        $tr->{_variation_effect_feature_cache}->{three_prime_utr} = $transferred->three_prime_utr();
    };
    if($@) {
        warning_msg($config, "Problem getting 3' UTR:".$@);
    }

    # peptide
    unless ($tr->{_variation_effect_feature_cache}->{peptide}) {
        my $translation = $tr->translate;
        $tr->{_variation_effect_feature_cache}->{peptide} = $translation ? $translation->seq : undef;
    }

    # protein features
    if(defined($config->{domains}) || defined($config->{write_cache})) {
        my $pfs = $tr->translation ? $tr->translation->get_all_ProteinFeatures : [];

        # clean them to save cache space
        foreach my $pf(@$pfs) {

            # remove everything but the coord, analysis and ID fields
            foreach my $key(keys %$pf) {
                delete $pf->{$key} unless
                    $key eq 'start' ||
                    $key eq 'end' ||
                    $key eq 'analysis' ||
                    $key eq 'hseqname';
            }

            # remove everything from the analysis but the display label
            foreach my $key(keys %{$pf->{analysis}}) {
                delete $pf->{analysis}->{$key} unless $key eq '_display_label';
            }
        }

        $tr->{_variation_effect_feature_cache}->{protein_features} = $pfs;
    }

    # codon table
    unless ($tr->{_variation_effect_feature_cache}->{codon_table}) {
        # for mithocondrial dna we need to to use a different codon table
        my $attrib = $tr->slice->get_all_Attributes('codon_table')->[0];

        $tr->{_variation_effect_feature_cache}->{codon_table} = $attrib ? $attrib->value : 1;
    }

    # sift/polyphen
    if(defined($config->{pfpma}) && defined($tr->{_variation_effect_feature_cache}->{peptide})) {
        my @a = qw(sift);

        # full build wants both polyphen scores
        if(defined($config->{build})) {
          push @a, ('polyphen_humvar', 'polyphen_humdiv');
        }
        # otherwise just fetch requested
        else {
          push @a, 'polyphen_'.$config->{polyphen_analysis};
        }

        foreach my $a(@a) {
            next unless defined($config->{(split "_", $a)[0]});
            $tr->{_variation_effect_feature_cache}->{protein_function_predictions}->{$a} ||= $config->{pfpma}->fetch_by_analysis_translation_md5($a, md5_hex($tr->{_variation_effect_feature_cache}->{peptide}));
            delete $tr->{_variation_effect_feature_cache}->{protein_function_predictions}->{$a}->{adaptor};
        }
    }

    # translation seq Edits
    if(my $tl = $tr->translation) {
      $tr->{_variation_effect_feature_cache}->{seq_edits} = $tl->get_all_SeqEdits();
    }

    # gene symbol
    if(defined $config->{symbol}) {

        # get from gene cache if found already
        if(defined($tr->{_gene}->{_symbol})) {
            $tr->{_gene_symbol} = $tr->{_gene}->{_symbol};
            $tr->{_gene_symbol_source} = $tr->{_gene}->{_symbol_source};
            $tr->{_gene_hgnc_id} = $tr->{_gene}->{_hgnc_id}
        }
        else {
            $tr->{_gene_symbol} ||= undef;
            $tr->{_gene_symbol_source} ||= undef;

            if(my $xref = $tr->{_gene}->display_xref) {
                $tr->{_gene_symbol} = $xref->display_id;
                $tr->{_gene_symbol_source} = $xref->dbname;
                $tr->{_gene_hgnc_id} = $xref->primary_id if $xref->dbname eq 'HGNC';
            }

            else {
                my ($entry) = @{$tr->{_gene}->get_all_DBEntries('RefSeq_gene_name')};
                $tr->{_gene_symbol} = $entry->display_id if $entry;
            }

            # cache it on the gene object too
            $tr->{_gene}->{_symbol} = $tr->{_gene_symbol};
            $tr->{_gene}->{_symbol_source} = $tr->{_gene_symbol_source};
            $tr->{_gene}->{_hgnc_id} = $tr->{_gene_hgnc_id} if defined($tr->{_gene_hgnc_id});
        }
    }

    # CCDS
    my @entries = grep {$_->database eq 'CCDS'} @{$tr->get_all_DBEntries};
    $tr->{_ccds} = $entries[0]->display_id if scalar @entries;
    $tr->{_ccds} ||= '-';

    # refseq
    @entries = grep {$_->database eq 'RefSeq_mRNA'} @{$tr->get_all_DBEntries};
    if(scalar @entries) {
        $tr->{_refseq} = join ",", map {$_->display_id} @entries;
    }
    else {
        $tr->{_refseq} = '-';
    }

    # Uniprot
    if(defined($tr->translation)) {
        @entries = grep {$_->database eq 'Uniprot/SWISSPROT'} @{$tr->translation->get_all_DBEntries};
        if(scalar @entries) {
            $tr->{_swissprot} = join ",", map {$_->primary_id} @entries;
        }
        else {
            $tr->{_swissprot} = '-';
        }

        @entries = grep {$_->database eq 'Uniprot/SPTREMBL'} @{$tr->translation->get_all_DBEntries};
        if(scalar @entries) {
            $tr->{_trembl} = join ",", map {$_->primary_id} @entries;
        }
        else {
            $tr->{_trembl} = '-';
        }

        @entries = grep {$_->database eq 'UniParc'} @{$tr->translation->get_all_DBEntries};
        if(scalar @entries) {
            $tr->{_uniparc} = join ",", map {$_->primary_id} @entries;
        }
        else {
            $tr->{_uniparc} = '-';
        }
    }

    # protein stable ID
    $tr->{_protein} = $tr->translation ? $tr->translation->stable_id : '-';

    return $tr;
}

sub get_dump_file_name {
    my $config = shift;
    my $chr    = shift;
    my $region = shift;
    my $type   = shift;

    $type ||= 'transcript';

    if($type eq 'transcript') {
        $type = '';
    }
    else {
        $type = '_'.$type;
    }

    #my ($s, $e) = split /\-/, $region;
    #my $subdir = int($s / 1e6);
    #
    #my $dir = $config->{dir}.'/'.$chr.'/'.$subdir;

    my $dir = $config->{dir}.'/'.$chr;
    my $dump_file = $dir.'/'.$region.$type.($config->{sereal} && $type ne '_var' ? '.sereal' : '.gz');

    # make directory if it doesn't exist
    if(defined($config->{write_cache}) && !(-e $dir)) {
        mkpath($dir);
    }

    return $dump_file;
}

sub serialize_to_file {
  my ($config, $file, $obj) = @_;

  if($config->{sereal}) {
    $config->{encoder} ||= Sereal::Encoder->new({compress => 1});

    open OUT, ">".$file or die("ERROR: Could not write to dump file $file");
    print OUT $config->{encoder}->encode($obj);
    close OUT;
  }

  else {
    open my $fh, "| gzip -9 -c > ".$file or die "ERROR: Could not write to dump file $file";
    nstore_fd($obj, $fh);
    close $fh;
  }
}

sub deserialize_from_file {
  my ($config, $file) = @_;

  my $obj;

  if($config->{sereal}) {
    $config->{decoder} ||= Sereal::Decoder->new();
    open IN, $file;
    $obj = $config->{decoder}->decode(join('', <IN>));
    close IN;
  }

  else {
    open my $fh, $config->{compress}." ".$file." |" or die "ERROR: $!";
    $obj = fd_retrieve($fh);
    close $fh;
  }

  return $obj;
}

# dumps out transcript cache to file
sub dump_transcript_cache {
    my $config = shift;
    my $tr_cache = shift;
    my $chr = shift;
    my $region = shift;

    debug("Dumping cached transcript data") unless defined($config->{quiet});

    # clean the slice adaptor before storing
    clean_slice_adaptor($config);

    strip_transcript_cache($config, $tr_cache);

    $config->{reg}->disconnect_all;
    delete $config->{sa}->{dbc}->{_sql_helper};

    my $dump_file = get_dump_file_name($config, $chr, $region, 'transcript');

    debug("Writing to $dump_file") unless defined($config->{quiet});

    serialize_to_file($config, $dump_file, $tr_cache);
}

# loads in dumped transcript cache to memory
sub load_dumped_transcript_cache {
    my $config = shift;
    my $chr = shift;
    my $region = shift;

    my $dump_file = get_dump_file_name($config, $chr, $region, 'transcript');

    return undef unless -e $dump_file;

    debug("Reading cached transcript data for chromosome $chr".(defined $region ? "\:$region" : "")." from dumped file") unless defined($config->{quiet});

    my $tr_cache = deserialize_from_file($config, $dump_file);

    # reattach adaptors
    foreach my $t(@{$tr_cache->{$chr}}) {
        if(defined($t->{translation})) {
            $t->{translation}->{adaptor} = $config->{tra} if defined $config->{tra};
            $t->{translation}->{transcript} = $t;
            weaken($t->{translation}->{transcript});
        }

        $t->{slice}->{adaptor} = $config->{sa};

        $_->{slice} ||= $t->{slice} for @{$t->{_trans_exon_array}};
    }

    return $tr_cache;
}

# strips cache before writing to disk
sub strip_transcript_cache {
    my $config = shift;
    my $cache = shift;

    foreach my $chr(keys %$cache) {
        foreach my $tr(@{$cache->{$chr}}) {
            foreach my $exon(@{$tr->{_trans_exon_array}}) {
                delete $exon->{slice}->{adaptor};

                for(qw(adaptor created_date modified_date is_current version is_constitutive _seq_cache dbID slice)) {
                    delete $exon->{$_};
                }
            }

            delete $tr->{adaptor};
            delete $tr->{slice}->{adaptor};
            delete $tr->{translation}->{adaptor} if defined($tr->{translation});
        }
    }
}

# cleans slice adaptor before storing in cache
sub clean_slice_adaptor{
    my $config = shift;

    # clean some stuff off the slice adaptor
    delete $config->{sa}->{asm_exc_cache};
    $config->{sa}->{sr_name_cache} = {};
    $config->{sa}->{sr_id_cache} = {};
    delete $config->{sa}->{db}->{seq_region_cache};
    delete $config->{sa}->{db}->{name_cache};
}

# dumps cached variations to disk
sub dump_variation_cache {
  my $config = shift;
  my $v_cache = shift;
  my $chr = shift;
  my $region = shift;

  my $dump_file = get_dump_file_name($config, $chr, $region, 'var');

  open DUMP, "| gzip -9 -c > ".$dump_file or die "ERROR: Could not write to adaptor dump file $dump_file";

  # get freqs from VCFs?
  if(defined($config->{freq_vcf}) && scalar @{$config->{freq_vcf}}) {
    freqs_from_vcf($config, $v_cache->{$chr}, $chr);
  }

  foreach my $pos(keys %{$v_cache->{$chr}}) {
    foreach my $v(@{$v_cache->{$chr}->{$pos}}) {
      my @tmp = (
        $v->{variation_name},
        $v->{failed} == 0 ? '' : $v->{failed},
        $v->{somatic} == 0 ? '' : $v->{somatic},
        $v->{start},
        $v->{end} == $v->{start} ? '' : $v->{end},
        $v->{allele_string},
        $v->{strand} == 1 ? '' : $v->{strand},
        $v->{minor_allele} || '',
        defined($v->{minor_allele_freq}) && $v->{minor_allele_freq} =~ /^[0-9\.]+$/ ? sprintf("%.4f", $v->{minor_allele_freq}) : '',
        $v->{clin_sig} || '',
        $v->{phenotype_or_disease} == 0 ? '' : $v->{phenotype_or_disease},
      );

      if(have_pubmed($config) && defined($config->{pubmed})) {
        push @tmp, $config->{pubmed}->{$v->{variation_name}} || '';
      }

      if(defined($config->{freq_vcf}) && scalar @{$config->{freq_vcf}}) {
        foreach my $pop(map {@{$_->{prefixed_pops} || $_->{pops}}} @{$config->{freq_vcf}}) {
          push @tmp, $v->{$pop} || '';
        }
      }

      if(defined($config->{freqs})) {
        if($config->{freqs}->{$v->{variation_name}}) {
          push @tmp, $config->{freqs}->{$v->{variation_name}};
        }
        else {
          push @tmp, '' for @{$config->{just_file_pops}};
        }
      }

      print DUMP join(" ", @tmp);
      print DUMP "\n";
    }
  }

  close DUMP;
}

# loads dumped variation cache
sub load_dumped_variation_cache {
    my $config = shift;
    my $chr = shift;
    my $region = shift;

    my $dump_file = get_dump_file_name($config, $chr, $region, 'var');

    return undef unless -e $dump_file;

    open DUMP, $config->{compress}." ".$dump_file." |" or die "ERROR: $!";

    # initialise so empty files don't spit warnings
    my $v_cache = { $chr => {} };

    while(<DUMP>) {
      chomp;
      my $v = parse_variation($config, $_);
      push @{$v_cache->{$chr}->{$v->{start}}}, $v if $v && $v->{start};
    }

    close DUMP;

    return $v_cache;
}

sub parse_variation {
  my $config = shift;
  my $line = shift;

  my @cols = @{get_variation_columns($config)};
  my $delim = defined($config->{'cache_var_type'}) && $config->{'cache_var_type'} eq 'tabix' ? "\t" : qr/ /;
  my @data = split $delim, $line;

  # assumption fix for old cache files
  if(scalar @data > scalar @cols) {
    push @cols, ('AFR', 'AMR', 'ASN', 'EUR');
  }

  # this switcher is a bit of a hack, should fix cache generation really
  my %v = map {$cols[$_] => $data[$_] eq '.' ? undef : $data[$_]} (0..(@data > @cols ? $#cols : $#data));

  $v{$_} ||= 0 for qw(failed somatic phenotype_or_disease);
  $v{end}     ||= $v{start};
  $v{strand}  ||= 1;

  # hack for odd frequency data
  if(defined($config->{old_maf})) {
    foreach my $pop(grep {defined($v{$_})} qw(AFR AMR ASN EUR)) {
     $v{$pop} =~ s/^.+?\://;
     $v{$pop} =~ s/\,.+//g;
     $v{$pop} = 1 - $v{$pop} if $v{$pop} =~ /\d+/ && $v{$pop} > 0.5;
    }
  }

  # sanity check frequency data
  foreach my $pop(grep {defined($v{$_})} qw(AFR AMR ASN EAS EUR SAS AA EA)) {
    $v{$pop} = undef unless $v{$pop} =~ /^([ACGTN-]+\:)?(0|0\.\d+|1)$/;
  }

  return \%v;
}

# gets variation cache columns
sub get_variation_columns {
    my $config = shift;

    if(!defined($config->{cache_variation_cols})) {
        my @copy = @VAR_CACHE_COLS;
        $config->{cache_variation_cols} = \@copy;
        push @{$config->{cache_variation_cols}}, 'pubmed' if have_pubmed($config) && defined($config->{pubmed});
        push @{$config->{cache_variation_cols}}, @{$config->{freq_file_pops}} if defined($config->{freq_file_pops});
    }

    return $config->{cache_variation_cols};
}

# gets frequency data from VCF during --build
sub freqs_from_vcf {
  my $config = shift;
  my $v_cache = shift;
  my $chr = shift;

  # sort by pos and exclude somatic vars
  my @list =
    grep {!$_->{somatic}}
    sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}}
    map {@{$v_cache->{$_}}}
    keys %$v_cache;
  return unless scalar @list;

  # create region string to pass to tabix
  # seems getting the whole region is faster than getting smaller chunks
  my $region_string = $chr.':'.($list[0]->{start} - 1).'-'.($list[-1]->{end} + 1);

  my %match;

  # iterate over each VCF file in the config
  foreach my $vcf_conf(@{$config->{freq_vcf}}) {
    my $file = $vcf_conf->{file};
    next unless -e $file;

    my $prefix = $vcf_conf->{prefix} || '';
    $prefix .= '_' if $prefix && $prefix !~ /\_$/;

    open VCF, "tabix -f $file $region_string 2>&1 |"
      or die "\nERROR: Could not open tabix pipe for $file\n";

    while(<VCF>) {
      chomp;

      # use parse_vcf to process the line
      # this allows us to match coordinates to input
      # without worrying about indel weirdness
      my ($vcf_vf) = @{parse_vcf($config, $_)};

      # check if we have a match
      if($v_cache->{$vcf_vf->{start}}) {

        # make sure to ignore somatic again here
        foreach my $v(grep {!$_->{somatic}} @{$v_cache->{$vcf_vf->{start}}}) {

          # several ways to match, start with simplest for speed
          my $match = 0;

          # name matches
          if($vcf_vf->{variation_name} && $v->{variation_name} eq $vcf_vf->{variation_name}) {
            $match = 1;
            $match{name}++;
          }

          # allele string and coords match
          elsif($vcf_vf->{allele_string} && $v->{allele_string} eq $vcf_vf->{allele_string} && $v->{start} == $vcf_vf->{start} && $v->{end} == $vcf_vf->{end}) {
            $match = 1;
            $match{allele_string}++;
          }

          # coords match and VCF allele string is contained in DB allele string
          elsif($v->{start} == $vcf_vf->{start} && $v->{end} == $vcf_vf->{end}) {
            my @v_alleles = split('/', $v->{allele_string});
            my @vcf_alleles = split('/', $vcf_vf->{allele_string});

            # ref allele must match
            if($v_alleles[0] eq $vcf_alleles[0]) {

              # check if the VCF alleles exist in the allele string
              my %h = map {$_ => 1} @v_alleles;
              my @m = grep {$h{$_}} @vcf_alleles;
              if(scalar @m == scalar @vcf_alleles) {
                $match = 1;
                $match{part_allele_string}++;
              }
            }
          }

          if($match) {

            # get the alleles from the VCF entry
            my @vcf_alleles = split('/', $vcf_vf->{allele_string});
            shift @vcf_alleles;

            # have to process ExAC differently from 1KG and ESP
            if($prefix =~ /exac/i) {

              my $freq_hash = {};

              # match to parts of the INFO field from the VCF
              while(m/A([CFN])\_?([A-Za-z]+)?\=([0-9\.\,\-e]+)/g) {
                my ($pop) = grep {($2 || '') eq $_} @{$vcf_conf->{pops}};

                if(defined($pop)) {

                  if($1 eq 'N') {
                    $freq_hash->{$pop}->{N} = $3;
                  }

                  else {
                    my @split = split(',', $3);

                    $freq_hash->{$pop}->{$1}->{$vcf_alleles[$_]} = $split[$_] for 0..$#split;
                  }
                }
              }

              foreach my $pop(keys %$freq_hash) {
                my @f;

                foreach my $allele(
                  keys %{{
                    map {$_ => 1} (
                      keys %{$freq_hash->{$pop}->{C}},
                      keys %{$freq_hash->{$pop}->{F}}
                    )
                  }}
                ) {
                  unless(exists($freq_hash->{$pop}->{F}->{$allele})) {
                    $freq_hash->{$pop}->{F}->{$allele} = sprintf(
                      '%.4g',
                      $freq_hash->{$pop}->{C}->{$allele} / $freq_hash->{$pop}->{N}
                    ) if defined $freq_hash->{$pop}->{C}->{$allele} && defined $freq_hash->{$pop}->{N} && $freq_hash->{$pop}->{N} > 0;
                  }

                  if(exists($freq_hash->{$pop}->{F}->{$allele})) {
                    push @f, $allele.':'.$freq_hash->{$pop}->{F}->{$allele};
                  }
                }

                my $store_name = $prefix.$pop;
                $store_name =~ s/\_$//;

                $v->{$store_name} = join(',', @f);
              }
            }

            else {
              # match to parts of the INFO field from the VCF
              while(m/([A-Z]+)\_A([CF])\=([0-9\.\,]+)/g) {

                # check the matched population name is one we're interested in
                my ($pop) = grep {$1 eq $_} @{$vcf_conf->{pops}};

                if($pop) {
                  my @f;

                  # have count, e.g AA_AC from ESP
                  if($2 eq 'C') {

                    # get total count
                    my @c = split(',', $3);
                    my $total = 0;
                    $total += $_ for @c;

                    # ESP VCFs include REF as last allele so remove it
                    pop @c;

                    my $i = 0;
                    foreach my $c(@c) {
                      my $f = sprintf('%.4f', ($c / $total));
                      $f =~ s/\.?0+$//;
                      push @f, $vcf_alleles[$i++].':'.($f || '0');
                    }
                  }

                  # have freq, e.g. AFR_AF from 1KG
                  else {
                    my $i = 0;
                    foreach my $f(split(',', $3)) {
                      push @f, $vcf_alleles[$i++].':'.($f || 0);
                    }
                  }

                  # add freqs to object
                  $v->{$prefix.$pop} = join(',', @f);
                }
              }
            }
          }
        }
      }
    }

    close VCF;
  }
}

# caches regulatory features
sub cache_reg_feats {
    my $config = shift;
    my $include_regions = shift;

    my $rf_cache;
    my $i;

    debug("Caching regulatory features") unless defined($config->{quiet});

    foreach my $chr(keys %$include_regions) {

        my $slice = get_slice($config, $chr, undef, 1);

        next unless defined $slice;

        # get a seq_region_Slice as for patch regions $slice won't cover the whole seq_region
        my $sr_slice = $slice->seq_region_Slice();

        # prefetch some things
        $slice->is_circular;

        # no regions?
        if(!scalar @{$include_regions->{$chr}}) {
            my $start = 1;
            my $end = $config->{cache_region_size};

            while($start < $slice->end) {
                push @{$include_regions->{$chr}}, $start.'-'.$end;
                $start += $config->{cache_region_size};
                $end += $config->{cache_region_size};
            }
        }

        my $region_count;

        if(scalar keys %$include_regions == 1) {
            my ($chr) = keys %$include_regions;
            $region_count = scalar @{$include_regions->{$chr}};
            debug("Caching transcripts for chromosome $chr") unless defined($config->{quiet});
        }

        foreach my $region(@{$include_regions->{$chr}}) {
            progress($config, $i++, $region_count || $config->{region_count});

            my ($s, $e) = split '-', $region;

            # adjust relative to seq_region
            $s = ($s - $slice->start) + 1;
            $e = ($e - $slice->start) + 1;

            # sanity check start and end
            $s = 1 if $s < 1;
            $e = $slice->length if $e > $slice->length;

            # get sub-slice
            my $sub_slice = $slice->sub_Slice($s, $e);
            next unless defined($sub_slice);

            $sub_slice->{coord_system}->{adaptor} = $config->{csa};

            foreach my $type(@REG_FEAT_TYPES) {
                my $features = $config->{$type.'_adaptor'}->fetch_all_by_Slice($sub_slice);
                next unless defined($features);

                # cell types
                if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {
                    foreach my $rf(@$features) {

                        my %cl;



                        # get cell type using regulatory_activity objects
                        if($type eq 'RegulatoryFeature') {
                          %cl =
                            map {$_->[0] => $_->[1]}
                            map {$_->[0] =~ s/ /\_/g; $_}
                            map {[$_->epigenome->display_label, $_->activity]}
                            grep {!$_->_is_multicell}
                            @{$rf->regulatory_activity};
                        }

                        # get cell type by fetching regfeats that contain this MotifFeature
                        elsif($type eq 'MotifFeature') {
                          %cl =
                            map {$_->[0] => $_->[1]}
                            map {$_->[0] =~ s/ /\_/g; $_}
                            map {[$_->epigenome->display_label, $_->activity]}
                            grep {!$_->_is_multicell}
                            map {@{$_->regulatory_activity}}
                            @{$config->{'RegulatoryFeature_adaptor'}->fetch_all_by_attribute_feature($rf)};
                        }

                        $rf->{cell_types} = \%cl;
                    }
                }

                push @{$rf_cache->{$chr}->{$type}},
                    map { clean_reg_feat($_) }
                    map { $_->transfer($sr_slice) }
                    @{$features};
            }
        }

        # delete reference to slice adaptor before we write to cache
        delete $slice->{adaptor} if defined($config->{write_cache});
        delete $sr_slice->{adaptor} if defined($config->{write_cache});
    }

    end_progress($config);

    return $rf_cache;
}


# cleans reg feats for caching
sub clean_reg_feat {
    my $rf = shift;

    foreach my $key(qw/adaptor binary_string bound_start bound_end attribute_cache feature_set analysis set _regulatory_activity _regulatory_build/) {
        delete $rf->{$key};
    }

    if(defined($rf->{binding_matrix})) {
        $rf->{_variation_effect_feature_cache}->{seq} = $rf->seq;

        foreach my $key(qw/adaptor feature_type analysis dbID/) {
            delete $rf->{binding_matrix}->{$key};
        }
    }

    $rf->{feature_type} = $rf->{feature_type}->{so_name} if $rf->{feature_type};

    return $rf;
}


# dumps out reg feat cache to file
sub dump_reg_feat_cache {
    my $config = shift;
    my $rf_cache = shift;
    my $chr = shift;
    my $region = shift;

    debug("Dumping cached reg feat data for $chr:$region") unless defined($config->{quiet});

    # clean the slice adaptor before storing
    clean_slice_adaptor($config);

    $config->{reg}->disconnect_all;
    delete $config->{sa}->{dbc}->{_sql_helper};

    foreach my $chr(keys %{$rf_cache}) {
        foreach my $type(keys %{$rf_cache->{$chr}}) {
            delete $_->{slice}->{coord_system}->{adaptor} for @{$rf_cache->{$chr}->{$type}};
        }
    }

    my $dump_file = get_dump_file_name($config, $chr, $region, 'reg');

    debug("Writing to $dump_file") unless defined($config->{quiet});

    serialize_to_file($config, $dump_file, $rf_cache);
}

#sub dump_reg_feat_cache_tabix {
#    my $config = shift;
#    my $rf_cache = shift;
#    my $chr = shift;
#    my $region = shift;
#
#    debug("Dumping cached reg feat data") unless defined($config->{quiet});
#
#    # clean the slice adaptor before storing
#    clean_slice_adaptor($config);
#
#    $config->{reg}->disconnect_all;
#    delete $config->{sa}->{dbc}->{_sql_helper};
#
#    $config->{reg}->disconnect_all;
#
#    my $dump_file = get_dump_file_name($config, $chr, $region, 'reg');
#
#    debug("Writing to $dump_file") unless defined($config->{quiet});
#
#    use Storable qw(nfreeze);
#    use MIME::Base64 qw(encode_base64);
#    open NEW, "| gzip -9 -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
#
#    foreach my $type(keys %{$rf_cache->{$chr}}) {
#        foreach my $rf(sort {$a->start <=> $b->start} @{$rf_cache->{$chr}->{$type}}) {
#            print NEW join "\t", (
#                $chr,
#                $rf->start,
#                $rf->end,
#                $type,
#                encode_base64(freeze($rf), "")
#            );
#            print NEW "\n";
#        }
#    }
#    close NEW;
#}

# loads in dumped transcript cache to memory
sub load_dumped_reg_feat_cache {
    my $config = shift;
    my $chr = shift;
    my $region = shift;

    my $dump_file = get_dump_file_name($config, $chr, $region, 'reg');

    return undef unless -e $dump_file;

    debug("Reading cached reg feat data for chromosome $chr".(defined $region ? "\:$region" : "")." from dumped file") unless defined($config->{quiet});

    my $rf_cache = deserialize_from_file($config, $dump_file);

    # reattach adaptors
    $_->{slice}->{adaptor} = $config->{sa} for map {@{$rf_cache->{$chr}->{$_}}} keys %{$rf_cache->{$chr}};

    return $rf_cache;
}



#sub load_dumped_reg_feat_cache_tabix {
#    my $config = shift;
#    my $chr = shift;
#    my $region = shift;
#    my $trim_regions = shift;
#
#    my $dump_file = get_dump_file_name($config, $chr, $region, 'reg');
#
#    #print STDERR "Reading from $dump_file\n";
#
#    return undef unless -e $dump_file;
#
#    debug("Reading cached reg feat data for chromosome $chr".(defined $region ? "\:$region" : "")." from dumped file") unless defined($config->{quiet});
#
#    my $rf_cache;
#
#    use MIME::Base64 qw(decode_base64);
#    use Storable qw(thaw);
#
#    my ($s, $e) = split /\-/, $region;
#    my @regions = grep {overlap($s, $e, (split /\-/, $_))} @{$trim_regions->{$chr}};
#    my $regions = "";
#    $regions .= " $chr\:$_" for @regions;
#
#    #print STDERR "tabix $dump_file $regions |\n";
#    #open IN, "tabix $dump_file $regions |";
#    open IN, "gzip -dc $dump_file |";
#    while(<IN>) {
#        my ($chr, $start, $end, $type, $blob) = split /\t/, $_;
#        next unless grep {overlap($start, $end, (split /\-/, $_))} @regions;
#        my $rf = thaw(decode_base64($blob));
#        push @{$rf_cache->{$chr}->{$type}}, $rf;
#    }
#    close IN;
#
#    $rf_cache->{$chr}->{$_} ||= [] for @REG_FEAT_TYPES;
#
#    return $rf_cache;
#}


# get custom annotation for a region
sub cache_custom_annotation {
    my $config = shift;
    my $include_regions = shift;
    my $chr = shift;

    #$include_regions = merge_regions($include_regions, $config, 1);

    my $annotation = {};

    my $total = scalar @{$config->{custom}} * scalar @{$include_regions->{$chr}};
    my $counter = 0;

    my $max_regions_per_tabix = 1000;

    debug("Caching custom annotations") unless defined($config->{quiet});

    foreach my $custom(@{$config->{custom}}) {

        my @regions = @{$include_regions->{$chr}};

        while(scalar @regions) {
            my $got_features = 0;

            my @tmp_regions = splice @regions, 0, $max_regions_per_tabix;

            progress($config, $counter, $total);
            $counter += scalar @tmp_regions;

            # some files may have e.g. chr10 instead of 10
            for my $tmp_chr($chr, 'chr'.$chr) {

                # bigwig needs to use bigWigToWig utility
                if($custom->{format} eq 'bigwig') {
                    my @tmp_files;

                    die "\nERROR: Could not find temporary directory ".$config->{tmpdir}." - use --tmpdir [dir] to define an existing directory\n" unless -d $config->{tmpdir};

                    foreach my $region(@tmp_regions) {
                        my ($s, $e) = split '-', $region;
                        my $tmp_file = $config->{tmpdir}.'/vep_tmp_'.$$.'_'.$tmp_chr.'_'.$s.'_'.$e;
                        push @tmp_files, $tmp_file;
                        my $bigwig_file = $custom->{file};
                        my $bigwig_output = `bigWigToWig -chrom=$tmp_chr -start=$s -end=$e $bigwig_file $tmp_file 2>&1`;

                        die "\nERROR: Problem using bigwig file $bigwig_file\n$bigwig_output" if $bigwig_output;
                    }

                    # concatenate all the files together
                    my $string = join(" ", @tmp_files);
                    my $tmp_file = $config->{tmpdir}.'/vep_tmp_'.$$;
                    `cat $string > $tmp_file`;
                    open CUSTOM, $tmp_file
                        or die "\nERROR: Could not read from temporary WIG file $tmp_file\n";

                    # unlink smaller files
                    unlink($_) for @tmp_files;
                }

                # otherwise use tabix
                else {
                    # tabix can fetch multiple regions, so construct a string
                    my $region_string = join " ", map {$tmp_chr.':'.$_} @tmp_regions;

                    open CUSTOM, "tabix -f ".$custom->{file}." $region_string 2>&1 |"
                        or die "\nERROR: Could not open tabix pipe for ".$custom->{file}."\n";
                }

                # set an error flag so we don't have to check every line
                my $error_flag = 1;

                # create a hash for storing temporary params (used by bigWig)
                my %tmp_params = ();

                while(<CUSTOM>) {
                    chomp;

                    # check for errors
                    if($error_flag) {
                        die "\nERROR: Problem using annotation file ".$custom->{file}."\n$_\n" if /invalid pointer|tabix|get_intv/;
                        $error_flag = 0;
                    }

                    my @data = split "\t", $_;

                    my $feature;

                    if($custom->{format} eq 'bed') {
                        $feature = {
                            chr    => $chr,
                            start  => $data[1],
                            end    => $data[2],
                            name   => $data[3],
                        };
                    }

                    elsif($custom->{format} eq 'vcf') {
                        my $tmp_vf = parse_vcf($config, $_)->[0];

                        $feature = {
                            chr    => $chr,
                            start  => $tmp_vf->{start},
                            end    => $tmp_vf->{end},
                            name   => $tmp_vf->{variation_name} || '.',
                        };

                        foreach my $field(@{$custom->{fields}}) {
                          if(m/$field\=(.+?)(\;|\s|$)/) {
                            $feature->{$field} = $1;
                          }
                        }

                        $feature->{$_} =~ s/\;/\,/g for keys %$feature;
                    }

                    elsif($custom->{format} eq 'gff' || $custom->{format} eq 'gtf') {

                        my $name;

                        # try and get a feature name from the attributes column
                        foreach my $attrib(split /\s*\;\s*/, $data[8]) {
                            my ($key, $value) = split /\=/, $attrib;
                            $name = $value if $key eq 'ID';
                        }

                        $name ||= $data[2]."_".$data[0].":".$data[3]."-".$data[4];

                        $feature = {
                            chr   => $chr,
                            start => $data[3],
                            end   => $data[4],
                            name  => $name,
                        };
                    }

                    elsif($custom->{format} eq 'bigwig') {

                        # header line from wiggle file
                        if(/^(fixed|variable)Step|^\#bedGraph/i) {
                            my @split = split /\s+/;
                            $tmp_params{type} = shift @split;
                            $tmp_params{type} =~ s/^\#//g;

                            foreach my $pair(@split) {
                                my ($key, $value) = split '=', $pair;
                                $tmp_params{$key} = $value;
                            }

                            # default to span of 1
                            $tmp_params{span} ||= 1;
                        }

                        elsif(defined($tmp_params{type})) {
                            if($tmp_params{type} eq 'fixedStep') {
                                $feature = {
                                    chr   => $chr,
                                    start => $tmp_params{start},
                                    end   => ($tmp_params{start} + $tmp_params{span}) - 1,
                                    name  => $data[0],
                                };

                                $tmp_params{start} += $tmp_params{step};
                            }
                            elsif($tmp_params{type} eq 'variableStep') {
                                $feature = {
                                    chr   => $chr,
                                    start => $data[0],
                                    end   => ($data[0] + $tmp_params{span}) - 1,
                                    name  => $data[1]
                                };
                            }
                            elsif($tmp_params{type} eq 'bedGraph') {
                                $feature = {
                                    chr   => $chr,
                                    start => $data[1] + 1,
                                    end   => $data[2],
                                    name  => $data[3]
                                };
                            }
                        }

                        else {
                          die("ERROR: Cannot parse line from bigWigtoWig output: \n$_\n");
                        }
                    }

                    if(defined($feature)) {
                        $got_features = 1;

                        if(!defined($feature->{name}) || $feature->{name} eq '.' || $custom->{coords}) {
                            $feature->{name} = $feature->{chr}.":".$feature->{start}."-".$feature->{end};
                        }

                        # add the feature to the cache
                        $annotation->{$chr}->{$custom->{name}}->{$feature->{start}}->{$feature->{name}} = $feature;
                    }
                }
                close CUSTOM;

                # unlink temporary wig files
                unlink($config->{tmpdir}.'/vep_tmp_'.$$) if $custom->{format} eq 'bigwig';

                # no need to fetch e.g. "chr21" features if just "21" worked
                last if $got_features;
            }
        }
    }

    end_progress($config);

    return $annotation;
}

# builds a full cache for this species
sub build_full_cache {
  my $config = shift;

  mkpath($config->{dir}) if !-d $config->{dir};

  my @slices = @{$config->{sa}->fetch_all('toplevel')};
  push @slices, map {$_->alternate_slice} map {@{$_->get_all_AssemblyExceptionFeatures}} @slices;
  push @slices, @{$config->{sa}->fetch_all('lrg', undef, 1, undef, 1)} if defined($config->{lrg});

  open SYN, ">".$config->{dir}."/chr_synonyms.txt" or die "ERROR: Could not write to synonyms file\n";
  foreach my $slice(@slices) {
    print SYN $slice->seq_region_name."\t".$_->name."\n" for @{$slice->get_all_synonyms};
    delete($slice->{synonym});
  }
  close SYN;

  if(lc($config->{build}) ne 'all') {
    my @i;

    foreach my $r(split(',', $config->{build})) {
      my ($f, $t) = split('-', $r);
      push @i, ($t ? ($f..$t) : $f);
    }

    my %inc = %{{map {$_ => 1} @i}};
    @slices = grep {$inc{$_->seq_region_name}} @slices;
  }

  # now do a filtering step to catch slices with the same name
  # this happens for Y in human; we only want the longest one
  my %by_name;
  $by_name{$_->seq_region_name}++ for @slices;
  # @slices = ();
  #
  # $DB::single = 1;
  #
  # foreach my $name(keys %by_name) {
  #   my @sorted = sort {$a->length <=> $b->length} @{$by_name{$name}};
  #   push @slices, $sorted[-1];
  # }

  debug("Going to dump features from ".(scalar @slices)." regions") unless defined($config->{quiet});

  # check and load pubmed
  $config->{pubmed} = get_pubmed($config) if have_pubmed($config);

  # check for unfinished build status
  my %build_done;

  if(-e $config->{dir}.'/.build') {
    open IN, $config->{dir}.'/.build';
    while(<IN>) { chomp; $build_done{$_} = 1; }
    close IN;
  }

  # open build status file
  open DONE, '>> '.$config->{dir}.'/.build' or die("ERROR: Unable to open build status file\n");

  # given build range?
  my ($range_from, $range_to) = (0, 1e12);
  if($config->{build_range}) {
    ($range_from, $range_to) = split /\D/, $config->{build_range};
  }

  foreach my $slice(@slices) {
    my $chr = $slice->seq_region_name;

    # check for features, we don't want a load of effectively empty dirs
    my $dbc = $config->{sa}->db->dbc;
    my $sth = $dbc->prepare("SELECT COUNT(*) FROM transcript WHERE seq_region_id = ?");
    $sth->execute($slice->get_seq_region_id);

    my $count;
    $sth->bind_columns(\$count);
    $sth->fetch;
    $sth->finish;

    # now check var feats
    if($config->{vfa} && $config->{vfa}->db) {
      $sth = $config->{vfa}->db->dbc->prepare("SELECT COUNT(*) FROM variation_feature WHERE seq_region_id = ?");
      $sth->execute($slice->get_seq_region_id);

      my $v_count;
      $sth->bind_columns(\$v_count);
      $sth->fetch;
      $sth->finish;

      $count += $v_count;
    }

    next unless $count > 0;

    my $regions;

    # for progress
    my $region_count = int($slice->end / $config->{cache_region_size}) + 1;
    my $counter = 0;

    # initial region
    my $start = 1 + ($config->{cache_region_size} * int($slice->start / $config->{cache_region_size}));
    my $end   = ($start - 1) + $config->{cache_region_size};

    debug((defined($config->{rebuild}) ? "Rebuild" : "Creat")."ing cache for chromosome $chr (".$slice->name.")") unless defined($config->{quiet});

    # cache slice
    $config->{slice_cache}->{$chr} = $slice;

    while($start < $slice->end) {

      if($end < $range_from || $start > $range_to) {
        $start += $config->{cache_region_size};
        $end += $config->{cache_region_size};
        next;
      }

      progress($config, $counter++, $region_count);

      # store quiet status
      my $quiet = $config->{quiet};
      $config->{quiet} = 1;

      # spoof regions
      $regions->{$chr} = [$start.'-'.$end];

      # store transcripts
      if($config->{build_parts} =~ /t/) {

        my $file = get_dump_file_name($config, $chr, $start.'-'.$end, 'transcript');

        if(!exists($build_done{$file}) || $by_name{$chr} > 1) {

          my $tmp_cache;

          # load and merge if we've already done this one
          if(exists($build_done{$file})) {
            my $old = load_dumped_transcript_cache($config, $chr, $start.'-'.$end);
            my $new = cache_transcripts($config, $regions);

            # merge
            my %merged = map {$_->dbID => $_} (@{$old->{$chr}}, @{$new->{$chr}});

            $tmp_cache->{$chr} = [sort {$a->start <=> $b->start || $b->end <=> $b->end} values %merged];
          }

          else {
            $tmp_cache = (defined($config->{rebuild}) ? load_dumped_transcript_cache($config, $chr, $start.'-'.$end) : cache_transcripts($config, $regions));
          }

          $tmp_cache->{$chr} ||= [];

          #(defined($config->{tabix}) ? dump_transcript_cache_tabix($config, $tmp_cache, $chr, $start.'-'.$end) : dump_transcript_cache($config, $tmp_cache, $chr, $start.'-'.$end));
          dump_transcript_cache($config, $tmp_cache, $chr, $start.'-'.$end);
          undef $tmp_cache;

          # restore slice adaptor
          $slice->{adaptor} ||= $config->{sa};

          print DONE "$file\n";
          $build_done{$file} = 1;
        }
      }


      # store reg feats
      if($config->{build_parts} =~ /r/ && defined($config->{regulatory})) {

        my $file = get_dump_file_name($config, $chr, $start.'-'.$end, 'reg');

        if(!exists($build_done{$file}) || $by_name{$chr} > 1) {

          my $rf_cache;

          # merge?
          if(exists($build_done{$file})) {
            my $old = load_dumped_reg_feat_cache($config, $chr, $start.'-'.$end);
            my $new = cache_reg_feats($config, $regions);

            # merge each RF type
            foreach my $type(keys %{{map {$_ => 1} (keys %{$old->{$chr}}, keys %{$new->{$chr}})}}) {

              # delete adaptors since they get re-added by load_dumped_reg_feat_cache
              if($old->{$chr}->{$type}->[0]) {
                delete $old->{$chr}->{$type}->[0]->{slice}->{adaptor};
                delete $old->{$chr}->{$type}->[0]->{slice}->{coord_system}->{adaptor};
              }

              my %merged = map {$_->dbID => $_} (@{$old->{$chr}->{$type} || []}, @{$new->{$chr}->{$type} || []});

              $rf_cache->{$chr}->{$type} = [sort {$a->start <=> $b->start || $b->end <=> $b->end} values %merged];
            }
          }
          else {
            $rf_cache = cache_reg_feats($config, $regions);
          }
          $rf_cache->{$chr} ||= {};

          delete $slice->{adaptor};
          delete $slice->{coord_system}->{adaptor};

          dump_reg_feat_cache($config, $rf_cache, $chr, $start.'-'.$end);
          #(defined($config->{tabix}) ? dump_reg_feat_cache_tabix($config, $rf_cache, $chr, $start.'-'.$end) : dump_reg_feat_cache($config, $rf_cache, $chr, $start.'-'.$end));
          undef $rf_cache;

          # this gets cleaned off but needs to be there for the next loop
          $slice->{coord_system}->{adaptor} = $config->{csa};

          # restore slice adaptor
          $slice->{adaptor} ||= $config->{sa};

          print DONE "$file\n";
          $build_done{$file} = 1;
        }
      }

      # store variations
      if($config->{build_parts} =~ /v/) {

        my $file = get_dump_file_name($config, $chr, $start.'-'.$end, 'var');

        if(!exists($build_done{$file}) || $by_name{$chr} > 1) {

          my $variation_cache;

          # merge?
          if(exists($build_done{$file})) {
            my $old = load_dumped_variation_cache($config, $chr, $start.'-'.$end);
            my $new = get_variations_in_region($config, $chr, $start.'-'.$end);

            foreach my $pos(keys %{{map {$_ => 1} (keys %{$old->{$chr}}, keys %{$new})}}) {
              my %merged = map {join("", sort grep {defined($_)} values %$_) => $_} (@{$old->{$chr}->{$pos} || []}, @{$new->{$pos} || []});
              $variation_cache->{$chr}->{$pos} = [sort {$a->{start} <=> $b->{start} || $a->{start} <=> $b->{start}} values %merged];
            }
          }
          else {
            $variation_cache->{$chr} = get_variations_in_region($config, $chr, $start.'-'.$end);
          }

          $variation_cache->{$chr} ||= {};

          dump_variation_cache($config, $variation_cache, $chr, $start.'-'.$end);
          undef $variation_cache;

          print DONE "$file\n";
          $build_done{$file} = 1;
        }
      }

      # restore quiet status
      $config->{quiet} = $quiet;

      # increment by cache_region_size to get next region
      $start += $config->{cache_region_size};
      $end += $config->{cache_region_size};
    }

    end_progress($config);

    undef $regions;
  }

  # remove build status file
  close DONE;
  unlink($config->{dir}.'/.build') or die("ERROR: Unable to remove build status file\n");

  write_cache_info($config);
}

# write an info file that defines what is in the cache
sub write_cache_info {
  my $config = shift;

  my $info_file = $config->{dir}.'/info.txt';

  open OUT, ">>$info_file" or die "ERROR: Could not write to cache info file $info_file\n";

  print OUT "# CACHE UPDATED ".get_time()."\n";

  foreach my $param(qw(
    species
    assembly
    host
    port
    user
    build
    regulatory
    sift
    polyphen
  )) {
    print OUT "$param\t".(defined $config->{$param} ? $config->{$param} : '-')."\n";
  }

  # version data
  my $version_data = get_version_data($config);
  print OUT "source\_$_\t".$version_data->{$_}."\n" for keys %$version_data;

  # cell types
  if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {

    my $regulatory_build_adaptor = $config->{RegulatoryFeature_adaptor}->db->get_RegulatoryBuildAdaptor();
    my $regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;
    my @cell_types = 
      sort
      map {s/ /\_/g; $_}
      map {$_->display_label}
      @{$regulatory_build->get_all_Epigenomes};

    print OUT "cell_types\t".(join ",", @cell_types);
    print OUT "\n";
  }

  # variation columns
  print OUT "variation_cols\t".(join ",", @{get_variation_columns($config)});
  print OUT "\n";

  close OUT;
}

# reads in cache info file
sub read_cache_info {
  my $config = shift;

  my $info_file = $config->{dir}.'/info.txt';

  open IN, $info_file or return 0;

  while(<IN>) {
    next if /^#/;
    chomp;
    my ($param, $value) = split "\t";

    if($param =~ s/^source_//) {
      $config->{version_data}->{$param} = $value;
    }
    elsif($param =~ /variation_col/) {
      $config->{'cache_'.$param} = [split ',', $value];
    }
    else {
      $config->{'cache_'.$param} = $value unless defined $value && $value eq '-';
    }
  }

  close IN;

  $config->{version_data} ||= {};

  return 1;
}

# gets list of chromosomes in cache as hashref
sub get_cache_chromosomes {
  my $config = shift;

  if(!exists($config->{cache_chromosomes})) {
    my %chrs = ();

    if(opendir DIR, $config->{dir}) {
      %chrs = map {$_ => 1} grep {!/^\./ && -d $config->{dir}.'/'.$_} readdir DIR;
      closedir DIR;
    }

    $config->{cache_chromosomes} = \%chrs;
  }

  return $config->{cache_chromosomes};
}

sub get_version_data {
  my $config = shift;

  if(!exists($config->{version_data}) && !defined($config->{offline})) {
    my %version_data = ();

    # sift/polyphen versions
    my $var_mca = $config->{reg}->get_adaptor($config->{species}, 'variation', 'metacontainer');

    if($var_mca) {
      foreach my $tool(qw(sift polyphen)) {
        if(defined($config->{$tool})) {

          my $sth = $var_mca->db->dbc->prepare(qq{
            SELECT meta_value
            FROM meta
            WHERE meta_key = ?
          });
          $sth->execute($tool.'_version');

          my $version;
          $sth->bind_columns(\$version);
          $sth->fetch();
          $sth->finish();
          $version_data{$tool} = $version if defined($version);

          # my $version = $var_mca->list_value_by_key($tool.'_version');
          # $version_data{$tool} = $version->[0] if defined($version) and scalar @$version;
        }
      }
    }

    # variation source versions
    my $sa = $config->{reg}->get_adaptor($config->{species}, 'variation', 'source');

    if($sa) {
      foreach my $source_name(qw(dbSNP COSMIC ClinVar ESP HGMD-PUBLIC)) {
        my $version = $sa->get_source_version($source_name);
        $version_data{$source_name} = $version if $version;
      }
    }

    # core source versions
    my $core_mca = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');

    if($core_mca) {
      foreach my $meta_key(qw(assembly.name gencode.version genebuild.initial_release_date)) {
        my $version = $core_mca->list_value_by_key($meta_key);
        $version_data{$meta_key} = $version->[0] if defined($version) && scalar @$version;
      }
    }

    if(defined($config->{refseq})) {
      if(my $refseq_mca = $config->{reg}->get_adaptor($config->{species}, 'otherfeatures', 'metacontainer')) {
        my $sth = $refseq_mca->db->dbc->prepare(qq{
          SELECT CONCAT(db_version, ' - ', db_file) FROM analysis WHERE logic_name = 'refseq_import' 
        });
        $sth->execute;

        my $version;
        $sth->bind_columns(\$version);
        $sth->fetch;
        $sth->finish;

        $version_data{refseq} = $version if defined($version);
      }
    }

    # funcgen versions
    if(defined($config->{regulatory})) {
      my $fg_mca = $config->{reg}->get_adaptor($config->{species}, 'funcgen', 'metacontainer');

      if($fg_mca) {
        foreach my $meta_key(qw(regbuild.version)) {
          my $version = $fg_mca->list_value_by_key($meta_key);
          $version_data{$meta_key} = $version->[0] if defined($version) && scalar @$version;
        }

        # from 85 version is in regulatory_build table
        unless($version_data{'regbuild.version'}) {
          my $sth = $fg_mca->db->dbc->prepare(qq{
            SELECT version FROM regulatory_build 
          });
          $sth->execute;

          my $version;
          $sth->bind_columns(\$version);
          $sth->fetch;
          $sth->finish;

          $version_data{'regbuild.version'} = $version if defined($version);
        }
      }
    }

    # remove extraneous meta gumpf
    foreach my $key(keys %version_data) {
      my $new_key = $key;
      $new_key =~ s/\..+//;
      $version_data{$new_key} = delete $version_data{$key};
    }

    $config->{version_data} = \%version_data;
  }

  $config->{version_data} ||= {};

  return $config->{version_data};
}

# format coords for printing
sub format_coords {
    my ($start, $end) = @_;

    if(defined($start)) {
        if(defined($end)) {
            if($start > $end) {
                return $end.'-'.$start;
            }
            elsif($start == $end) {
                return $start;
            }
            else {
                return $start.'-'.$end;
            }
        }
        else {
            return $start.'-?';
        }
    }
    elsif(defined($end)) {
        return '?-'.$end;
    }
    else  {
        return '-';
    }
}




# METHODS TO FIND CO-LOCATED / EXISTING VARIATIONS
##################################################

# compare a new vf to one from the cache / DB
sub is_var_novel {
    my $config = shift;
    my $existing_var = shift;
    my $new_var = shift;

    my $is_novel = 1;

    $is_novel = 0 if $existing_var->{start} == $new_var->start && $existing_var->{end} == $new_var->end;

    if(defined($config->{check_alleles})) {
        my %existing_alleles;

        $existing_alleles{$_} = 1 for split '\/', $existing_var->{allele_string};

        my $seen_new = 0;
        foreach my $a(split '\/', ($new_var->allele_string || "")) {
            reverse_comp(\$a) if $new_var->strand ne $existing_var->{strand};
            $seen_new = 1 unless defined $existing_alleles{$a};
        }

        $is_novel = 1 if $seen_new;
    }

    return $is_novel;
}

# check frequencies of existing var against requested params
sub check_frequencies {
    my $config = shift;
    my $var = shift;

    my $var_name = $var->{variation_name};

    my $freq_pop      = $config->{freq_pop};
    my $freq_freq     = $config->{freq_freq};
    my $freq_gt_lt    = $config->{freq_gt_lt};

    my $pass = 0;
    my $checked_cache = 0;

    delete $config->{filtered_freqs};

    # if we can, check using cached frequencies as this is way quicker than
    # going to the DB
    if($freq_pop =~ /1kg|esp/i) {
        my $freq;
        my $sub_pop = uc((split '_', $freq_pop)[-1]);

        $freq = $var->{minor_allele_freq} if $sub_pop =~ /all/i;

        if(!defined($freq)) {
            $freq = (split(':', $var->{$sub_pop}))[-1] if defined($var->{$sub_pop});
        }

        if(defined($freq) && $freq =~ /\d/) {
            $pass = 1 if $freq >= $freq_freq and $freq_gt_lt eq 'gt';
            $pass = 1 if $freq <= $freq_freq and $freq_gt_lt eq 'lt';
            push @{$config->{filtered_freqs}}, $freq_pop.':'.$freq;
        }

        $checked_cache = 1;
    }

    if(defined($config->{va}) && $checked_cache == 0) {
        my $v = $config->{va}->fetch_by_name($var_name);

        my $freq_pop_name = (split '_', $freq_pop)[-1];
        $freq_pop_name = undef if $freq_pop_name =~ /1kg|hap|any/;

        foreach my $a(@{$v->get_all_Alleles}) {
            next unless defined $a->{population} || defined $a->{'_population_id'};
            next unless defined $a->frequency;
            next if $a->frequency > 0.5;

            my $pop_name = $a->population->name;

            if($freq_pop =~ /1kg/) { next unless $pop_name =~ /^1000.+(low|phase).+/i; }
            if($freq_pop =~ /hap/) { next unless $pop_name =~ /^CSHL-HAPMAP/i; }
            if($freq_pop =~ /any/) { next unless $pop_name =~ /^(CSHL-HAPMAP)|(1000.+(low|phase).+)/i; }
            if(defined $freq_pop_name) { next unless $pop_name =~ /$freq_pop_name/i; }

            $pass = 1 if $a->frequency >= $freq_freq and $freq_gt_lt eq 'gt';
            $pass = 1 if $a->frequency <= $freq_freq and $freq_gt_lt eq 'lt';

            $pop_name =~ s/\:/\_/g;
            push @{$config->{filtered_freqs}}, $pop_name.':'.$a->frequency;

            #warn "Comparing allele ", $a->allele, " ", $a->frequency, " for $var_name in population ", $a->population->name, " PASS $pass";
        }
    }

    return 0 if $config->{freq_filter} eq 'exclude' and $pass == 1;
    return 0 if $config->{freq_filter} eq 'include' and $pass == 0;
    return 1;
}

# gets all variations in a region
sub get_variations_in_region {
    my $config = shift;
    my $chr = shift;
    my $region = shift;

    my ($start, $end) = split '-', $region;

    my %variations;

    if(defined($config->{vfa}->db)) {
        my $sr_cache = $config->{seq_region_cache};

        if(!defined($sr_cache)) {
            $sr_cache = cache_seq_region_ids($config);
            $config->{seq_region_cache} = $sr_cache;
        }

        # no seq_region_id?
        return {} unless defined($sr_cache) && defined($sr_cache->{$chr});

        my $maf_cols = have_maf_cols($config) ? 'vf.minor_allele, vf.minor_allele_freq' : 'NULL, NULL';

        my $phenotype_attrib_id = phenotype_attrib_id($config) || 0;

        my $sth = $config->{vfa}->db->dbc->prepare(qq{
            SELECT vf.variation_id, vf.variation_name, IF(fv.variation_id IS NULL, 0, 1), vf.somatic, vf.seq_region_start, vf.seq_region_end, vf.allele_string, vf.seq_region_strand, $maf_cols, REPLACE(vf.clinical_significance, " ", "_"), IF(FIND_IN_SET(?, evidence_attribs) > 0, 1, 0)
            FROM variation_feature vf
            LEFT JOIN failed_variation fv ON fv.variation_id = vf.variation_id
            WHERE vf.seq_region_id = ?
            AND vf.seq_region_start >= ?
            AND vf.seq_region_start <= ?
        });

        $sth->execute($phenotype_attrib_id, $sr_cache->{$chr}, $start, $end);

        my %v;
        $v{$_} = undef for @VAR_CACHE_COLS;

        my ($var_id, %vars_by_id);
        $sth->bind_col(1, \$var_id);

        if(have_maf_cols($config)) {
            $sth->bind_col($_+2, \$v{$VAR_CACHE_COLS[$_]}) for (0..$#VAR_CACHE_COLS);
        }
        else {
            $sth->bind_col($_+2, \$v{$VAR_CACHE_COLS[$_]}) for (0..4);
        }

        while($sth->fetch) {
            my %v_copy = %v;
            $v_copy{allele_string} =~ s/\s+/\_/g;
            push @{$variations{$v{start}}}, \%v_copy;

            # store by var_id too to get stuff from variation table
            $vars_by_id{$var_id} = \%v_copy;
        }

        $sth->finish();
    }

    return \%variations;
}

sub cache_seq_region_ids {
    my $config = shift;

    my (%cache, $chr, $id);

    my $sth = $config->{vfa}->db->dbc->prepare(qq{
        SELECT seq_region_id, name FROM seq_region
    });

    $sth->execute();
    $sth->bind_columns(\$id, \$chr);
    $cache{$chr} = $id while $sth->fetch();
    $sth->finish;

    return \%cache;
}

sub have_maf_cols {
    my $config = shift;

    if(!defined($config->{have_maf_cols})) {

        if(defined($config->{vfa}) && defined($config->{vfa}->db)) {
            my $sth = $config->{vfa}->db->dbc->prepare(qq{
                DESCRIBE variation_feature
            });

            $sth->execute();
            my @cols = map {$_->[0]} @{$sth->fetchall_arrayref};
            $sth->finish();

            $config->{have_maf_cols} = 0;
            $config->{have_maf_cols} = 1 if grep {$_ eq 'minor_allele'} @cols;
        }
        else {
            $config->{have_maf_cols} = 0;
        }
    }

    return $config->{have_maf_cols};
}

sub have_pubmed {
  my $config = shift;
  return 0 if defined($config->{build_test});

  if(!defined($config->{have_pubmed})) {

    if(defined($config->{vfa}) && defined($config->{vfa}->db)) {

      my $sth = $config->{vfa}->db->dbc->prepare(qq{
        SELECT COUNT(*) FROM variation_citation
      });
      $sth->execute;

      my $count;
      $sth->bind_columns(\$count);
      $sth->fetch();
      $sth->finish();

      $config->{have_pubmed} = $count;
    }
    else {
      $config->{have_pubmed} = 0;
    }
  }

  return $config->{have_pubmed};
}

sub get_pubmed {
  my $config = shift;

  my $sth = $config->{vfa}->db->dbc->prepare(qq{
    SELECT v.name, GROUP_CONCAT(p.pmid)
    FROM variation v, variation_citation c, publication p
    WHERE v.variation_id = c.variation_id
    AND c.publication_id = p.publication_id
    AND p.pmid IS NOT NULL
    GROUP BY v.variation_id
  });
  $sth->execute;

  my ($v, $p, %pm);
  $sth->bind_columns(\$v, \$p);
  $pm{$v} = $p while $sth->fetch();
  $sth->finish();

  return \%pm;
}

sub phenotype_attrib_id {
  my $config = shift;

  if(!exists($config->{phenotype_attrib_id})) {
    my $sth = $config->{vfa}->db->dbc->prepare(qq{
      SELECT attrib_id FROM attrib WHERE value = 'Phenotype_or_Disease';
    });

    my $a;
    $sth->execute();
    $sth->bind_columns(\$a);
    $sth->fetch();
    $sth->finish();

    $config->{phenotype_attrib_id} = $a;
  }

  return $config->{phenotype_attrib_id};
}

sub lc_rf_type {
  my $type = shift;
  my $return = $type;
  $return =~ s/Feature/_feature/;
  return lc($return);
}

sub merge_hashes {
    my ($x, $y, $add) = @_;

    foreach my $k (keys %$y) {
        if (!defined($x->{$k})) {
            $x->{$k} = $y->{$k};
        } else {
            if(ref($x->{$k}) eq 'ARRAY') {
                $x->{$k} = merge_arrays($x->{$k}, $y->{$k});
            }
            elsif(ref($x->{$k}) eq 'HASH') {
                $x->{$k} = merge_hashes($x->{$k}, $y->{$k}, $add);
            }
            else {
                $x->{$k} = ($add && $x->{$k} =~ /^[0-9\.]+$/ && $y->{$k} =~ /^[0-9\.]+$/ ? $x->{$k} + $y->{$k} : $y->{$k});
            }
        }
    }
    return $x;
}

sub merge_arrays {
    my ($x, $y) = @_;

    my %tmp = map {$_ => 1} (@$x, @$y);

    return [keys %tmp];
}




# DEBUG AND STATUS METHODS
##########################

# gets time
sub get_time() {
    my @time = localtime(time());

    # increment the month (Jan = 0)
    $time[4]++;

    # add leading zeroes as required
    for my $i(0..4) {
        $time[$i] = "0".$time[$i] if $time[$i] < 10;
    }

    # put the components together in a string
    my $time =
        ($time[5] + 1900)."-".
        $time[4]."-".
        $time[3]." ".
        $time[2].":".
        $time[1].":".
        $time[0];

    return $time;
}

# prints debug output with time
sub debug {
    my $text = (@_ ? (join "", @_) : "No message");
    my $time = get_time;

    print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}

# prints warning messages to STDERR or a log file
sub warning_msg {
  my $config = shift;
  my $text = shift;

  $text = 'WARNING: '.$text unless $text =~ /^warn/i;
  $text = $text."\n" unless $text =~ /\n$/;

  if(!defined($config->{warning_fh})) {
    if($config->{warning_file} && $config->{warning_file} =~ /^stderr$/i) {
      $config->{warning_fh} = *STDERR;
    }

    else {
      $config->{warning_fh} = FileHandle->new();
      $config->{warning_file} ||= ($config->{output_file} || 'vep').'_warnings.txt';
      my $file = $config->{warning_file};
      $config->{warning_fh}->open(">".$file) or die("ERROR: Could not write to warnings file $file\n");
    }
  }

  $config->{warning_count}++;

  my $fh = $config->{warning_fh};

  print $fh $text;

  if(!defined($config->{quiet})) {
    print STDERR $config->{no_progress} ? $text : "\n$text";
  }
}

# finds out memory usage
sub memory {
    my @mem;

    open IN, "ps -o rss,vsz $$ |";
    while(<IN>) {
        next if $_ =~ /rss/i;
        chomp;
        @mem = split;
    }
    close IN;

    return \@mem;
}

sub mem_diff {
    my $config = shift;
    my $mem = memory();
    my @diffs;

    if(defined($config->{memory})) {
        for my $i(0..(scalar @{$config->{memory}} - 1)) {
            push @diffs, $mem->[$i] - $config->{memory}->[$i];
        }
    }
    else {
        @diffs = @$mem;
    }

    $config->{memory} = $mem;

    return \@diffs;
}

# update or initiate progress bar
sub progress {
    my ($config, $i, $total) = @_;

    return if defined($config->{quiet}) || defined($config->{no_progress});

    $i = $total if $i > $total;

    my $width = $config->{terminal_width} || 60;
    my $percent = int(($i/$total) * 100);
    my $numblobs = int((($i/$total) * $width) - 2);

    # this ensures we're not writing to the terminal too much
    return if(defined($config->{prev_prog})) && $numblobs.'-'.$percent eq $config->{prev_prog};
    $config->{prev_prog} = $numblobs.'-'.$percent;

    #printf("\r%s of %s", $i, $total);
    printf("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]");
}

# end progress bar
sub end_progress {
    my $config = shift;
    return if defined($config->{quiet}) || defined($config->{no_progress});
    progress($config, 1,1);
    print "\n";
    delete $config->{prev_prog};
}

1;

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::Variation::Utils::ParseLine
#
#

=head1 NAME

Bio::EnsEMBL::Variation::Utils::ParseLine - Methods used by the VCFCollection class for parsing

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::Utils::ParseLine qw(parse_line);

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Utils::ParseLine;

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



1;

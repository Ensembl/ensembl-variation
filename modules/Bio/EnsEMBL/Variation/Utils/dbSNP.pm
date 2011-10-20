package Bio::EnsEMBL::Variation::Utils::dbSNP;

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(decode_bitfield);

use constant ENCODING_VERSION => 5;

# an example string, with the fields and offsets

#  F0 F1   F2   F3 F4 F5 F6 F7 F8 F9 
#  05 0160 000a 01 05 05 12 11 01 01
#  0  2 4  6 8  10 12 14 16 18 20 22

# offsets into the string

use constant {
    F0      => 0,
    F1_1    => 2,
    F1_2    => 4,
    F2_1    => 6,
    F2_2    => 8,
    F3      => 10,
    F4      => 12,
    F5      => 14,
    F6      => 16,
    F7      => 18,
    F8      => 20,
    F9      => 22,
};

# masks to retrieve the required bits

my %encoding = (

    version     => [F0, 0b111],

    trace_archive       => [F1_1, 0b1000_0000],
    assembly_archive    => [F1_1, 0b100_0000],
    entrez_geo          => [F1_1, 0b10_0000],
    probe_db            => [F1_1, 0b1_0000],
    entrez_gene         => [F1_1, 0b1000],
    entrez_sts          => [F1_1, 0b100],
    has_structure       => [F1_1, 0b10],
    submitter_link_out  => [F1_1, 0b1],

    clinical            => [F1_2, 0b100_0000],
    precious            => [F1_2, 0b10_0000],
    provisional_tpa     => [F1_2, 0b1_0000],
    pubmed              => [F1_2, 0b1000],
    sra                 => [F1_2, 0b100],
    organism_db_link    => [F1_2, 0b10],
    mgc_clone           => [F1_2, 0b1],
    
    utr_3       => [F2_1, 0b1000_0000],
    utr_5       => [F2_1, 0b100_0000],
    acceptor_ss => [F2_1, 0b10_0000],
    donor_ss    => [F2_1, 0b1_0000],
    intron      => [F2_1, 0b1000],
    region_3    => [F2_1, 0b100],
    region_5    => [F2_1, 0b10],
    in_gene     => [F2_1, 0b1],
     
    stop_loss   => [F2_2, 0b10_0000],
    frameshift  => [F2_2, 0b1_0000],
    missense    => [F2_2, 0b1000],
    stop_gain   => [F2_2, 0b100],
    has_ref     => [F2_2, 0b10],
    has_syn     => [F2_2, 0b1],
   
    has_other_snp           => [F3, 0b1_0000],
    has_assembly_conflict   => [F3, 0b1000],
    is_assembly_specific    => [F3, 0b100],
    weight                  => [F3, 0b11],
    
    is_mutation     => [F4, 0b1000],
    is_validated    => [F4, 0b100],
    maf_all_pops    => [F4, 0b10,],
    maf_some_pops   => [F4, 0b1],
    
    marker_high_density         => [F5, 0b100],
    in_haplotype_tagging_set    => [F5, 0b10],
    genotypes_available         => [F5, 0b1],
    
    tgp_2010_production     => [F6, 0b100_0000],
    tgp_validated           => [F6, 0b10_0000],
    tgp_2010_pilot          => [F6, 0b1_0000],
    tgp_2009_pilot          => [F6, 0b1000],
    hm_phase_3_genotyped    => [F6, 0b100],
    hm_phase_2_genotyped    => [F6, 0b10],
    hm_phase_1_genotyped    => [F6, 0b1],

    has_mesh            => [F7, 0b1000_0000],
    clinical_assay      => [F7, 0b100_0000],
    has_tf              => [F7, 0b10_0000],
    lsdb                => [F7, 0b1_0000],
    dbgap_significant   => [F7, 0b1000],
    dbgap_lod_score     => [F7, 0b100],
    third_party_annot   => [F7, 0b10],
    omim                => [F7, 0b1],

    var_class   => [F8, 0b1111],

    is_suspect                  => [F9, 0b100_0000],
    is_somatic                  => [F9, 0b10_0000],
    contig_allele_not_present   => [F9, 0b1_0000],
    withdrawn                   => [F9, 0b1000],
    cluster_no_overlap          => [F9, 0b100],
    strain_specific             => [F9, 0b10],
    genotype_conflict           => [F9, 0b1],

);

my %var_class = (
    0b0001  => 'snp',
    0b0010  => 'dips',
    0b0011  => 'heterozygous',
    0b0100  => 'microsatellite',
    0b0101  => 'named',
    0b0110  => 'no_variation',
    0b0111  => 'mixed',
    0b1000  => 'multi_base',
);

sub decode_bitfield {
    my $s = shift;

    my %res;

    for my $enc (keys %encoding) {

        my ($offset, $mask) = @{ $encoding{$enc} };
        $res{$enc} = hex(substr($s, $offset, 2)) & $mask;
        
        # check that the version matches ours
        if ($enc eq 'version' && $res{$enc} != ENCODING_VERSION) {
            warn "Version field does not match the expected version";
            return undef;
        }

        # lookup the class description 
        $res{$enc} = $var_class{$res{$enc}} if $enc eq 'var_class';
        
        # get rid of anything set to 0
        delete $res{$enc} unless $res{$enc};

    }

    return \%res;
}

1;


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

    version     => [F0, [3,2,1]],

    trace_archive       => [F1_1, 8],
    assembly_archive    => [F1_1, 7],
    entrez_geo          => [F1_1, 6],
    probe_db            => [F1_1, 5],
    entrez_gene         => [F1_1, 4],
    entrez_sts          => [F1_1, 3],
    has_structure       => [F1_1, 2],
    submitter_link_out  => [F1_1, 1],

    clinical            => [F1_2, 7],
    precious            => [F1_2, 6],
    provisional_tpa     => [F1_2, 5],
    pubmed              => [F1_2, 4],
    sra                 => [F1_2, 3],
    organism_db_link    => [F1_2, 2],
    mgc_clone           => [F1_2, 1],
    
    utr_3       => [F2_1, 8],
    utr_5       => [F2_1, 7],
    acceptor_ss => [F2_1, 6],
    donor_ss    => [F2_1, 5],
    intron      => [F2_1, 4],
    region_3    => [F2_1, 3],
    region_5    => [F2_1, 2],
    in_gene     => [F2_1, 1],
     
    stop_loss   => [F2_2, 6],
    frameshift  => [F2_2, 5],
    missense    => [F2_2, 4],
    stop_gain   => [F2_2, 3],
    has_ref     => [F2_2, 2],
    has_syn     => [F2_2, 1],
   
    has_other_snp           => [F3, 5],
    has_assembly_conflict   => [F3, 4],
    is_assembly_specific    => [F3, 3],
    weight                  => [F3, [1,2]],
    
    is_mutation     => [F4, 4],
    is_validated    => [F4, 3],
    maf_all_pops    => [F4, 2],
    maf_some_pops   => [F4, 1],
    
    marker_high_density         => [F5, 3],
    in_haplotype_tagging_set    => [F5, 2],
    genotypes_available         => [F5, 1],
    
    tgp_2010_production     => [F6, 7],
    tgp_validated           => [F6, 6],
    tgp_2010_pilot          => [F6, 5],
    tgp_2009_pilot          => [F6, 4],
    hm_phase_3_genotyped    => [F6, 3],
    hm_phase_2_genotyped    => [F6, 2],
    hm_phase_1_genotyped    => [F6, 1],

    has_mesh            => [F7, 8],
    clinical_assay      => [F7, 7],
    has_tf              => [F7, 6],
    lsdb                => [F7, 5],
    dbgap_significant   => [F7, 4],
    dbgap_lod_score     => [F7, 3],
    third_party_annot   => [F7, 2],
    omim                => [F7, 1],

    var_class   => [F8, [4,3,2,1]],

    is_suspect                  => [F9, 7],
    is_somatic                  => [F9, 6],
    contig_allele_not_present   => [F9, 5],
    withdrawn                   => [F9, 4],
    cluster_no_overlap          => [F9, 3],
    strain_specific             => [F9, 2],
    genotype_conflict           => [F9, 1],

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

        my ($offset, $bits) = @{ $encoding{$enc} };

        # if bits isn't an array, put the single bit into an array
        $bits = [$bits] unless ref $bits eq 'ARRAY';

        # OR together all the bits to give us our mask
        my $mask;

        for my $bit (@$bits) {
            $mask |= 2**($bit-1);
        }

        $res{$enc} = hex(substr($s, $offset, 2)) & $mask;
        
        # check that the version matches ours
        if ($enc eq 'version' && $res{$enc} != ENCODING_VERSION) {
            print $res{$enc}, "\n";
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


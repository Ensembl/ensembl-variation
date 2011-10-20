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

# offsets into the string for each byte

my %offsets = (
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
);

# a hash mapping the values encoded in each field to the bits used encode them

my %fields= (

    F0 => {
        version => [3,2,1],
    },

    F1_1 => {
        trace_archive       => 8,
        assembly_archive    => 7,
        entrez_geo          => 6,
        probe_db            => 5,
        entrez_gene         => 4,
        entrez_sts          => 3,
        has_structure       => 2,
        submitter_link_out  => 1,
    },
   
    F1_2 => {
        clinical            => 7,
        precious            => 6,
        provisional_tpa     => 5,
        pubmed              => 4,
        sra                 => 3,
        organism_db_link    => 2,
        mgc_clone           => 1,
    },

    F2_1 => { 
        utr_3       => 8,
        utr_5       => 7,
        acceptor_ss => 6,
        donor_ss    => 5,
        intron      => 4,
        region_3    => 3,
        region_5    => 2,
        in_gene     => 1,
    },

    F2_2 => {
        stop_loss   => 6,
        frameshift  => 5,
        missense    => 4,
        stop_gain   => 3,
        has_ref     => 2,
        has_syn     => 1,
    },

    F3 => {
        has_other_snp           => 5,
        has_assembly_conflict   => 4,
        is_assembly_specific    => 3,
        weight                  => [1,2],
    },
    
    F4 => {
        is_mutation     => 4,
        is_validated    => 3,
        maf_all_pops    => 2,
        maf_some_pops   => 1,
    },
    
    F5 => {
        marker_high_density         => 3,
        in_haplotype_tagging_set    => 2,
        genotypes_available         => 1,
    },

    F6 => {
        tgp_2010_production     => 7,
        tgp_validated           => 6,
        tgp_2010_pilot          => 5,
        tgp_2009_pilot          => 4,
        hm_phase_3_genotyped    => 3,
        hm_phase_2_genotyped    => 2,
        hm_phase_1_genotyped    => 1,
    },

    F7 => {
        has_mesh            => 8,
        clinical_assay      => 7,
        has_tf              => 6,
        lsdb                => 5,
        dbgap_significant   => 4,
        dbgap_lod_score     => 3,
        third_party_annot   => 2,
        omim                => 1,
    },

    F8 => {
        var_class   => [4,3,2,1],
    },

    F9 => {
        is_suspect                  => 7,
        is_somatic                  => 6,
        contig_allele_not_present   => 5,
        withdrawn                   => 4,
        cluster_no_overlap          => 3,
        strain_specific             => 2,
        genotype_conflict           => 1,
    },
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

    for my $field (keys %fields) {

        for my $value (keys %{ $fields{$field} }) {

            my $bits = $fields{$field}->{$value};

            # if bits isn't an array, put the single bit into an array
            $bits = [$bits] unless ref $bits eq 'ARRAY';

            # OR together all the bits to give us our mask
            my $mask;

            for my $bit (@$bits) {
                $mask |= 2**($bit-1);
            }
            
            # extract the relevant characters from the string, 
            # convert them to an integer, and apply our mask
            $res{$value} = hex(substr($s, $offsets{$field}, 2)) & $mask;
        
            # check that the version matches ours
            if ($value eq 'version' && $res{$value} != ENCODING_VERSION) {
                warn "Version field does not match the expected version (".$res{$value}." vs ".ENCODING_VERSION.")";
                return undef;
            }

            # lookup the class description 
            $res{$value} = $var_class{$res{$value}} if $value eq 'var_class';
            
            # get rid of anything set to 0
            delete $res{$value} unless $res{$value};
        }
    }

    return \%res;
}

1;


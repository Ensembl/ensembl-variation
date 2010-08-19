package Bio::EnsEMBL::Variation::Utils::VariationEffect;

use strict;
use warnings;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&transcript_effect &regulatory_region_effect &get_all_variation_effects);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Data::Dumper;

sub overlap {
    my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
    return ($f1_end >= $f2_start and $f1_start <= $f2_end);
}

sub affects_transcript {
    my ($vf, $tran) = @_;
    return not overlap($vf->start, $vf->end, $tran->start - 5000, $tran->end + 5000);
}

sub within_transcript {
    my ($vf, $tran) = @_;
    return overlap($vf->start, $vf->end, $tran->start, $tran->end);
}

sub _before_start {
    my ($vf, $tran, $dist) = @_;
    #return 0 if within_transcript($vf, $tran);
    return ($vf->end >= ($tran->start - $dist) and $vf->end < $tran->start);
}

sub _after_end {
    my ($vf, $tran, $dist) = @_;
    #return 0 if within_transcript($vf, $tran);
    return ($vf->start <= ($tran->end + $dist) and $vf->start > $tran->end);
}

sub upstream {
    my ($vf, $tran, $dist) = @_;
    return $tran->strand == 1 ? _before_start($vf, $tran, $dist) : _after_end($vf, $tran, $dist);
}

sub downstream {
    my ($vf, $tran, $dist) = @_;
    return $tran->strand == 1 ? _after_end($vf, $tran, $dist) : _before_start($vf, $tran, $dist);
}

sub upstream_5KB {
    my ($vf, $tran) = @_;
    return upstream($vf, $tran, 5000);
}

sub downstream_5KB {
    my ($vf, $tran) = @_;
    return downstream($vf, $tran, 5000);
}

sub upstream_2KB {
    my ($vf, $tran) = @_;
    return upstream($vf, $tran, 2000);
}

sub downstream_2KB {
    my ($vf, $tran) = @_;
    return downstream($vf, $tran, 2000);
}

sub within_nmd_transcript {
    my ($vf, $tran) = @_;
    return (within_transcript($vf, $tran) and $tran->biotype eq 'nonsense_mediated_decay');
}

sub _start_splice_site {
    my ($vf, $tran) = @_;

    for my $intron (@{ $tran->get_all_Introns }) {
        if (overlap($vf->start, $vf->end, $intron->start, $intron->start+2)) {
            return 1;
        }
    }

    return 0;
}

sub _end_splice_site {
    my ($vf, $tran) = @_;
    
    for my $intron (@{ $tran->get_all_Introns }) {
        if (overlap($vf->start, $vf->end, $intron->end-2, $intron->end)) {
            return 1;
        }
    }

    return 0;
}

sub donor_splice_site {
    my ($vf, $tran) = @_;
    return $tran->strand == 1 ? _start_splice_site($vf, $tran) : _end_splice_site($vf, $tran);
}

sub acceptor_splice_site {
    my ($vf, $tran) = @_;
    return $tran->strand == 1 ? _end_splice_site($vf, $tran) : _start_splice_site($vf, $tran); 
}

sub essential_splice_site {
    my ($vf, $tran) = @_;
    return (acceptor_splice_site($vf, $tran) or donor_splice_site($vf, $tran));
}

sub splice_region {
    my ($vf, $tran) = @_;
    
    my ($into_exon, $into_intron) = $tran->strand == 1 ? (8, 3) : (3, 8);

    for my $intron (@{ $tran->get_all_Introns }) {

        print "VF: ",$vf->start," - ",$vf->end,"\n";
        print "IN: ",$intron->start," - ",$intron->end,"\n";
        print "into_exon = $into_exon, into_intron = $into_intron\n";
        print "IN start: ",$intron->start-$into_exon," - ",$intron->start+$into_intron,"\n";
        print "IN end: ",$intron->end-$into_intron," - ",$intron->end+$into_exon,"\n";

        if ( overlap($vf->start, $vf->end, $intron->start-3, $intron->start+8) or
             overlap($vf->start, $vf->end, $intron->end-8, $intron->end+3) ) {
            return 1;
        }
    }

    return 0;
}

sub within_intron {
    my ($vf, $tran) = @_;

    for my $intron (@{ $tran->get_all_Introns }) {
        if (overlap($vf->start, $vf->end, $intron->start, $intron->end)) {
            return 1;
        }
    }

    return 0;
}

sub _start_utr {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    return overlap($vf->start, $vf->end, $tran->start, $tran->coding_region_start-1);    
}

sub _end_utr {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    return overlap($vf->start, $vf->end, $tran->coding_region_end+1, $tran->end);    
}

sub within_5_prime_utr {
    my ($vf, $tran) = @_;
    return $tran->strand == 1 ? _start_utr($vf, $tran) : _end_utr($vf, $tran); 
}

sub within_3_prime_utr {
    my ($vf, $tran) = @_;
    return $tran->strand == 1 ? _end_utr($vf, $tran) : _start_utr($vf, $tran); 
}

sub map_coords {
    my ($vf, $tran) = @_;
    my $tm = $tran->get_TranscriptMapper;
    print "vf: ".$vf->start," - ",$vf->end," ",$vf->strand,"\n";
    $tran->{_cdna_coords} ||= [ $tm->genomic2cdna($vf->start, $vf->end, $vf->strand) ];
    $tran->{_cds_coords}  ||= [ $tm->genomic2cds($vf->start, $vf->end, $vf->strand) ];
    $tran->{_pep_coords}  ||= [ $tm->genomic2pep($vf->start, $vf->end, $vf->strand) ];
}

sub complex_indel {
    my ($vf, $tran) = @_;
    
    return 0 unless $vf->var_class eq 'in-del';

    map_coords($vf, $tran);

    print "num cds coords: ".scalar(@{$tran->{_cds_coords}})."\n";
    print "num cdna coords: ".scalar(@{$tran->{_cdna_coords}})."\n";
    print "num pep coords: ".scalar(@{$tran->{_pep_coords}})."\n";

    my $tm = $tran->get_TranscriptMapper;

    my @cds_coords = $tm->genomic2cds($vf->start, $vf->end, $vf->strand);

    print "CDS:\n";
    #for my $e (@{$tran->{_cds_coords}}) {
    for my $e (@cds_coords) {
        print ref $e, "\n";
    }

    print "CDNA:\n";
    for my $e (@{$tran->{_cdna_coords}}) {
        print ref $e, "\n";
    }

    print "T cds: ".$tran->coding_region_start." - ".$tran->coding_region_end."\n";

    #return @{ $tran->{_cds_coords} } > 1;
    return (scalar(@cds_coords) > 1);
}

sub within_coding_region {
    my ($vf, $tran) = @_;
    return overlap($vf->start, $vf->end, $tran->coding_region_start, $tran->coding_region_end);
}



sub synonymous_coding {
    my ($vf, $tran) = @_;

    map_coords($vf, $tran);

    return 0 if (@{ $tran->{_pep_coords} } != 1);

#    my $mrna_seq = $tran->translateable_seq;
#    
#    my $mrna = Bio::Seq->new(
#        -seq        => $mrna_seq,
#        -moltype    => 'dna',
#        -alphabet   => 'dna'
#    );
#
#    my ($attrib) = @{ $tran->slice->get_all_Attributes('codon_table') };
#
#    my $codon_table = $attrib ? $attrib->value : 1; # default to the standard vertebrate codon table
#
#    my $peptide = $mrna->translate(undef, undef, undef, $codon_table)->seq;
    
    my $pep_coord = $tran->{_pep_coords}->[0];
    my $cds_coord = $tran->{_cds_coords}->[0];
    
    return 0 if $pep_coord->isa('Bio::EnsEMBL::Mapper::Gap');

    my $peptide = $tran->translate;

    my $pep_len = $pep_coord->end - $pep_coord->start + 1;

    my $ref_aa = substr($peptide, $pep_coord->start - 1, $pep_len);

    my @alleles = split /\//, $vf->allele_string;

    my $var_len = $vf->end - $vf->start + 1;

    for my $allele (@alleles) {
        $allele =~ s/-//; # replace insertions with an empty string

        my $allele_len = length($allele);

        if (abs($allele_len - $var_len) % 3) {
            # frameshift_coding
        }

    }


    for my $e (@{$tran->{_cds_coords}}) {
        
    }
}

# NB: most specific rules should go first

my $effect_rules = {
    'no_effect'                     => \&affects_transcript,
    '5KB_upstream_variant'          => \&upstream_5KB,
    '5KB_downstream_variant'        => \&downstream_5KB,
    '2KB_upstream_variant'          => \&upstream_2KB,
    '2KB_downstream_variant'        => \&downstream_2KB,
    'splice_donor_variant'          => \&donor_splice_site,
    'splice_acceptor_variant'       => \&acceptor_splice_site,
    'splice_region_variant'         => \&splice_region,
    'intron_variant'                => \&within_intron,
    '5_prime_UTR_variant'           => \&within_5_prime_utr,
    '3_prime_UTR_variant'           => \&within_3_prime_utr,
#    'complex_change_in_transcript'  => \&complex_indel,
};

sub transcript_effect {

    my ($vf, $tran) = @_;

    my $result = {};

    print "T: ",$tran->start, "-", $tran->end," ",$tran->strand,"\n";
    print "V: ",$vf->start, "-", $vf->end," ",$vf->strand,"\n";
    print $tran->stable_id,"\n";
    print $vf->variation_name,"\n";

    for my $effect (keys %$effect_rules) {
        my $rule = $effect_rules->{$effect};
        if ($rule->($vf, $tran)) {
            print "$effect holds\n";
        }
    }
}

sub regulatory_region_effect {
    my ($reg_region, $vf) = @_;

}

sub get_all_variation_effects {

}

1;


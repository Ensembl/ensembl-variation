package Bio::EnsEMBL::Variation::Utils::VariationEffect;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::VariationFeatureOverlap;
use Bio::EnsEMBL::Variation::TranscriptVariationNew;
use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;
use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
use Bio::EnsEMBL::Variation::OverlapConsequence;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&transcript_effect &regulatory_region_effect &get_all_variation_effects &overlap);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);

use Data::Dumper;

use Inline C => <<'END_C';

int overlap (int f1_start, int f1_end, int f2_start, int f2_end) {
    return (f1_end >= f2_start && f1_start <= f2_end);
}

END_C

#sub overlap {
#    my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
#    return ($f1_end >= $f2_start and $f1_start <= $f2_end);
#}

sub overlaps_transcript {
    my ($vf, $tran) = @_;
    
    return 0 unless $tran->isa('Bio::EnsEMBL::Transcript');
    
    return overlap($vf->start, $vf->end, $tran->start - 5000, $tran->end + 5000);
}

sub within_transcript {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    
    return overlap($vf->start, $vf->end, $tran->start, $tran->end);
}

sub _before_start {
    my ($vf, $tran, $dist) = @_;
    return ( ($vf->end >= ($tran->start - $dist)) and ($vf->end < $tran->start) );
}

sub _after_end {
    my ($vf, $tran, $dist) = @_;
    return ( ($vf->start <= ($tran->end + $dist)) and ($vf->start > $tran->end) );
}

sub _upstream {
    my ($vf, $tran, $dist) = @_;
    return $tran->strand == 1 ? _before_start($vf, $tran, $dist) : _after_end($vf, $tran, $dist);
}

sub _downstream {
    my ($vf, $tran, $allele, $dist) = @_;
    return $tran->strand == 1 ? _after_end($vf, $tran, $allele, $dist) : _before_start($vf, $tran, $allele, $dist);
}

sub upstream_5KB {
    my $tva    = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;

    return _upstream($vf, $tran, 5000);
}

sub downstream_5KB {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript; 

    return _downstream($vf, $tran, 5000);
}

sub upstream_2KB {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript; 

    return _upstream($vf, $tran, 2000);
}

sub downstream_2KB {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript; 

    return _downstream($vf, $tran, 2000);
}

sub within_nmd_transcript {
    my $tva     = shift;
    my $tran    = $tva->transcript; 

    return ( within_transcript($tva) and ($tran->biotype eq 'nonsense_mediated_decay') );
}

sub within_non_coding_gene {
    my $tva     = shift;
    my $tran    = $tva->transcript;
    
    return ( within_transcript($tva) and (not $tran->translation) );
}

sub within_miRNA {
    my $tva     = shift;
    my $tv      = $tva->transcript_variation;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    
    if ($tran->biotype eq 'miRNA') {
        my ($attribute) = @{ $tran->get_all_Attributes('miRNA') };
        
        if (defined $attribute && $attribute->value =~ /(\d+)-(\d+)/) { 
            for my $coord ($tv->mapper->cdna2genomic($1, $2, $tran->strand)) {
                if ($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                    if (overlap($vf->start, $vf->end, $coord->start, $coord->end)) {
                        return 1;
                    }
                }
            }
        }
    }
    
    return 0;
}

sub donor_splice_site {
    my $tva     = shift;
    my $tran    = $tva->transcript;
    
    return $tran->strand == 1 ? 
        $tva->transcript_variation->intron_effects->{start_splice_site} :
        $tva->transcript_variation->intron_effects->{end_splice_site};
}

sub acceptor_splice_site {
    my $tva     = shift;
    my $tran    = $tva->transcript;
    
    return $tran->strand == 1 ? 
        $tva->transcript_variation->intron_effects->{end_splice_site} :
        $tva->transcript_variation->intron_effects->{start_splice_site};
}

sub essential_splice_site {
    my $tva = shift;
    return ( acceptor_splice_site($tva) or donor_splice_site($tva) );
}

sub splice_region {
    my $tva    = shift;

    return $tva->transcript_variation->intron_effects->{splice_region};
}

sub within_intron {
    my $tva    = shift;

    return $tva->transcript_variation->intron_effects->{intronic};
}

sub within_coding_region {
    my $tva    = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    
    return 0 unless $tran->translation;
    
    return overlap($vf->start, $vf->end, $tran->coding_region_start, $tran->coding_region_end);
}

sub _before_coding {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    return overlap($vf->start, $vf->end, $tran->start, $tran->coding_region_start-1);    
}

sub _after_coding {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    return overlap($vf->start, $vf->end, $tran->coding_region_end+1, $tran->end);    
}

sub within_5_prime_utr {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    
    my $five_prime_of_coding = 
        $tran->strand == 1 ? 
            _before_coding($vf, $tran) : 
            _after_coding($vf, $tran);
    
    return ( $five_prime_of_coding and (not within_intron($tva)) );
}

sub within_3_prime_utr {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript; 
    
    my $three_prime_of_coding = 
        $tran->strand == 1 ? 
            _after_coding($vf, $tran) : 
            _before_coding($vf, $tran);
    
    return ( $three_prime_of_coding and (not within_intron($tva)) );
}

sub complex_indel {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    
    return 0 unless $vf->var_class eq 'in-del';

    return @{ $tva->transcript_variation->cds_coords } > 1;
}

sub synonymous_coding {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    return 0 unless $tva->affects_cds;
    
    my $alt_aa = $tva->amino_acid;
    my $ref_aa = $tv->reference_allele->amino_acid;
    
    my $alt_seq = $tva->seq;
    my $ref_seq = $tv->reference_allele->seq;
    
    return 0 if ($alt_seq eq '-' or $ref_seq eq '-');
    
    return ( $alt_aa eq $ref_aa );
}

sub non_synonymous_coding {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    return 0 unless $tva->affects_cds;

    my $alt_aa = $tva->amino_acid;
    my $ref_aa = $tv->reference_allele->amino_acid;
    
    my $alt_seq = $tva->seq;
    my $ref_seq = $tv->reference_allele->seq;
    
    return 0 if ($alt_seq eq '-' or $ref_seq eq '-');
    
    return ( $alt_aa ne $ref_aa );
}

sub stop_gained {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    return 0 unless $tva->affects_cds;
    
    my $alt_aa = $tva->amino_acid;
    my $ref_aa = $tv->reference_allele->amino_acid;

    return ( ($alt_aa =~ /\*/) and ($ref_aa !~ /\*/) );
}

sub stop_lost {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    return 0 unless $tva->affects_cds;
    
    my $alt_aa = $tva->amino_acid;
    my $ref_aa = $tv->reference_allele->amino_acid;

    return ( ($alt_aa !~ /\*/) and ($ref_aa =~ /\*/) );
}

sub frameshift {
    my $tva = shift;

    return 0 unless $tva->affects_cds;
    
    return 0 if partial_codon($tva);

    my $tv = $tva->transcript_variation;
    my $var_len = $tv->cds_end - $tv->cds_start + 1;

    my $seq = $tva->seq;
    
    $seq = '' if $seq eq '-';

    return abs( length($seq) - $var_len ) % 3;
}

sub partial_codon {
    my $tva = shift;
    
    return 0 unless $tva->affects_cds;

    my $tv = $tva->transcript_variation;

    my $cds_length = length $tv->translateable_seq;

    my $codon_cds_start = ($tv->pep_start * 3) - 2;

    my $last_codon_length = $cds_length - ($codon_cds_start - 1);
    
    return ( $last_codon_length < 3 and $last_codon_length > 0 );
}

sub within_coding_frameshift_intron {
    my $tva = shift;
    
    return (within_coding_region($tva) and 
        $tva->transcript_variation->intron_effects->{within_frameshift_intron});
}

my @conseq_hashes = ( 
    {
        SO_term         => '5KB_upstream_variant',
        predicate       => \&upstream_5KB,
        ensembl_term    => 'UPSTREAM',
        SO_id           => 'SO:0001635',
        is_definitive   => 1,
        rank            => 1,
    },
    {
        SO_term         => '5KB_downstream_variant',
        predicate       => \&downstream_5KB,
        ensembl_term    => 'DOWNSTREAM',
        SO_id           => 'SO:0001633',
        is_definitive   => 1,
    },
    {
        SO_term         => '2KB_upstream_variant',
        predicate       => \&upstream_2KB,
        ensembl_term    => 'UPSTREAM',
        SO_id           => 'SO:0001636',
        NCBI_term       => 'near-gene-5',
        is_definitive   => 1,
    },
    {
        SO_term         => '2KB_downstream_variant',
        predicate       => \&downstream_2KB,
        ensembl_term    => 'DOWNSTREAM',
        SO_id           => 'SO:0001634',
        NCBI_term       => 'near-gene-3',
        is_definitive   => 1,
    },
    {
        SO_term         => 'splice_donor_variant',
        predicate       => \&donor_splice_site,
        ensembl_term    => 'ESSENTIAL_SPLICE_SITE',
        SO_id           => 'SO:0001575',
        NCBI_term       => 'splice-5',
        is_definitive   => 0,
    },
    {
        SO_term         => 'splice_acceptor_variant',
        predicate       => \&acceptor_splice_site,
        ensembl_term    => 'ESSENTIAL_SPLICE_SITE',
        SO_id           => 'SO:0001574',
        NCBI_term       => 'splice-3',
        is_definitive   => 0,
    },
    {
        SO_term         => 'splice_region_variant',
        predicate       => \&splice_region,
        ensembl_term    => 'SPLICE_SITE',
        SO_id           => 'SO:0001630',
        is_definitive   => 0,
    },
    {
        SO_term         => 'intron_variant',
        predicate       => \&within_intron,
        ensembl_term    => 'INTRONIC',
        SO_id           => 'SO:0001627',
        NCBI_term       => 'intron',
        is_definitive   => 0,
    },
    {
        SO_term         => '5_prime_UTR_variant',
        predicate       => \&within_5_prime_utr,
        ensembl_term    => '5PRIME_UTR',
        SO_id           => 'SO:0001623',
        NCBI_term       => 'untranslated-5',
        is_definitive   => 0,
    }, 
    {
        SO_term         => '3_prime_UTR_variant',
        predicate       => \&within_3_prime_utr,
        ensembl_term    => '3PRIME_UTR',
        SO_id           => 'SO:0001624',
        NCBI_term       => 'untranslated-3',
        is_definitive   => 0,
    },
    {
        SO_term         => 'complex_change_in_transcript',
        predicate       => \&complex_indel,
        ensembl_term    => 'COMPLEX_INDEL',
        SO_id           => 'SO:0001577',
        is_definitive   => 0,
    },
    {
        SO_term         => 'synonymous_codon',
        predicate       => \&synonymous_coding,
        ensembl_term    => 'SYNONYMOUS_CODING',
        SO_id           => 'SO:0001588',
        NCBI_term       => 'cds-synon',
        is_definitive   => 0,
    },
    {
        SO_term         => 'non_synonymous_codon',
        predicate       => \&non_synonymous_coding,
        ensembl_term    => 'NON_SYNONYMOUS_CODING',
        SO_id           => 'SO:0001583',
        NCBI_term       => 'missense',
        is_definitive   => 0,
    },
    {
        SO_term         => 'stop_gained',
        predicate       => \&stop_gained,
        ensembl_term    => 'STOP_GAINED',
        SO_id           => 'SO:0001587',
        NCBI_term       => 'nonsense',
        is_definitive   => 0,
    },
    {
        SO_term         => 'stop_lost',
        predicate       => \&stop_lost,
        ensembl_term    => 'STOP_LOST',
        SO_id           => 'SO:0001578',
        is_definitive   => 0,
    },
    {
        SO_term         => 'frameshift_variant',
        predicate       => \&frameshift,
        ensembl_term    => 'FRAMESHIFT_CODING',
        SO_id           => 'SO:0001589',
        NCBI_term       => 'frameshift',
        is_definitive   => 0,
    },
    {
        SO_term         => 'incomplete_terminal_codon_variant',
        predicate       => \&partial_codon,
        ensembl_term    => 'PARTIAL_CODON',
        SO_id           => 'SO:0001626',
        is_definitive   => 1,
    },
    {
        SO_term         => 'NMD_transcript_variant',
        predicate       => \&within_nmd_transcript,
        ensembl_term    => 'NMD_TRANSCRIPT',
        SO_id           => 'SO:0001621',
        is_definitive   => 1,
    },
    {
        SO_term         => 'nc_transcript_variant',
        predicate       => \&within_non_coding_gene,
        ensembl_term    => 'WITHIN_NON_CODING_GENE',
        SO_id           => 'SO:0001619',
        is_definitive   => 0,
    },
    {
        SO_term         => 'mature_miRNA_variant',
        predicate       => \&within_miRNA,
        ensembl_term    => 'WITHIN_MATURE_miRNA',
        SO_id           => 'SO:0001620',
        is_definitive   => 0,
    },
    {
        SO_term         => 'coding_sequence_variant',
        predicate       => \&within_coding_frameshift_intron,
        ensembl_term    => 'CODING_UNKNOWN',
        SO_id           => 'SO:0001580',
        is_definitive   => 0,
    },
);

# create OverlapConsequence objects which store details about each consequence of
# interest along with a predicate that can be used to test if the consequence holds
# for a particular VariationFeatureOverlapAllele

# sort the list so that definitive consequences (i.e. those that rule out all other
# consequences) are tested first, to minimise redundant tests

my @consequences = sort { $b->is_definitive <=> $a->is_definitive } map {
    Bio::EnsEMBL::Variation::OverlapConsequence->new_fast($_)
} @conseq_hashes;

sub transcript_effect {

    my ($vf, $tran) = @_;

    # unless this vf overlaps the transcript, it has no effect

    return undef unless overlaps_transcript($vf, $tran);
 
    # create a TranscriptVariation object representing this overlap
 
    my $tv = Bio::EnsEMBL::Variation::TranscriptVariationNew->new_fast({            
        variation_feature  => $vf,
        feature            => $tran,
    });
    
    # we now look at each allele of the VariationFeature in turn and calculate
    # its effect on the transcript
    
    # get the allele string, expand it, and split it into separate alleles
    
    my $allele_string = $vf->allele_string;
    
    expand(\$allele_string);
    
    my @alleles = split /\//, $allele_string;

    # if the strands differ we need to reverse complement the allele sequences
    
    unless ($tran->strand == $vf->strand) {
        map { reverse_comp(\$_) } @alleles;
    }
    
    # create an object representing the reference allele
    
    my $ref_allele = shift @alleles;

    my $ref_tva = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
        variation_feature_overlap   => $tv,
        seq                         => $ref_allele,
        is_reference                => 1,
    });
    
    $tv->reference_allele($ref_tva);

    # create objects representing the alternate alleles, and calculate their 
    # effect on the transcript using the rules above

    for my $allele (@alleles) {
        
        my $tva = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
            variation_feature_overlap   => $tv,
            seq                         => $allele,
        });
        
        $tva->calc_consequences(\@consequences);
       
        $tv->alt_alleles($tva);
    }

    # finally, return the TranscriptVariation object that now stores all of this information

    return $tv;
}

sub regulatory_region_effect {
    my ($reg_region, $vf) = @_;

}

sub get_all_variation_effects {

    my ($vf) = @_;

    # get all features that overlap this variation
    #
    #

}

1;


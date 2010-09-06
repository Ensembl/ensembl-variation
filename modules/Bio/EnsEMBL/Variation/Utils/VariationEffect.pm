package Bio::EnsEMBL::Variation::Utils::VariationEffect;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::VariationFeatureOverlap;
use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;
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
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature;
    
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
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    return _upstream($vf, $tran, 5000);
}

sub downstream_5KB {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    return _downstream($vf, $tran, 5000);
}

sub upstream_2KB {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    return _upstream($vf, $tran, 2000);
}

sub downstream_2KB {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    return _downstream($vf, $tran, 2000);
}

sub within_nmd_transcript {
    my $vfoa    = shift;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    return ( within_transcript($vfoa) and ($tran->biotype eq 'nonsense_mediated_decay') );
}

sub within_non_coding_gene {
    my $vfoa    = shift;
    my $tran    = $vfoa->variation_feature_overlap->feature;
    
    return ( within_transcript($vfoa) and (not $tran->translation) );
}

sub within_miRNA {
    my $vfoa    = shift;
    my $vfo     = $vfoa->variation_feature_overlap;
    my $vf      = $vfo->variation_feature;
    my $tran    = $vfo->feature;
    
    if ($tran->biotype eq 'miRNA') {
        my ($attribute) = @{ $tran->get_all_Attributes('miRNA') };
        
        if (defined $attribute && $attribute->value =~ /(\d+)-(\d+)/) { 
            for my $coord ($vfo->mapper->cdna2genomic($1, $2, $tran->strand)) {
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
    my $vfoa    = shift;
    my $tran    = $vfoa->variation_feature_overlap->feature;
    
    return $tran->strand == 1 ? 
        $vfoa->variation_feature_overlap->intron_effects->{start_splice_site} :
        $vfoa->variation_feature_overlap->intron_effects->{end_splice_site};
}

sub acceptor_splice_site {
    my $vfoa    = shift;
    my $tran    = $vfoa->variation_feature_overlap->feature;
    
    return $tran->strand == 1 ? 
        $vfoa->variation_feature_overlap->intron_effects->{end_splice_site} :
        $vfoa->variation_feature_overlap->intron_effects->{start_splice_site};
}

sub essential_splice_site {
    my $vfoa = shift;
    return ( acceptor_splice_site($vfoa) or donor_splice_site($vfoa) );
}

sub splice_region {
    my $vfoa    = shift;

    return $vfoa->variation_feature_overlap->intron_effects->{splice_region};
}

sub within_intron {
    my $vfoa    = shift;

    return $vfoa->variation_feature_overlap->intron_effects->{intronic};
}

sub within_coding_region {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature;
    
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
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
    
    my $five_prime_of_coding = 
        $tran->strand == 1 ? 
            _before_coding($vf, $tran) : 
            _after_coding($vf, $tran);
    
    return ( $five_prime_of_coding and (not within_intron($vfoa)) );
}

sub within_3_prime_utr {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
    
    my $three_prime_of_coding = 
        $tran->strand == 1 ? 
            _after_coding($vf, $tran) : 
            _before_coding($vf, $tran);
    
    return ( $three_prime_of_coding and (not within_intron($vfoa)) );
}

sub complex_indel {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    
    return 0 unless $vf->var_class eq 'in-del';

    return @{ $vfoa->variation_feature_overlap->cds_coords } > 1;
}

sub synonymous_coding {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;
    
    my $alt_aa = $vfoa->aa;
    my $ref_aa = $vfoa->variation_feature_overlap->reference_allele->aa;
    
    my $alt_seq = $vfoa->seq;
    my $ref_seq = $vfoa->variation_feature_overlap->reference_allele->seq;
    
    return 0 if ($alt_seq eq '-' or $ref_seq eq '-');
    
    return ( $alt_aa eq $ref_aa );
}

sub non_synonymous_coding {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;

    my $alt_aa = $vfoa->aa;
    my $ref_aa = $vfoa->variation_feature_overlap->reference_allele->aa;
    
    my $alt_seq = $vfoa->seq;
    my $ref_seq = $vfoa->variation_feature_overlap->reference_allele->seq;
    
    return 0 if ($alt_seq eq '-' or $ref_seq eq '-');
    
    return ( $alt_aa ne $ref_aa );
}

sub stop_gained {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;
    
    my $alt_aa = $vfoa->aa;
    my $ref_aa = $vfoa->variation_feature_overlap->reference_allele->aa;

    return ( ($alt_aa =~ /\*/) and ($ref_aa !~ /\*/) );
}

sub stop_lost {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;
    
    my $alt_aa = $vfoa->aa;
    my $ref_aa = $vfoa->variation_feature_overlap->reference_allele->aa;

    return ( ($alt_aa !~ /\*/) and ($ref_aa =~ /\*/) );
}

sub frameshift {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;
    
    return 0 if partial_codon($vfoa);

    my $vfo = $vfoa->variation_feature_overlap;
    my $var_len = $vfo->cds_end - $vfo->cds_start + 1;

    my $seq = $vfoa->seq;
    
    $seq = '' if $seq eq '-';

    return abs( length($seq) - $var_len ) % 3;
}

sub partial_codon {
    my $vfoa = shift;
    
    return 0 unless $vfoa->affects_cds;

    my $vfo = $vfoa->variation_feature_overlap;

    my $cds_length = length $vfo->translateable_seq;

    my $codon_cds_start = ($vfo->pep_start * 3) - 2;

    my $last_codon_length = $cds_length - ($codon_cds_start - 1);
    
    return ( $last_codon_length < 3 and $last_codon_length > 0 );
}

sub within_coding_frameshift_intron {
    my $vfoa = shift;
    
    return (within_coding_region($vfoa) and 
        $vfoa->variation_feature_overlap->intron_effects->{within_frameshift_intron});
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
 
    # we cache some slow-to-compute features in the transcript object itself
    # (this means this data will get garbage collected along with the transcript
    # and we don't have to worry about clearing up our own cache) but if
    # this is the first time we have seen this transcript we will need to
    # build the cache first
 
    my $tran_features = $tran->{_variation_effect_feature_cache};
    
    unless ($tran_features) {
        
        $tran_features = {
            introns => $tran->get_all_Introns,
            mapper  => $tran->get_TranscriptMapper,
        };
        
        if (my $translation = $tran->translate) {
             $tran_features->{translateable_seq} = $tran->translateable_seq;
             $tran_features->{peptide} = $translation->seq;
        }

        $tran->{_variation_effect_feature_cache} = $tran_features;
    }

    # a VariationFeatureOverlap represents the combination of a VariationFeature
    # and (in this case) a Transcript, and stores various cached features which 
    # are used by the consequence predicates above

    my $vfo = Bio::EnsEMBL::Variation::VariationFeatureOverlap->new_fast({            
        variation_feature  => $vf,
        feature            => $tran,
        feature_type       => 'transcript',
        translateable_seq  => $tran_features->{translateable_seq},
        peptide            => $tran_features->{peptide},
        introns            => $tran_features->{introns},
        mapper             => $tran_features->{mapper},
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

    my $ref_vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
        variation_feature_overlap   => $vfo,
        seq                         => $ref_allele,
        is_reference                => 1,
    });
    
    $vfo->reference_allele($ref_vfoa);

    # create objects representing the alternate alleles, and calculate their 
    # effect on the transcript using the rules above

    for my $allele (@alleles) {
        
        my $vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
            variation_feature_overlap   => $vfo,
            seq                         => $allele,
        });
        
        $vfoa->calc_consequences(\@consequences);
       
        $vfo->alt_alleles($vfoa);
    }

    # finally, return the VariationFeatureOverlap object that now stores all of this information

    return $vfo;
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


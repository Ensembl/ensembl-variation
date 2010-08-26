package Bio::EnsEMBL::Variation::Utils::VariationEffect;

use strict;
use warnings;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&transcript_effect &regulatory_region_effect &get_all_variation_effects);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Data::Dumper;

sub overlap {
    my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
    return ($f1_end >= $f2_start and $f1_start <= $f2_end);
}

sub overlaps_transcript {
    my ($vf, $tran) = @_;
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
    return ($vf->end >= ($tran->start - $dist) and $vf->end < $tran->start);
}

sub _after_end {
    my ($vf, $tran, $dist) = @_;
    return ($vf->start <= ($tran->end + $dist) and $vf->start > $tran->end);
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
   
    return 0 unless $tran->isa('Bio::EnsEMBL::Transcript');

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
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

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
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
    
    return $tran->strand == 1 ? _start_splice_site($vf, $tran) : _end_splice_site($vf, $tran);
}

sub acceptor_splice_site {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
    
    return $tran->strand == 1 ? _end_splice_site($vf, $tran) : _start_splice_site($vf, $tran); 
}

sub essential_splice_site {
    my $vfoa    = shift;
    return (acceptor_splice_site($vfoa) or donor_splice_site($vfoa));
}

sub splice_region {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    for my $intron (@{ $tran->get_all_Introns }) {

        if ( overlap($vf->start, $vf->end, $intron->start-3, $intron->start+8) or
             overlap($vf->start, $vf->end, $intron->end-8, $intron->end+3) ) {
            return 1;
        }
    }

    return 0;
}

sub within_intron {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

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
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
    
    return $tran->strand == 1 ? _start_utr($vf, $tran) : _end_utr($vf, $tran); 
}

sub within_3_prime_utr {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
    
    return $tran->strand == 1 ? _end_utr($vf, $tran) : _start_utr($vf, $tran); 
}

sub complex_indel {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    
    return 0 unless $vf->var_class eq 'in-del';

    return @{ $vfoa->variation_feature_overlap->cds_coords } > 1;
}

my $_coding_effect_cache = {};

sub _calc_coding_effect {
    
    my $vfoa = shift;

    my $allele_effects = $_coding_effect_cache->{$tran}->{$vf} ||= {};

    unless (keys %$allele_effects) {

        my @pep_coords = _map_coords($vf, $tran, 'pep');
        my @cds_coords = _map_coords($vf, $tran, 'cds');
    
        return unless @pep_coords == 1;

        my $pep_coord = shift @pep_coords;
        my $cds_coord = shift @cds_coords;

        return if $pep_coord->isa('Bio::EnsEMBL::Mapper::Gap');

        my $peptide = $tran->translate->seq;
        my $cds     = $tran->translateable_seq;

        my $var_pep_len = $pep_coord->end - $pep_coord->start + 1;

        my $ref_aa = substr($peptide, $pep_coord->start - 1, $var_pep_len);

        print "peptide:\n$peptide\nvar_pep_len: $var_pep_len\nref aa: $ref_aa \n";
        print "pep_coords: ",$pep_coord->start," - ", $pep_coord->end, "\n";

        my @alleles = split /\//, $vf->allele_string;

        unless ($vf->strand == $tran->strand) {
            map { reverse_comp(\$_) } @alleles;
        }

        my $var_nt_len = $cds_coord->end - $cds_coord->start + 1;

        my $codon_cds_start = $pep_coord->start * 3 - 2;
        my $codon_cds_end   = $pep_coord->end * 3;
        my $codon_len       = $codon_cds_ends - $codon_cds_start + 1;

        my $ref_codon = substr($cds, $codon_cds_start-1, $codon_len);

        for my $allele (@alleles) {
            $allele =~ s/-//; # replace insertions with an empty string

            my $allele_len = length($allele);

            if (abs($allele_len - $var_nt_len) % 3) {
                # frameshift_coding
            }

            my $new_cds = $cds;

            substr($new_cds, $cds_coord->start-1, $allele_len) = $allele;

            my $new_codon = substr($new_cds, $codon_cds_start-1, $codon_len + ($allele_len - $var_nt_len));

            my $new_codon_seq = Bio::Seq->new(
                -seq        => $new_codon,
                -moltype    => 'dna',
                -alphabet   => 'dna'
            );

            my $new_aa = $new_codon_seq->translate;

            $new_aa = '-' if length($new_aa) < 1;

            if (lc($new_aa) ne lc($ref_aa)) {
                if ($new_aa =~ /\*/ and $ref_aa !~ /\*/) {
                    # stop gained
                }
                elsif ($ref_aa =~ /\*/ and $ref)
            }
            else {
                # synonymous coding
                return 1;
            }
        }
    }

    return $allele ? $allele_effect->{$allele} : $allele_effect;
}

sub synonymous_coding {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;
    
    return $vfoa->aa eq $vfoa->variation_feature_overlap->reference_allele->aa;
}

sub non_synonynmous_coding {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;

    return $vfoa->aa ne $vfoa->variation_feature_overlap->reference_allele->aa;
}

sub stop_gained {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;

    return $vfoa->aa =~ /\*/ and $vfoa->variation_feature_overlap->reference_allele->aa !~ /\*/;
}

sub stop_lost {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;

    return $vfoa->aa !~ /\*/ and $vfoa->variation_feature_overlap->reference_allele->aa =~ /\*/;
}

sub frameshift {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;

    my $vfo = $vfoa->variation_feature_overlap;
    my $var_len = $vfo->cds_end - $vfo->cds_start + 1;

    return abs( length($vfoa->seq) - $var_len ) % 3;
}

# NB: most specific rules should go first

my $effect_rules = {
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
    'complex_change_in_transcript'  => \&complex_indel,
    'synonymous_codon'              => \&synonymous_coding,
    'non_synonymous_codon'          => \&non_synonymous_coding,
    'frameshift_variant'            => \&frameshift,
};

my @consequences;

for my $effect (keys %$effect_rules) {
    my $predicate = $effect_rules->{$effect};
    push @consequences, Bio::EnsEMBL::Variation::Consequence->new_fast({
            SO_term     => $effect,
            predicate   => $predicate,
        });
}

sub transcript_effect {

    my ($vf, $tran) = @_;

    my $result = {};

    print "T: ",$tran->start, "-", $tran->end," ",$tran->strand,"\n";
    print "V: ",$vf->start, "-", $vf->end," ",$vf->strand,"\n";
    print $tran->stable_id,"\n";
    print $vf->variation_name,"\n";

    if (overlaps_transcript($vf, $tran)) {

        my $vfo = Bio::EnsEMBL::Variation::VariationFeatureOverlap->new_fast({            
            variation_feature  => $vf,
            feature            => $transcript,
            feature_type       => 'transcript',
        });
        
        my $tm = $tran->get_TranscriptMapper;

        $vfo->cdna_coords([$tm->genomic2cdna($vf->start, $vf->end, $tran->strand)]);
        $vfo->cds_coords([$tm->genomic2cds($vf->start, $vf->end, $tran->strand)]);
        $vfo->pep_coords([$tm->genomic2pep($vf->start, $vf->end, $tran->strand)]);

        my @alleles = split /\//, $vf->allele_string;

        my $ref_allele = shift @alleles;

        my $ref_vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
            variation_feature_overlap   => $vfo
            seq                         => $ref_allele,
            is_reference                => 1,
        });

        $vfo->reference_allele($ref_vfoa);

        for my $allele (@alleles) {
            my $vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
                variation_feature_overlap   => $vfo
                seq                         => $allele,
            });
            $vfoa->calc_consequences(\@consequences);
            $vfo->alleles($vfoa);
        }

        return $vfo;
    }
    else {
        return undef;
    }
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


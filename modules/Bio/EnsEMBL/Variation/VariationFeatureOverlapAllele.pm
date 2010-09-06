package Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

use strict;
use warnings;

use Scalar::Util qw(weaken);

sub new {
    my ($class) = @_;
    return bless {}, $class;
}

sub new_fast {
    my ($class, $hashref) = @_;
    my $self = bless $hashref, $class;
    # avoid a memory leak, because the vfo also has a reference to us
    weaken $self->{variation_feature_overlap} if $self->{variation_feature_overlap};
    return $self;
}

sub variation_feature_overlap {
    my ($self, $variation_feature_overlap) = @_;
    if ($variation_feature_overlap) {
        $self->{variation_feature_overlap} = $variation_feature_overlap;
        # avoid a memory leak, because the vfo also has a reference to us
        weaken $self->{variation_feature_overlap}
    }
    return $self->{variation_feature_overlap};
}

sub seq {
    my ($self, $seq) = @_;
    $self->{seq} = $seq if $seq;
    return $self->{seq};
}

sub aa {
    my ($self, $aa) = @_;
    $self->{aa} = $aa if $aa;
    return $self->{aa};
}

sub codon {
    my ($self, $codon) = @_;
    $self->{codon} = $codon if $codon;
    return $self->{codon};
}

sub affects_cds {
    my ($self) = @_;

    unless (defined $self->{affects_cds}) {
        
        my $vfo = $self->variation_feature_overlap;
       
        my @pep_coords = @{ $vfo->pep_coords };
       
        if (not $vfo->feature->isa('Bio::EnsEMBL::Transcript')) {
            $self->{affects_cds} = 0;
        }
        elsif (@pep_coords != 1) {
            $self->{affects_cds} = 0;
        }
        elsif ($pep_coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
            $self->{affects_cds} = 0;
        }
        else {

            # by now we have established that the vf falls entirely in the CDS, 
            # so store the various start and end coords we need

            my $tran = $vfo->feature;
            
            my ($cds, $codon_cds_start, $codon_cds_end, $codon_len);
            
            # there's some stuff we only need to calculate once per vf, rather
            # than once per allele
            
            if (defined $vfo->pep_start) {
                
                # but we do need these in any case
                
                $codon_cds_start = $vfo->pep_start * 3 - 2;
                $codon_cds_end   = $vfo->pep_end * 3;
                $codon_len       = $codon_cds_end - $codon_cds_start + 1;
                $cds             = $vfo->translateable_seq;
            }
            else {
                
                # this is the stuff we only need calculate once per vf
                
                my $exon_phase = $tran->start_Exon->phase;
    
                my $pep_coord = $vfo->pep_coords->[0];
                $vfo->pep_start($pep_coord->start);
                $vfo->pep_end($pep_coord->end);
    
                my $cds_coord = $vfo->cds_coords->[0];
                $vfo->cds_start($cds_coord->start + ($exon_phase > 0 ? $exon_phase : 0));
                $vfo->cds_end($cds_coord->end + ($exon_phase > 0 ? $exon_phase : 0));
    
                my $cdna_coord = $vfo->cdna_coords->[0];
                $vfo->cdna_start($cdna_coord->start);
                $vfo->cdna_end($cdna_coord->end);
                
                # find the reference amino acid and codon
    
                my $peptide = $vfo->peptide;
                $cds        = $vfo->translateable_seq;
 
                my $var_pep_len = $vfo->pep_end - $vfo->pep_start + 1;
    
                my $ref_aa = substr($peptide, $vfo->pep_start - 1, $var_pep_len);
            
                $ref_aa = '-' unless $ref_aa;
    
                #print "PEPTIDE:\n$peptide\nstart: ".$vfo->pep_start." end: ".$vfo->pep_end." len: ".$var_pep_len."\nREFERENCE AA: $ref_aa\n";
    
                $vfo->reference_allele->aa($ref_aa);
    
                $codon_cds_start = $vfo->pep_start * 3 - 2;
                $codon_cds_end   = $vfo->pep_end * 3;
                $codon_len       = $codon_cds_end - $codon_cds_start + 1;
    
                my $ref_codon = substr($cds, $codon_cds_start-1, $codon_len);
                
                $ref_codon = '-' unless $ref_codon;
                
                $vfo->reference_allele->codon($ref_codon);
            }
            
            # everything else does need to be calculated for each allele
            
            my $vf_nt_len = $vfo->cds_end - $vfo->cds_start + 1;
            
            # calculate the effect this allele has on the reference

            my $seq = $self->seq;

            $seq =~ s/-//;
            
            if ($seq) {
                my $new_cds = $cds;
                
                my $allele_len = length($seq);
    
                substr($new_cds, $vfo->cds_start-1, $allele_len) = $seq;
                
                my $new_codon = substr($new_cds, $codon_cds_start-1, $codon_len + ($allele_len - $vf_nt_len));
    
                $self->codon($new_codon);
    
                my $new_codon_seq = Bio::Seq->new(
                    -seq        => $new_codon,
                    -moltype    => 'dna',
                    -alphabet   => 'dna'
                );
    
                my $new_aa = $new_codon_seq->translate->seq;
    
                $new_aa = '-' unless $new_aa;
                
                $self->aa($new_aa);
            }
            else {
                $self->codon('-');
                $self->aa('-');
            }
          
            $self->{affects_cds} = 1;
        }
    }

    return $self->{affects_cds};
}

sub consequences {
    my ($self, @new_consequences) = @_;
    
    $self->{consequences} ||= [];

    if (@new_consequences) {
        my $consequences = $self->{consequences};
        push @$consequences, @new_consequences; 
    }
    
    return $self->{consequences};
}

sub is_reference {
    my ($self, $is_reference) = @_;
    $self->{is_reference} = $is_reference if defined $is_reference;
    return $self->{is_reference};
}

sub calc_consequences {
    my ($self, $consequences) = @_;

    for my $consequence (@$consequences) {
        if ($consequence->predicate->($self)) {
            $self->consequences($consequence);
            last if $consequence->is_definitive;
        }
    }
}

1;

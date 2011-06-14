=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;
    
    my $vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new(
        -variation_feature_overlap  => $vfo,
        -variation_feature_seq      => 'A',
        -is_reference               => 0,
    );

    print "sequence with respect to the feature: ", $vfoa->feature_seq, "\n";
    print "sequence with respect to the variation feature: ", $vfoa->variation_feature_seq, "\n";
    print "consequence SO terms: ", (join ",", map { $_->SO_term } @{ $vfoa->get_all_OverlapConsequences }), "\n";

=head1 DESCRIPTION

A VariationFeatureOverlapAllele object represents a single allele of a 
VariationFeatureOverlap. It is the super-class of various feature-specific allele
classes such as TranscriptVariationAllele and RegulatoryFeatureVariationAllele and 
contains methods not specific to any particular feature type. Ordinarily you will 
not create these objects yourself, but instead you would create e.g. a 
TranscriptVariation object which will then create VariationFeatureOverlapAlleles 
based on the allele string of the associated VariationFeature. 

=cut

package Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);
use Scalar::Util qw(weaken);

our $UNAMBIGUOUS_NUCLEOTIDES = qr/^[ACGT-]+$/i;

our $ALL_NUCLEOTIDES = qr/^[ACGTUMRWSYKVHDBXN-]+$/i;

our $SPECIFIED_LENGTH = qr /(\d+) BP (INSERTION|DELETION)/i;

=head2 new

  Arg [-VARIATION_FEATURE_OVERLAP] : 
    The Bio::EnsEMBL::VariationFeatureOverlap with which this allele is 
    associated

  Arg [-VARIATION_FEATURE_SEQ] :
    The allele sequence with respect to the associated VariationFeature

  Arg [-IS_REFERENCE] :
    A flag indicating if this allele is the reference allele or not

  Example : 
    my $vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new(
        -variation_feature_ovelap   => $vfo,
        -variation_feature_seq      => 'A',
        -is_reference               => 0
    );

  Description: Constructs a new VariationFeatureOverlapAllele instance given a 
               VariationFeatureOverlap and the sequence of the allele
  Returntype : A new Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele instance 
  Exceptions : throws unless both VARIATION_FEATURE_OVERLAP and VARIATION_FEATURE_SEQ 
               are supplied
  Status     : At Risk

=cut 

sub new {
    my $class = shift;

    my (
        $variation_feature_overlap, 
        $variation_feature_seq, 
        $is_reference, 
    ) = rearrange([qw(
            VARIATION_FEATURE_OVERLAP 
            VARIATION_FEATURE_SEQ 
            IS_REFERENCE
        )], @_);

    assert_ref($variation_feature_overlap, 'Bio::EnsEMBL::Variation::VariationFeatureOverlap');

    throw("Allele sequence required (variation "+$variation_feature_overlap->variation_feature->variation_name+")") 
        unless $variation_feature_seq;

    my $self = bless {
        variation_feature_overlap   => $variation_feature_overlap,
        variation_feature_seq       => $variation_feature_seq,
        is_reference                => $is_reference,
    }, $class;
    
    # avoid a memory leak, because the vfo also has a reference to us
    weaken $self->{variation_feature_overlap};

    return $self;
}

sub new_fast {
    my ($class, $hashref) = @_;
    my $self = bless $hashref, $class;
    # avoid a memory leak, because the vfo also has a reference to us
    weaken $self->{variation_feature_overlap} if $self->{variation_feature_overlap};
    return $self;
}

sub dbID {
    my ($self, $dbID) = @_;
    $self->{dbID} = $dbID if defined $dbID;
    return $self->{dbID};
}

=head2 variation_feature_overlap

  Description: Get/set the associated VariationFeatureOverlap
  Returntype : Bio::EnsEMBL::Variation::VariationFeatureOverlap
  Exceptions : throws if the argument is the wrong type
  Status     : At Risk

=cut

sub variation_feature_overlap {
    my ($self, $variation_feature_overlap) = @_;
    
    if ($variation_feature_overlap) {
        assert_ref($variation_feature_overlap, 'Bio::EnsEMBL::Variation::VariationFeatureOverlap');
        $self->{variation_feature_overlap} = $variation_feature_overlap;
        # avoid a memory leak, because the vfo also has a reference to us
        weaken $self->{variation_feature_overlap};
    }
    
    return $self->{variation_feature_overlap};
}

=head2 variation_feature

  Description: Get the associated VariationFeature
  Returntype : Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : none
  Status     : At Risk

=cut

sub variation_feature {
    my $self = shift;
    return $self->variation_feature_overlap->variation_feature;
}

=head2 feature

  Description: Get the associated Feature
  Returntype : Bio::EnsEMBL::Feature
  Exceptions : none
  Status     : At Risk

=cut

sub feature {
    my $self = shift;
    return $self->variation_feature_overlap->feature;
}

=head2 feature_seq

  Description: Get the sequence of this allele relative to the associated Feature.
               This will be the same as the variation_feature_seq when the associated
               VariationFeature is on the same strand as the Feature, or the reverse
               complement when the strands differ.
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub feature_seq {
    my $self = shift;
    
    unless ($self->{feature_seq}) {
        
        # check if we need to reverse complement the variation_feature_seq
        
        if (($self->variation_feature->strand != $self->feature->strand) && $self->seq_is_dna) {
            my $vf_seq = $self->variation_feature_seq;
            reverse_comp(\$vf_seq);
            $self->{feature_seq} = $vf_seq;
        }
        else {
            $self->{feature_seq} = $self->{variation_feature_seq};
        }
    }
    
    return $self->{feature_seq};
}

=head2 variation_feature_seq

  Args [1]   : The allele sequence relative to the VariationFeature
  Description: Get/set the sequence of this allele relative to the associated VariationFeature.
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub variation_feature_seq {
    # the sequence of this allele relative to the variation feature
    my ($self, $variation_feature_seq) = @_;
    $self->{variation_feature_seq} = $variation_feature_seq if $variation_feature_seq;
    return $self->{variation_feature_seq};
}

=head2 seq_is_unambiguous_dna

  Description: identify if the sequence of this allele is unambiguous DNA 
               i.e. if we can meaningfully translate it
  Returntype : bool
  Exceptions : none
  Status     : At Risk

=cut

sub seq_is_unambiguous_dna {
    my $self = shift;

    unless (defined $self->{seq_is_unambiguous_dna}) {
        $self->{seq_is_unambiguous_dna} = 
            $self->{variation_feature_seq} =~ /$UNAMBIGUOUS_NUCLEOTIDES/ ? 1 : 0;
    }
    
    return $self->{seq_is_unambiguous_dna};
}

=head2 seq_is_dna

  Description: identify if the sequence of this allele is DNA including ambiguity 
               codes, use seq_is_unambiguous_dna to check for alleles that do not
               include ambiguity codes
  Returntype : bool
  Exceptions : none
  Status     : At Risk

=cut

sub seq_is_dna {
    my $self = shift;

    unless (defined $self->{seq_is_dna}) {
        $self->{seq_is_dna} = 
            $self->{variation_feature_seq} =~ /$ALL_NUCLEOTIDES/ ? 1 : 0;
    }
    
    return $self->{seq_is_dna};
}

=head2 seq_length

  Description: return the length of this allele sequence, this is better than
               just using length($vfoa->feature_seq) because we check if the
               sequence is valid DNA, and also look for allele strings like 
               "(3 BP INSERTION)" to determine the length
  Returntype : int or undef if we cannot determine the length
  Exceptions : none
  Status     : At Risk

=cut

sub seq_length {
    my $self = shift;

    my $seq = $self->variation_feature_seq;

    if ($self->seq_is_dna) {
        if ($seq eq '-') {
            return 0;
        }
        else {
            return length($seq);
        }
    }
    elsif ($seq =~ /$SPECIFIED_LENGTH/) {
        return $1;
    }
    
    return undef;
}

=head2 is_reference

  Args [1]   : A boolean value 
  Description: Get/set a flag indicating if this allele is the reference allele
  Returntype : bool
  Exceptions : none
  Status     : At Risk

=cut

sub is_reference {
    my ($self, $is_reference) = @_;
    $self->{is_reference} = $is_reference if defined $is_reference;
    return $self->{is_reference};
}

=head2 allele_string

  Description: Return a '/' delimited string of the reference allele variation_feature_seq 
               and the variation_feature_seq of this allele
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub allele_string {
    my $self = shift;
    
    my $ref = $self->variation_feature_overlap->get_reference_VariationFeatureOverlapAllele->variation_feature_seq;
    
    # for the HGMDs and CNV probes where the alleles are artificially set to be
    # the same, just return the reference sequence
    
    if ($ref eq $self->variation_feature_seq) {
        return $ref;
    }
    else {
        return $ref.'/'.$self->variation_feature_seq;
    }
}

=head2 get_all_OverlapConsequences

  Description: Get a list of all the OverlapConsequences of this allele, calculating them 
               on the fly if necessary
  Returntype : listref of Bio::EnsEMBL::Variation::OverlapConsequence objects
  Exceptions : none
  Status     : At Risk

=cut

sub get_all_OverlapConsequences {
    my $self = shift;
    
    unless ($self->{overlap_consequences}) {

        # calculate consequences on the fly
        
        my $cons = [];
        
        for my $oc (values %OVERLAP_CONSEQUENCES) {
            if ($oc->feature_class eq ref $self->feature) {
                if ($oc->predicate->($self)) {
                    push @$cons, $oc;
                }
            }
        }

        $self->{overlap_consequences} = $cons;
    }
    
    return $self->{overlap_consequences};
}

=head2 add_OverlapConsequence

  Arg [1]    : Bio::EnsEMBL::Variation::OverlapConsequence instance
  Description: Add an OverlapConsequence to this allele's list 
  Returntype : none
  Exceptions : throws if the argument is the wrong type
  Status     : At Risk

=cut

sub add_OverlapConsequence {
    my ($self, $oc) = @_;
    assert_ref($oc, 'Bio::EnsEMBL::Variation::OverlapConsequence');
    $self->{overlap_consequences} ||= [];
    push @{ $self->{overlap_consequences} }, $oc;
}

sub SO_isa {
    my ($self, $query) = @_;
    
    if (my $adap = $self->variation_feature_overlap->{adaptor}) {
        if (my $ota = $adap->db->dnadb->get_OntologyTermAdaptor) {
            my $term = $ota->fetch_by_accession();
            my @parents = $ota->fetch_by_child_term($term);
        }
    }
    
    for my $cons (@{ $self->get_all_OverlapConsequences }) {
        if ($cons->SO_term eq $query) {
            return 1;
        }
    } 
}


sub _convert_to_sara {
    my $self = shift;
    
    my $oc = Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
        'label'        => 'SARA',
        'description'  => 'Same as reference allele',
        'rank'         => '99',
        'display_term' => 'SARA',
        'SO_term'      => 'SARA',
    });
    
    $self->add_OverlapConsequence($oc);
    
    return $self;
}

1;

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

=head1 NAME

Bio::EnsEMBL::Variation::OverlapConsequence

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::OverlapConsequence;
    
    my $oc = Bio::EnsEMBL::Variation::OverlapConsequence->new(
        -display_term       => 'NON_SYNONYMOUS_CODING',
        -SO_term            => 'non_synonymous_codon',
        -SO_accession       => 'SO:0001583',
        -NCBI_term          => 'missense',
        -feature_SO_term    => 'mRNA',
        -description        => 'In coding sequence and results in an amino acid change in the encoded peptide sequence',
        -predicate          => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::non_synonymous_codon',
        -label              => 'Non-synonymous coding',
        -rank               => 7,
        -feature_class      => 'Bio::EnsEMBL::Transcript',
        -impact             => 'MODERATE',
        -include            => {
          coding => 1,
        },
    );

    if ($oc->predicate($transcript_variation_allele)) {
        print "This allele is: ", $oc->display_term, "\n";
    }

=head1 DESCRIPTION

An OverlapConsequence represents the consequence of an allele of a VariationFeature overlapping
some other Ensembl Feature (and therefore applies to VariationFeatureOverlapAllele objects as these
represent just such an event). It contains various values that represent the consequence type, such
as the Sequence Ontology (SO) term and accession (which should always be unique), the Ensembl
display_term (which will not always be unique), the relative rank of this consequence when compared 
to other consequences etc. It also contains a reference to a subroutine, referred to as the 
'predicate', which if a called with a VariationFeatureOverlapAllele (or a subclass) as the first and
only argument, will return a true or false value if this consequence type applies to this allele.

The list of OverlapConsequences used by Ensembl is defined in the Bio::EnsEMBL::Variation::Utils::Constants
module, and can be imported from there.

=cut

package Bio::EnsEMBL::Variation::OverlapConsequence;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Variation::Utils::VariationEffect;

=head2 new
    
  Arg [-SO_ACCESSION] : 
    The Sequence Ontology accession for this consequence type
  
  Arg [-SO_TERM] : 
    The Sequence Ontology term for this consequence type
  
  Arg [-FEATURE_SO_TERM] : 
    The Sequence Ontology term for the feature affected by this consequence type
  
  Arg [-FEATURE_CLASS] : 
    The Ensembl class that represents the feature affected by this consequence type

  Arg [-VARIANT_FEATURE_CLASS] : 
    The Ensembl class that represents the variation feature this consequence applies to

  Arg [-INCLUDE] :
    A hashref of predicate conditions
  
  Arg [-PREDICATE] : 
    A reference to a subroutine that checks if this consequence type holds for
    a given VariationFeatureOverlapAllele (or the name of such a subroutine)
  
  Arg [-RANK] : 
    The relative rank of this consequence type when compred to other OverlapConsequence
    objects
  
  Arg [-IMPACT] : 
    Impact rating, one of: MODIFIER, LOW, MODERATE, HIGH
  
  Arg [-DISPLAY_TERM] : 
    The Ensembl display term for this consequence type (used by default on the website)

  Arg [-NCBI_TERM] : 
    The NCBI term for this consequence type
  
  Arg [-DESCRIPTION] : 
    A freetext description of this consequence type (used on the website)
  
  Arg [-LABEL] : 
    A freetext label briefly describing this consequence type (used on the website)
  
  Arg [-IS_DEFAULT] : 
    A flag indicating if this is the default consequence type used when none other applies
    (in Ensembl this currently set on the intergenic OverlapConsequence)
  
  Example : 
    my $oc = Bio::EnsEMBL::Variation::OverlapConsequence->new(
        -display_term       => 'NON_SYNONYMOUS_CODING',
        -SO_term            => 'non_synonymous_codon',
        -SO_accession       => 'SO:0001583',
        -NCBI_term          => 'missense',
        -feature_SO_term    => 'mRNA',
        -description        => 'In coding sequence and results in an amino acid change in the encoded peptide sequence',
        -predicate          => 'Bio::EnsEMBL::Variation::Utils::VariationEffect::non_synonymous_codon',
        -label              => 'Non-synonymous coding',
        -rank               => 7,
        -impact             => 'MODERATE',
        -tier               => 1,
        -feature_class      => 'Bio::EnsEMBL::Transcript',
        -include            => {
          coding => 1,
        },
    );
  
  Description: Constructs a new OverlapConsequence instance
  Returntype : A new Bio::EnsEMBL::Variation::OverlapConsequence instance 
  Exceptions : none
  Status     : Stable

=cut 

sub new {
    my $class = shift;
    
    my (
        $SO_accession,
        $SO_term,
        $feature_SO_term,
        $feature_class,
        $variant_feature_class,
        $include,
        $predicate,
        $rank,
        $impact,
        $tier,
        $display_term,
        $NCBI_term,
        $description,
        $label,
        $is_default,
    ) = rearrange([qw(
            SO_ACCESSION
            SO_TERM
            FEATURE_SO_TERM
            FEATURE_CLASS
            VARIANT_FEATURE_CLASS
            INCLUDE
            PREDICATE
            RANK
            IMPACT
            TIER
            DISPLAY_TERM
            NCBI_TERM
            DESCRIPTION
            LABEL
            IS_DEFAULT
        )], @_);

    my $self = bless {
        SO_accession            => $SO_accession,
        SO_term                 => $SO_term,
        feature_SO_term         => $feature_SO_term,
        feature_class           => $feature_class,
        variant_feature_class   => $variant_feature_class,
        include                 => $include,
        predicate               => $predicate,
        rank                    => $rank,
        impact                  => $impact,
        tier                    => $tier,
        display_term            => $display_term,
        NCBI_term               => $NCBI_term,
        description             => $description,
        label                   => $label,
        is_default              => $is_default,
    }, $class;

    return $self;
}

sub new_fast {
    my ($class, $hashref) = @_;
    return bless $hashref, $class;
}

=head2 SO_accession

  Arg [1]    : (optional) accession to set
  Description: Get/set the Sequence Ontology accession for this consequence type
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut

sub SO_accession {
    my ($self, $SO_accession) = @_;
    $self->{SO_accession} = $SO_accession if $SO_accession;
    return $self->{SO_accession};
}

=head2 SO_term

  Arg [1]    : (optional) term to set
  Description: Get/set the Sequence Ontology term for this consequence type
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut

sub SO_term {
    my ($self, $SO_term) = @_;
    $self->{SO_term} = $SO_term if $SO_term;
    return $self->{SO_term};
}

=head2 feature_SO_term

  Arg [1]    : (optional) term to set
  Description: Get/set the Sequence Ontology term for the feature affected by this consequence type
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut

sub feature_SO_term {
    my ($self, $feature_SO_term) = @_;
    $self->{feature_SO_term} = $feature_SO_term if $feature_SO_term;
    return $self->{feature_SO_term};
}

=head2 feature_class

  Arg [1]    : (optional) class to set
  Description: Get/set the Ensembl class representing the feature affected by this consequence type
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut

sub feature_class {
    my ($self, $feature_class) = @_;
    $self->{feature_class} = $feature_class if $feature_class;
    return $self->{feature_class} || '';
}

=head2 include

  Arg [1]    : (optional) hashref of include conditions
  Description: Get/set include conditions for predicate to be executed
  Returntype : hashref
  Exceptions : none
  Status     : Stable

=cut

sub include {
  my ($self, $include) = @_;
  $self->{include} = $include if $include;
  return $self->{include} || {};
}

=head2 predicate

  Arg [1]    : (optional) reference to subroutine (or the name of a subroutine)
  Description: Get/set the predicate used to check if this consequence type applies 
               to a given VariationFeatureOverlapAllele. Currently, if you supply 
               a name (rather than a coderef), this subroutine must be found in the
               Bio::EnsEMBL::Variation::Utils::VariationEffect module.
  Returntype : coderef
  Exceptions : throws if a name is supplied and the subroutine cannot be found in 
               the expected module
  Status     : Stable

=cut
 
sub predicate {
    my ($self, $predicate) = @_;
    
    $self->{predicate} = $predicate if $predicate;
    
    if ($self->{predicate} && ref $self->{predicate} ne 'CODE') {
        my $name = $self->{predicate};

        if (defined &$name && $name =~ /^Bio::EnsEMBL::Variation::Utils::VariationEffect/) {
            $self->{predicate} = \&$name;
        }
        else {
            throw("Can't find a subroutine called $name in the VariationEffect module?");
        }
    }
    
    return $self->{predicate};
}

=head2 rank

  Arg [1]    : (optional) rank to set
  Description: Get/set the relative rank of this OverlapConsequence when compared to other
               OverlapConsequence objects. This is used, for example, to determine the most 
               severe consequence of a VariationFeature. 
  Returntype : integer
  Exceptions : none
  Status     : Stable

=cut

sub rank {
    my ($self, $rank) = @_;
    $self->{rank} = $rank if $rank;
    return $self->{rank};
}

=head2 impact

  Arg [1]    : (optional) impact level to set
  Description: Get/set the impact level of this OverlapConsequence. One of MODIFIER, LOW,
               MODERATE, HIGH 
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut

sub impact {
    my ($self, $impact) = @_;
    $self->{impact} = $impact if $impact;
    return $self->{impact};
}

=head2 tier

  Arg [1]    : (optional) tier to set
  Description: Get/set the tier this OverlapConsequence belongs to. Variations will be
               assigned consequences in tier order; if a tier 1 consequence is assigned,
               no tier 2 consequences will be checked/assigned. 
  Returntype : integer
  Exceptions : none
  Status     : At Risk

=cut

sub tier {
    my ($self, $tier) = @_;
    $self->{tier} = $tier if $tier;
    return $self->{tier};
}


=head2 display_term

  Arg [1]    : (optional) term to set
  Description: Get/set the Ensembl display term for this consequence type. This is
               used by default on the website.
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub display_term {
    my ($self, $display_term) = @_;
    $self->{display_term} = $display_term if $display_term;
    return $self->{display_term} || $self->SO_term;
}

=head2 NCBI_term

  Arg [1]    : (optional) term to set
  Description: Get/set the NCBI term for this consequence type
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub NCBI_term {
    my ($self, $NCBI_term) = @_;
    $self->{NCBI_term} = $NCBI_term if $NCBI_term;
    return $self->{NCBI_term} || $self->SO_term;
}

=head2 description

  Arg [1]    : (optional) description to set
  Description: Get/set the description for this consequence type. This is used on the
               website and is intended to be a freetext description of this consequence.
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut

sub description {
    my ($self, $description) = @_;
    $self->{description} = $description if $description;
    return $self->{description};
}

=head2 label

  Arg [1]    : (optional) label to set
  Description: Get/set the label for this consequence type. This is used on the
               website and is intended to be a short description of this consequence.
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut

sub label {
    my ($self, $label) = @_;
    $self->{label} = $label if $label;
    return $self->{label};
}

=head2 variant_feature_class

  Arg [1]    : (optional) class as a atring
  Description: Get/set the class of variant features that this consequence can apply to
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut

sub variant_feature_class {
    my ($self, $class) = @_;
    $self->{variant_feature_class} = $class if $class;
    return $self->{variant_feature_class};
}

=head2 is_default

  Arg [1]    : (optional) flag
  Description: Get/set a flag indicating if this is the default consequence type.
               There should only be one default OverlapConsequence, in Ensembl this
               flag is only set on the INTERGENIC OverlapConsequence object.
  Returntype : bool
  Exceptions : none
  Status     : Stable

=cut

sub is_default {
    my ($self, $is_default) = @_;
    $self->{is_default} = $is_default if defined $is_default;
    return $self->{is_default};
}

sub get_all_parent_SO_terms {
    my ($self) = @_;
    
    if (my $adap = $self->{adaptor}) {
        if (my $goa = $adap->db->get_SOTermAdaptor) {
            
        }
    }
}

1;

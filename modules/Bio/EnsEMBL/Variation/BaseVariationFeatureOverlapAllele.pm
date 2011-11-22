package Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);

sub new_fast {
    my ($class, $hashref) = @_;
    return bless $hashref, $class;
}

sub base_variation_feature_overlap {
    my ($self, $bvfo) = @_;
    $self->{base_variation_feature_overlap} = $bvfo if defined $bvfo;
    return $self->{base_variation_feature_overlap};
}

sub base_variation_feature {
    my $self = shift;
    return $self->base_variation_feature_overlap->base_variation_feature(@_);
}

sub feature {
    my $self = shift;
    return $self->base_variation_feature_overlap->feature(@_);
}

sub length {
    my ($self, $length) = @_;
    $self->{length} = $length if defined $length;
    return $self->{length};
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
            
            #warn "variant feature class: ".$oc->variant_feature_class."\n";
            #warn "base variation feature class: ".(ref $self->base_variation_feature)."\n";
            #warn "feature class: ".$oc->feature_class."\n\n";
            
            if ($oc->variant_feature_class && $self->base_variation_feature->isa($oc->variant_feature_class)) {
                if ($self->feature->isa($oc->feature_class)) {
                    if ($oc->predicate->($self)) {
                        push @$cons, $oc;
                    }
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

1;

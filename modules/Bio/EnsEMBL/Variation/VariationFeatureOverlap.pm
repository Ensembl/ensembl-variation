package Bio::EnsEMBL::Variation::VariationFeatureOverlap;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(expand);
use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

sub new {
    my ($class, $hashref) = @_;
    
    my $self = bless $hashref, $class;
    
    # now look at each allele of the VariationFeature in turn
    
    # get the allele string, expand it, and split it into separate alleles
    
    my $vf   = $self->{variation_feature};
    my $tran = $self->{transcript};
    
    my $allele_string = $vf->allele_string;
    
    expand(\$allele_string);
    
    my @alleles = split /\//, $allele_string;
  
    # create an object representing the reference allele
    
    my $ref_allele = shift @alleles;

    my $ref_vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
        variation_feature_overlap   => $self,
        variation_feature_seq       => $ref_allele,
        is_reference                => 1,
    });
    
    $self->reference_allele($ref_vfoa);

    # create objects representing the alternate alleles
    
    my @alt_alleles;
    
    for my $allele (@alleles) {
        
        my $vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
            variation_feature_overlap   => $self,
            variation_feature_seq       => $allele,
        });
        
        push @alt_alleles, $vfoa,
    }
    
    $self->alt_alleles(@alt_alleles);
    
    return $self;
}

sub new_fast {
    my ($class, $hashref) = @_;
    return bless $hashref, $class;
}

sub dbID {
    my ($self, $dbID) = @_;
    $self->{dbID} = $dbID if defined $dbID;
    return $self->{dbID};
}

sub feature_type_id {
    my ($self, $feature_type_id) = @_;
    $self->{feature_type_id} = $feature_type_id if defined $feature_type_id;
    # XXX: find the correct feature type id
    return $self->{feature_type_id} || 1;
}

sub variation_feature {
    my ($self, $variation_feature) = @_;
    
    $self->{variation_feature} = $variation_feature if $variation_feature;
    
    if (my $vf_id = $self->{_variation_feature_id}) {
        
        # lazy-load the VariationFeature
        
        if (my $adap = $self->{adaptor}) {
            if (my $vfa = $adap->db->get_VariationFeatureAdaptor) {
                if (my $vf = $vfa->fetch_by_dbID($vf_id)) {
                    $self->{variation_feature} = $vf;
                    delete $self->{_variation_feature_id};
                }
            }
        }
    }
    
    return $self->{variation_feature};
}

sub feature {
    my ($self, $feature) = @_;
    
    $self->{feature} = $feature if $feature;
    
    # we can't lazy-load the feature here because only sub-classes
    # know which class their feature is, so the lazy-loading has to
    # be implemented in the subclass method
 
    return $self->{feature};
}

sub reference_allele {
    my ($self, $reference_allele) = @_;
    $self->{reference_allele} = $reference_allele if $reference_allele;
    return $self->{reference_allele};
}

sub alt_alleles {
    my ($self, @new_alt_alleles) = @_;
    
    my $alt_alleles = $self->{alt_alleles} ||= [];

    push @$alt_alleles, @new_alt_alleles if @new_alt_alleles;

    return $alt_alleles;
}

sub alleles {
    my ($self) = @_;
    return [ $self->reference_allele, @{ $self->alt_alleles } ];
}

1;



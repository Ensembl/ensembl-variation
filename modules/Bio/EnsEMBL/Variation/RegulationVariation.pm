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

package Bio::EnsEMBL::Variation::RegulationVariation;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlap);

sub new {
    my $class = shift;
    return $class->SUPER::new(@_);
}

sub feature_label {
    my ($self, $feature_label) = @_;
    $self->{feature_label} = $feature_label if $feature_label;
    return $self->{feature_label};
}

sub target_feature {
    # XXX: fetch the target feature
}

sub target_feature_stable_id {
    my ($self, $target_feature_stable_id) = @_;
    $self->{target_feature_stable_id} = $target_feature_stable_id if $target_feature_stable_id;
    return $self->{target_feature_stable_id};
}

sub _fake_stable_id {
    my $self = shift;
    my $f = $self->feature;
    return $f->stable_id if $f->can('stable_id');
    return join '_', $f->display_label, $f->seq_region_start, $f->seq_region_end
}

1;

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

package Bio::EnsEMBL::Variation::Pipeline::InitRegulationEffect;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub fetch_input {
    my $self = shift;
    my $species = $self->param_required('species');
    $self->warning($species); 

    my $cdba = $self->get_adaptor($species, 'core');
    my $fdba = $self->get_adaptor($species, 'funcgen');
    my $vdba = $self->get_adaptor($species, 'variation');

    # clear tables
    my $vdbc = $vdba->dbc();

    unless ($self->param('only_regulatory_feature')) {
        $vdbc->do('TRUNCATE TABLE motif_feature_variation');
        $vdbc->do('ALTER TABLE motif_feature_variation DISABLE KEYS');
    }
    unless ($self->param('only_motif_feature')) {
#        $vdbc->do('TRUNCATE TABLE regulatory_feature_variation');
#        $vdbc->do('ALTER TABLE regulatory_feature_variation DISABLE KEYS');
    }

    # get regulation object ids
    my $rfa = $fdba->get_RegulatoryFeatureAdaptor or die 'Failed to get RegulatoryFeatureAdaptor';
    my $mfa = $fdba->get_MotifFeatureAdaptor or die 'Failed to get MotifFeatureAdaptor';
    my $fsa = $fdba->get_FeatureSetAdaptor or die 'Failed to get FeatureSetAdaptor';

    my $sa = $cdba->get_SliceAdaptor or die 'Failed to get SliceAdaptor';

    my $slices = $sa->fetch_all('toplevel', undef, 0, 1);

    if ($self->param('debug')) {
        my $slice = $sa->fetch_by_region('chromosome', 12);
        $slices = [];
        push @$slices, $slice;
    }

#    my $regulatory_feature_set = $fsa->fetch_by_name('RegulatoryFeatures:MultiCell');
#    my @external_feature_sets = @{$fsa->fetch_all_by_type('external')};

    foreach my $slice (@$slices) {
        # get all RegulatoryFeatures
        my @feature_ids = ();
        unless ($self->param('only_motif_feature')) {
#          my $it = $rfa->fetch_Iterator_by_Slice($slice);
          $self->warning($slice->seq_region_name);
          my @rfs = @{$rfa->fetch_all_by_Slice($slice) || []};
          $self->warning(scalar @rfs);
          foreach my $rf (@rfs) {
              push @feature_ids, { feature_id => $rf->stable_id,
                                   feature_type => 'regulatory_feature',
                                   species => $species, };
          }
        }
        # get all MotifFeatures
        unless ($self->param('only_regulatory_feature')) {
          $self->warning($slice->seq_region_name);

            my @mfs = @{$mfa->fetch_all_by_Slice($slice) || []};
            foreach my $mf (@mfs) {
                push @feature_ids, { feature_id => $mf->dbID,
                                     feature_type => 'motif_feature',
                                     species => $species, };  
            }
        }
        # get all ExternalFeatures
#        if ($self->param('include_external_features')) {
#            foreach my $external_fset (@external_feature_sets) {
#                my $feature_set = $fsa->fetch_by_name($external_fset->name);
#                foreach my $external_feature (@{$feature_set->get_Features_by_Slice($slice)}) {
#                    push @feature_ids, { feature_id => $external_feature->dbID,
#                                         feature_type => 'external_feature', };
#                }
#            }
#        }
        $self->dataflow_output_id(\@feature_ids, 2);
    }
}


sub write_output {
    my $self = shift;
    return;
}
1;

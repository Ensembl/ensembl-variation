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

package Bio::EnsEMBL::Variation::Pipeline::FinishRegulationEffect;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {
    my $self = shift;

    my $consequence_types = "'intergenic_variant','splice_acceptor_variant','splice_donor_variant','stop_lost','coding_sequence_variant','missense_variant','stop_gained','synonymous_variant','frameshift_variant','non_coding_transcript_variant','non_coding_transcript_exon_variant','mature_miRNA_variant','NMD_transcript_variant','5_prime_UTR_variant','3_prime_UTR_variant','incomplete_terminal_codon_variant','intron_variant','splice_region_variant','downstream_gene_variant','upstream_gene_variant','start_lost','stop_retained_variant','inframe_insertion','inframe_deletion','transcript_ablation','transcript_fusion','transcript_amplification','transcript_translocation','TFBS_ablation','TFBS_fusion','TFBS_amplification','TFBS_translocation','regulatory_region_ablation','regulatory_region_fusion','regulatory_region_amplification','regulatory_region_translocation','feature_elongation','feature_truncation','protein_altering_variant','regulatory_region_variant','TF_binding_site_variant'";

  if ($self->param('update_vf') || $self->param('only_update_vf')) {
    my $vdba = $self->get_species_adaptor('variation');
    my $dbc = $vdba->dbc;

    # pre-clean-up:
    foreach my $table (qw/regulatory_region_consequences variation_feature_overlap_regulation variation_feature_consequences/) {
      $dbc->do(qq{DROP TABLE IF EXISTS $table;})
    }

    my @regulatory_tables = ('motif_feature_variation', 'regulatory_feature_variation');

    $self->warning('Collect variation features overlapping regulatory features');

    $dbc->do(qq{
      CREATE TABLE IF NOT EXISTS regulatory_region_consequences(
      variation_feature_id int(10), 
      consequence_types set($consequence_types) NOT NULL DEFAULT 'intergenic_variant',
      key variation_feature_idx(variation_feature_id)
    );});

    for my $table (@regulatory_tables) {
      $dbc->do(qq{
        INSERT INTO regulatory_region_consequences(variation_feature_id, consequence_types)
        SELECT variation_feature_id, consequence_types 
        FROM $table;});
    }

    $self->warning('Completed collect variation features overlapping regulatory features.');
    $self->warning('Collect overlap with variation_feature table.');

    $dbc->do(qq{
      CREATE TABLE IF NOT EXISTS variation_feature_overlap_regulation(
      variation_feature_id int(10), 
      consequence_types set($consequence_types) NOT NULL DEFAULT 'intergenic_variant',
      key variation_feature_idx(variation_feature_id)
    );});

    $dbc->do(qq{
      INSERT INTO variation_feature_overlap_regulation(variation_feature_id, consequence_types)
      SELECT vf.variation_feature_id, vf.consequence_types
      FROM variation_feature vf, regulatory_region_consequences rrc
      WHERE rrc.variation_feature_id = vf.variation_feature_id;});

    $dbc->do(qq{
      INSERT INTO regulatory_region_consequences 
      SELECT * FROM variation_feature_overlap_regulation;});

    # combine consequences for all variation_features in regulatory_region_consequences
    my $tmp_table = 'variation_feature_consequences';

    $dbc->do(qq{CREATE TABLE IF NOT EXISTS $tmp_table LIKE variation_feature;});

    $dbc->do(qq{
     INSERT INTO $tmp_table (variation_feature_id, consequence_types)
     SELECT variation_feature_id, GROUP_CONCAT(DISTINCT(consequence_types))
     FROM regulatory_region_consequences
     GROUP BY variation_feature_id;}) or die "Populating tmp table failed";

    $self->warning('Final update of variation_feature table.');

    $tmp_table = 'variation_feature_consequences';
    # update variation feature
    $dbc->do(qq{
      UPDATE variation_feature vf, $tmp_table vfc
      SET vf.consequence_types = vfc.consequence_types
      WHERE vf.variation_feature_id = vfc.variation_feature_id;}) or die "Failed to update vf table";

    $self->warning('Completed update of variation_feature table.');

    # post-clean-up:
    foreach my $table (qw/regulatory_region_consequences variation_feature_overlap_regulation variation_feature_consequences/) {
      $dbc->do(qq{DROP TABLE IF EXISTS $table;})
    }
    foreach my $table (qw/motif_feature_variation regulatory_feature_variation/) {
      $dbc->do(qq{ALTER TABLE $table ENABLE KEYS;});
    } 

  }
}

1;

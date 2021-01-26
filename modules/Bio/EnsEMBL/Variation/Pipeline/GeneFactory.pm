=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::Pipeline::GeneFactory;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use File::Path qw(mkpath rmtree);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;
my $DEBUG_GENE_IDS = 0;
sub fetch_input {
   
    my $self = shift;

    my $include_lrg = $self->param('include_lrg');
    my $biotypes = $self->param('limit_biotypes');

    my $core_dba = $self->get_species_adaptor('core');
    my $var_dba = $self->get_species_adaptor('variation');
    
    my $dbc = $var_dba->dbc();

    my $ga = $core_dba->get_GeneAdaptor or die "Failed to get gene adaptor";

    my @genes;

    if ( grep {defined($_)} @$biotypes ) {  # If array is not empty  
       # Limiting genes to specified biotypes 
       @genes = map { @{$ga->fetch_all_by_logic_name($_)} } @$biotypes;

    } else { 
       # fetch all genes 
       @genes = @{ $ga->fetch_all };
    }  

    if ($include_lrg) {
        # fetch the LRG genes as well
        push @genes, @{ $ga->fetch_all_by_biotype('LRG_gene') }
    }

    my @gene_output_ids; 
    my $gene_count = 0;

    for my $gene (@genes) {
      $gene_count++;
      push @gene_output_ids, {
        gene_stable_id  => $gene->stable_id,
      };
      
      if ($DEBUG) {
          last if $gene_count >= 500;
      }
    }

    if ($DEBUG_GENE_IDS) {
      my @gene_ids = qw/ENSG00000078328 ENSG00000149972 ENSG00000182185/; # slow: ENSG00000078328, fast: ENSG00000276644
      foreach my $gene_id (@gene_ids) {
        if (!grep {$_->{gene_stable_id} eq $gene_id } @gene_output_ids) {
          push @gene_output_ids, {
            gene_stable_id => $gene_id,
          };
        }
      }
    } 


    $self->param('gene_output_ids', \@gene_output_ids);

}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id($self->param('gene_output_ids'), 2);
  return;
}

1;

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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
use FileHandle;
use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;
my $DEBUG_GENE_IDS = 0;
sub fetch_input {
   
    my $self = shift;

    my $mtmp = $self->param('mtmp_table');

    my $include_lrg = $self->param('include_lrg');
    my $biotypes = $self->param('limit_biotypes');

    my $core_dba = $self->get_species_adaptor('core');
    my $var_dba = $self->get_species_adaptor('variation');
    
    my $dbc = $var_dba->dbc;

    my $ga = $core_dba->get_GeneAdaptor or die "Failed to get gene adaptor";

    my @genes;
    my @gene_output_ids; 
    my $gene_count = 0;
    my @delete_transcripts = ();
    my @old_genes = ();
    my %genes_hash;
    my @all_vf;

    if (-e $self->param('update_diff')) {

      my $file = $self->param('update_diff');
      open (DIFF, $file) or die "Can't open file $file: $!";

      while (<DIFF>) {
        chomp;
        next if /^transcript_id/;
        my ($transcript_id, $status, $gene_id) = split(/\t/);
        $transcript_id = "\'$transcript_id\'";
        $genes_hash{$gene_id} .= $transcript_id . ",";
      }

      foreach my $gene (keys %genes_hash) {
        my $core_dba = $self->get_species_adaptor('core');
        my $ga = $core_dba->get_GeneAdaptor or die "Failed to get gene adaptor";
        if(!defined( $ga->fetch_by_stable_id($gene) ) ) {
          my $transcripts = $genes_hash{$gene};
          $transcripts =~s/,$//;
          my @vf_ids = $dbc->do(qq{
              SELECT DISTINCT(variation_feature_id) FROM  transcript_variation
              WHERE feature_stable_id IN ($transcripts);
          });
          push (@all_vf, @vf_ids);
        }
      }

      open (DIFF, $file) or die "Can't open file $file: $!";

      while (<DIFF>) {
        chomp;
        next if /^transcript_id/;
        my ($transcript_id, $status, $gene_id, $other_info) = split(/\t/);
        if ($status ne "deleted") {
          push @gene_output_ids, {gene_stable_id  => $gene_id,}
        }
        if ($status eq "deleted") {
          push @delete_transcripts, $transcript_id;
        }
        # Remove Deleted transcripts
        if (@delete_transcripts > 500){
          my $joined_ids = '"' . join('", "', @delete_transcripts) . '"';
           
          $dbc->do(qq{
                     DELETE FROM  transcript_variation
                     WHERE   feature_stable_id IN ($joined_ids);
          }) or die "Deleting stable ids failed";

          $dbc->do(qq{
                      DELETE FROM  MTMP_transcript_variation
                      WHERE   feature_stable_id IN ($joined_ids);
           });
           # Reset delete_transcripts list
           @delete_transcripts = ();
        }
      }
	      my $tvdel_fh = FileHandle->new();
        $tvdel_fh->open(">>" .$self->param('pipeline_dir'). "/del_log/deleted_transcripts.txt") or die "Cannot open dump file " . $!;
        $self->dump_deleted_tv(@all_vf, $tvdel_fh); 
        $tvdel_fh->close();

        $include_lrg = 0; #Switch off as tends to be set to 1 in setup
    }
    elsif ( grep {defined($_)} @$biotypes ) {  # If array is not empty  
       # Limiting genes to specified biotypes 
       @genes = map { @{$ga->fetch_all_by_logic_name($_)} } @$biotypes;
    } 

    else { 
       # fetch all genes 
       @genes = @{ $ga->fetch_all };
    }  

    if ($include_lrg) {
        # fetch the LRG genes as well
        push @genes, @{ $ga->fetch_all_by_biotype('LRG_gene') }
    }

    for my $gene (@genes) {
      $gene_count++;
      push @gene_output_ids, {gene_stable_id  => $gene->stable_id,};
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

    # Remove Deleted transcripts - IS THIS NEEDED?
    if (-e $self->param('update_diff')){
        my $joined_ids = '"' . join('", "', @delete_transcripts) . '"';
        return if $joined_ids eq "";
        $dbc->do(qq{
                  DELETE FROM  transcript_variation
                  WHERE   feature_stable_id IN ($joined_ids);
        }) or die "Deleting stable ids failed";

       $dbc->do(qq{
                      DELETE FROM  MTMP_transcript_variation
                      WHERE   feature_stable_id IN ($joined_ids);
            });
    }
}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id($self->param('gene_output_ids'), 2);
  return;
}

sub dump_deleted_tv {
  my ($self, @deleted_vf, $fh) = @_;
  if (@deleted_vf) {
    foreach(@deleted_vf) {
      print $fh $_,"\n";
    }
  }
}

1;

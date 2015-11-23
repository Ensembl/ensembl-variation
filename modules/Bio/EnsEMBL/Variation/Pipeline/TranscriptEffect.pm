=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::TranscriptVariation;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT overlap);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG   = 0;

sub run {
  my $self = shift;

  my $gene_id = $self->required_param('gene_stable_id'); 

  my $disambiguate_sn_alleles = 
    $self->param('disambiguate_single_nucleotide_alleles');

  my $mtmp = $self->param('mtmp_table');
  
  my $variations_to_include;
  
  # if (my $vars = $self->param('variations_to_include')) {
  #   # turn the list of variation names into a hash to speed up checking
  #   $variations_to_include = { map { $_ => 1 } @$vars };
  # }

  my $core_dba = $self->get_species_adaptor('core');
  my $var_dba = $self->get_species_adaptor('variation');
  
  my $ga = $core_dba->get_GeneAdaptor;
  my $sa = $core_dba->get_SliceAdaptor;
  
  my $tva = $var_dba->get_TranscriptVariationAdaptor;

  # we need to include failed variations
  $tva->db->include_failed_variations(1);

  if((my $fasta = $self->param('fasta')) && !$Bio::EnsEMBL::Slice::_fasta_redefined && !$Bio::EnsEMBL::Slice::fasta_db) {

    # we need to find the assembly version to tell it about PARs
    my ($highest_cs) = @{$core_dba->get_CoordSystemAdaptor->fetch_all()};
    my $assembly = $highest_cs->version();

    setup_fasta(-FASTA => $fasta, -ASSEMBLY => $assembly);
  }

  my $gene = $ga->fetch_by_stable_id($gene_id) 
    or die "failed to fetch gene for stable id: $gene_id";


  # create summary tables for web index building
  my $genelu_ins_sth = $var_dba->dbc->prepare(qq[
    insert ignore into variation_genename
    (variation_id, gene_name) values (?,?)
  ]);

  my $hgvslu_ins_sth = $var_dba->dbc->prepare(qq[
    insert ignore into variation_hgvs
    (variation_id, hgvs_name) values (?,?)
  ]);

  my $slice = $sa->fetch_by_gene_stable_id(
    $gene->stable_id, 
    MAX_DISTANCE_FROM_TRANSCRIPT
  ) or die "failed to get slice around gene: ".$gene->stable_id;
  
  # call seq here to help cache
  $slice->seq();

  $gene = $gene->transfer($slice);

  my @vfs = (
    @{ $slice->get_all_VariationFeatures },
    @{ $slice->get_all_somatic_VariationFeatures }
  );

  my @write_data;

  for my $transcript (@{ $gene->get_all_Transcripts }) {

    for my $vf(@vfs) {

      if (defined $variations_to_include) {
        next unless $variations_to_include->{$vf->variation_name};
      }

      next unless overlap($vf->start, $vf->end, $transcript->start - MAX_DISTANCE_FROM_TRANSCRIPT, $transcript->end + MAX_DISTANCE_FROM_TRANSCRIPT);

      my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
        -transcript     => $transcript,
        -variation_feature  => $vf,
        -adaptor      => $tva,
        -disambiguate_single_nucleotide_alleles => $disambiguate_sn_alleles,
        -no_transfer    => 1,
      );

      # if the variation has no effect on the transcript $tv will be undef

      if ($tv && ( scalar(@{ $tv->consequence_type }) > 0) ) {

        # store now or save to store later? Uncomment out the behaviour you want
        # save to store later uses more memory but means you don't have to sort human TV after the run
        
        # $tva->store($tv, $mtmp);
        
        push @write_data, @{$tva->_get_write_data($tv)};
      
        ## populate tables for website index building

        my $var_id = $vf->get_Variation_dbID();

        $genelu_ins_sth->execute( $var_id, $gene->display_xref->display_id ) 
          if defined $gene->display_xref();

        for my $allele (@{ $tv->get_all_alternate_TranscriptVariationAlleles }) {

          next unless defined $allele->hgvs_transcript();

          my $hgvs_transcript = (split/\:/, $allele->hgvs_transcript())[1];

          $hgvslu_ins_sth->execute( $var_id, $hgvs_transcript) if defined $hgvs_transcript;
          
          next unless defined $allele->hgvs_protein();

          my $hgvs_protein  = (split/\:/, $allele->hgvs_protein())[1];

          $hgvslu_ins_sth->execute( $var_id, $hgvs_protein) 
            if defined $hgvs_protein && $hgvs_protein =~/^p/; ## don't store synonymous
        }
      }                       
    }
  }

  if(@write_data) {
    $var_dba->dbc->do("LOCK TABLES transcript_variation WRITE, MTMP_transcript_variation WRITE");
    $tva->_store_write_data(\@write_data, 1);
    $tva->_store_mtmp_write_data($tva->_get_mtmp_write_data_from_tv_write_data(\@write_data), 1) if $mtmp;
    $var_dba->dbc->do("UNLOCK TABLES");
  }

  return;
}

sub write_output {
  my $self = shift;
}

1;

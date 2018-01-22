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

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::TranscriptHaplotypeAdaptor;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::CDSHaplotype;
use Bio::EnsEMBL::Variation::ProteinHaplotype;
use Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

our %HAPLOTYPE_TYPES = (
  'Transcript' => 1,
  'CDS' => 1,
  'Protein' => 1,
);

# stores a TranscriptHaplotype object
# sub store {
# 	my ($self, $th) = @_;
	
#   assert_ref($th, 'Bio::EnsEMBL::Variation::TranscriptHaplotype');
  
# 	my $dbh = $self->dbc->db_handle;
	
# 	my $sth = $dbh->prepare_cached(qq{
# 		INSERT INTO transcript_haplotype(
# 			transcript_stable_id,
# 			type,
#       md5,
#       sequence,
#       differences
# 		) VALUES (?,?,?,?,?)
# 	});
	
#   $sth->execute(
#     $th->transcript->stable_id(),
#     $th->type(),
#     $th->hex(),
#     $th->has_indel ? $th->seq : undef,
#     join(",", map {$_->{diff}} @{$th->get_all_diffs}),
#   );
  
#   my $th_id = $dbh->last_insert_id(undef, undef, 'transcript_haplotype', 'transcript_haplotype_id');
  
# 	$sth->finish;
  
#   # add individuals
#   $sth = $dbh->prepare_cached(qq{
#     INSERT INTO individual_transcript_haplotype(
#       transcript_haplotype_id,
#       individual_id,
#       count
#     ) VALUES (?,?,?)
#   });
  
#   # prefetch individual dbIDs
#   my $ia = $self->db->get_IndividualAdaptor();
#   my %map = map {$_->name() => $_->dbID()} @{$ia->fetch_all_by_name_list([keys %{$th->{individuals}}])};
  
#   foreach my $ind(keys %{$th->{individuals}}) {
#     $sth->execute(
#       $th_id,
#       $map{$ind},
#       $th->{individuals}->{$ind}
#     );
#   }
  
# 	$sth->finish;
	
# 	return 1;
# }


=head2 fetch_all_by_Transcript

  Arg[1]     : Bio::EnsEMBL::Transcript $transcript
  Arg[2]     : (optional) Bio::EnsEMBL::Variation::Sample $sample
               OR Bio::EnsEMBL::Variation::Population $population
  Example    : my $ths = $tha->fetch_all_by_Transcript($transcript)
  Description: Get TranscriptHaplotypes for a transcript, optionally limited
               to a sample or population
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Transcript {
  my $self = shift;
  return $self->_generic_fetch_all_by_Transcript('Transcript', @_);
}


=head2 fetch_all_CDSHaplotypes_by_Transcript

  Arg[1]     : Bio::EnsEMBL::Transcript $transcript
  Arg[2]     : (optional) Bio::EnsEMBL::Variation::Sample $sample
               OR Bio::EnsEMBL::Variation::Population $population
  Example    : my $ths = $tha->fetch_all_CDSHaplotypes_by_Transcript($transcript)
  Description: Get CDSHaplotypes for a transcript, optionally limited
               to a sample or population
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_CDSHaplotypes_by_Transcript {
  my $self = shift;
  return $self->_generic_fetch_all_by_Transcript('CDS', @_);
}


=head2 fetch_all_ProteinHaplotypes_by_Transcript

  Arg[1]     : Bio::EnsEMBL::Transcript $transcript
  Arg[2]     : (optional) Bio::EnsEMBL::Variation::Sample $sample
               OR Bio::EnsEMBL::Variation::Population $population
  Example    : my $ths = $tha->fetch_all_ProteinHaplotypes_by_Transcript($transcript)
  Description: Get ProteinHaplotypes for a transcript, optionally limited
               to a sample or population
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_ProteinHaplotypes_by_Transcript {
  my $self = shift;
  return $self->_generic_fetch_all_by_Transcript('Protein', @_);
}


=head2 get_TranscriptHaplotypeContainer_by_Transcript

  Arg[1]     : Bio::EnsEMBL::Transcript $transcript
  Example    : my $thc = $tha->get_TranscriptHaplotypeContainer_by_Transcript($transcript)
  Description: Get TranscriptHaplotypeContainer for a transcript
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_TranscriptHaplotypeContainer_by_Transcript {
  my ($self, $tr, $filters) = @_;

  assert_ref($tr, 'Bio::EnsEMBL::Transcript');
  assert_ref($filters, 'HASH') if $filters;

  # THC gets cached on a transcript, but we want to update it if filters used have changed
  if($filters || $tr->{_cached_th_filter_key}) {
    delete $tr->{_transcript_haplotype_container};
    delete $tr->{_cached_th_filter_key};
  }

  $tr->{_cached_th_filter_key} = $filters ? 1 : 0;
  
  # we cache the container on the transcript so we don't fetch it more than once
  if(!exists($tr->{_transcript_haplotype_container})) {

    # get VCF Collection Adaptor
    $self->db->use_vcf(1);
    
    my $vca = $self->db->get_VCFCollectionAdaptor();

    # doing this saves a lot of time!
    $_->_add_Populations_to_Samples for @{$vca->fetch_all};
    
    my @gts;

    # we don't want variants in introns
    my $tr_slice = $tr->slice;

    foreach my $exon(@{$tr->get_all_Exons}) {
      push @gts,
        map {@{$_->get_all_SampleGenotypeFeatures_by_Slice($exon->feature_Slice, undef, 1)}}
        @{$vca->fetch_all};
    }

    # transfer attached VFs
    # but only do it once per unique VF, otherwise we duplicate effort and get issues downstream
    my %transferred_vfs;
    foreach my $gt(@gts) {
      my $vf = $gt->variation_feature;
      my $vf_dbID = $vf->dbID;

      my $tr_vf = $transferred_vfs{$vf_dbID} ||= $vf->transfer($tr_slice);
      $gt->variation_feature($tr_vf);
    }
    
    $tr->{_transcript_haplotype_container} = Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer->new(
      -transcript => $tr,
      -genotypes  => \@gts,
      -samples    => [map {@{$_->get_all_Samples}} @{$vca->fetch_all}],
      -db         => $self->db,
      -filters    => $filters,
    );
  }

  return $tr->{_transcript_haplotype_container};
}

# sub fetch_all_by_transcript_stable_id {
#   my $self = shift;
#   my $id   = shift;
#   return $self->generic_fetch(qq{th.transcript_stable_id = "$id"});
# }

# generic fetch method
sub _generic_fetch_all_by_Transcript {
  my ($self, $type, $tr, $sample_or_population) = @_;
  
  # sanity checking
  throw("ERROR: Type must be one of ".join(",", sort keys %HAPLOTYPE_TYPES)) unless $type && $HAPLOTYPE_TYPES{$type};
  assert_ref($tr, 'Bio::EnsEMBL::Transcript');

  my $container = $self->get_TranscriptHaplotypeContainer_by_Transcript($tr);
  return [] unless $container;

  if($sample_or_population) {
    my $method = sprintf('get_all_%sHaplotypes_by_Sample', $type);

    if($sample_or_population->isa('Bio::EnsEMBL::Variation::Sample')) {
      return $container->$method($sample_or_population);
    }
    elsif($sample_or_population->isa('Bio::EnsEMBL::Variation::Population')) {
      
      my @samples;

      if($sample_or_population->adaptor && $sample_or_population->adaptor->db) {
        @samples = @{$sample_or_population->get_all_Samples};
      }
      else {
        foreach my $sample(@{$container->get_all_Samples}) {
          push @samples, $sample if grep {$_->name eq $sample_or_population->name} @{$sample->get_all_Populations};
        }
      }

      return [values %{{map {$_->_hex => $_} map {@{$container->$method($_)}} @samples}}];      
    }
    else {
      throw("ERROR: Expected Bio::EnsEMBL::Variation::Sample or Bio::EnsEMBL::Variation::Population, got ".ref($sample_or_population));
    }
  }
  else {
    my $method = sprintf('get_all_%sHaplotypes', $type);
    return $container->$method;
  }
}

# sub _tables{
#   my $self = shift;
  
#   return (['transcript_haplotype','th'], ['individual_transcript_haplotype', 'ith']);
# }

# sub _columns{
#   return qw(th.transcript_haplotype_id th.transcript_stable_id th.type th.md5 th.sequence th.differences ith.individual_id ith.count);
# }


# sub _left_join {
#   return (
#     ['individual_transcript_haplotype', 'ith.transcript_haplotype_id = th.transcript_haplotype_id'],
#     #['population_transcript_haplotype', 'pth.transcript_haplotype_id = th.transcript_haplotype_id']
#   );
# }

# sub _objs_from_sth{
#   my $self = shift;
#   my $sth = shift;
  
# 	my ($transcript_haplotype_id, $transcript_stable_id, $type, $md5, $sequence, $differences, $individual_id, $count);
	
# 	$sth->bind_columns(\$transcript_haplotype_id, \$transcript_stable_id, \$type, \$md5, \$sequence, \$differences, \$individual_id, \$count);
  
#   my $by_id = {};
#   my %individuals = ();
  
#   while($sth->fetch) {
    
#     # create CDSHaplotype
#     if($type eq 'cds') {
#       $by_id->{$transcript_haplotype_id} ||= Bio::EnsEMBL::Variation::CDSHaplotype->new(
#         -dbID => $transcript_haplotype_id,
#         -type => $type,
#         -hex => $md5,
#         -sequence => $sequence,
#         -diffs => [map {{diff => $_}} split(",", $differences || '')],
#       );
#     }
    
#     # create ProteinHaplotype
#     elsif($type eq 'protein') {
#       $by_id->{$transcript_haplotype_id} ||= Bio::EnsEMBL::Variation::ProteinHaplotype->new(
#         -dbID => $transcript_haplotype_id,
#         -type => $type,
#         -hex => $md5,
#         -sequence => $sequence,
#         -diffs => [map {{diff => $_}} split(",", $differences || '')],
#       );
#     }
    
#     else {
#       throw('Invalid TranscriptHaplotype type "'.$type.'"');
#     }
    
#     # add individuals
#     $by_id->{$transcript_haplotype_id}->{individuals}->{$individual_id} = $count;
#     $individuals{$individual_id} = 1;
#   }
  
#   # fetch individuals
#   my $ia = $self->db->get_IndividualAdaptor();
#   my %inds_by_id = map {$_->dbID => $_} @{$ia->fetch_all_by_dbID_list([keys %individuals])};
  
#   # convert the ids to names
#   foreach my $th_id(keys %$by_id) {
#     my %new_hash = map {$inds_by_id{$_}->name => $by_id->{$th_id}->{individuals}->{$_}} keys %{$by_id->{$th_id}->{individuals}};
#     $by_id->{$th_id}->{individuals} = \%new_hash;
#   }
  
#   return [values %$by_id];
# }

1;

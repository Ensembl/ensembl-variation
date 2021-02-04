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

package Bio::EnsEMBL::Variation::Pipeline::AncestralAlleles::PostProcessing;

use strict;
use warnings;

use FileHandle;
use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {
  my $self = shift;
  my $species_data = $self->param('post_process_species');
  foreach my $data (@$species_data) {
    my $species_name = $data->{species_name};
    my $species_dir = $data->{species_dir};
    $self->param('species', $species_name);
    my $vdba = $self->get_species_adaptor('variation');
    my $dbc = $vdba->dbc;
    $self->update_ancestral_alleles($species_dir, $dbc);
    if ($self->param('create_stats')) {
      $self->compare_previous_release_stats($species_dir, $dbc);
    }
    $self->update_meta($dbc);
  }
}

sub update_ancestral_alleles {
  my $self = shift;
  my $species_dir = shift;
  my $dbc = shift;

  opendir(my $dh, $species_dir) || die "Can't opendir $species_dir: $!";
  my @update_files = grep { $_ =~ /\.out$/ } readdir($dh);
  closedir $dh;

  foreach my $update_file (@update_files) {
    my $fh = FileHandle->new("$species_dir/$update_file", 'r');
    while (<$fh>) {
      chomp;
      $dbc->do($_) or die $!;
    }
    $fh->close;
  } 
}

sub compare_previous_release_stats {
  my $self = shift;
  my $species_dir = shift;
  my $dbc = shift;
  my $non_dbSNP_only = $self->param('non_dbSNP_only');
  my $sql = qq/SELECT ancestral_allele, COUNT(*) FROM variation_feature GROUP BY ancestral_allele;/;
  if ($non_dbSNP_only) {
    my $dbSNP_source_id = $self->get_source_id('dbSNP');
    $sql = qq/SELECT ancestral_allele, COUNT(*) FROM variation_feature WHERE source_id != $dbSNP_source_id GROUP BY ancestral_allele;/;
  }
  my $previous_release_stats = {};
  my $fh = FileHandle->new("$species_dir/previous_release_stats", 'r');
  while (<$fh>) {
    chomp;
    my ($allele, $count) = split/\s/;
    $previous_release_stats->{$allele} = $count;
  }
  $fh->close;
  my $new_alleles = {};
  my $ancestral_allele_counts = $dbc->sql_helper()->execute( -SQL => $sql);
  $fh = FileHandle->new("$species_dir/cmp_previous_release_stats", 'w');
  print $fh "allele old_count new_count new_count/old_count\n";
  foreach (sort {$b->[1] <=> $a->[1]} @$ancestral_allele_counts) {
    my $allele = $_->[0] || 'NULL';
    my $new_count = $_->[1];
    $new_alleles->{$allele} = $new_count;
    my $old_count = $previous_release_stats->{$allele};
    if ($old_count) {
      my $diff = $new_count / $old_count;
      print $fh "$allele $old_count $new_count $diff\n";
    } else {
      print $fh "$allele 0 $new_count 0\n";
    }
  }

  foreach my $allele (keys %$previous_release_stats) {
    if (!$new_alleles->{$allele}) {
      my $old_count = $previous_release_stats->{$allele};
      print $fh "$allele $old_count 0 0\n";
    }
  }

  $fh->close;
}

sub update_meta {
  my $self = shift;
  my $dbc = shift;
  my $update_meta_sth = $dbc->prepare(qq[ insert ignore into meta( meta_key, meta_value) values (?,?)]);
  $update_meta_sth->execute('AncestralAlleles_run_date', $self->run_date() );
}


1;

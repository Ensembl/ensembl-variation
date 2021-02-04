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

package Bio::EnsEMBL::Variation::Pipeline::AncestralAlleles::Init;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);
use Data::Dumper;
use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {
  my $self = shift;

  my $run_pipeline = {};
 
  my $compara_dir = $self->param('compara_dir'); 
  opendir(my $dh, $compara_dir) || die "Can't opendir $compara_dir: $!";
  my @ancestral_files = grep { $_ ne '.' and $_ ne '..' } readdir($dh);
  closedir $dh;

  die "There are no ancestral files in $compara_dir" if (scalar @ancestral_files == 0);

  my $registry = 'Bio::EnsEMBL::Registry';
  my $registry_file = $self->param('ensembl_registry');
  $registry->load_all($self->param('ensembl_registry'));
  $self->param('registry', $registry);

  my $vdbas = $registry->get_all_DBAdaptors(-group => 'variation');
  my $species = {};
  foreach my $vdba (@$vdbas) {
      my $species_name = $vdba->species();
      my $cdba =  $registry->get_DBAdaptor($species_name, 'core');
      if ($cdba) {
        $self->param('species', $species_name);
        my $assembly = $self->get_assembly;
        my ($ancestral_file) = grep { $_ =~ /$species_name/ && $_ =~ /$assembly/ } @ancestral_files;
        if ($ancestral_file) {
          $run_pipeline->{$species_name} = { file => $ancestral_file, assembly => $assembly};
        }
      } else {
        $self->warning("No core database for $species_name");
      }
  }

  my @input = ();
  my @post_processing_input = ();
  foreach my $species_name (keys %$run_pipeline) {
    my $species_dir = $self->create_species_dir($species_name);
    if ($self->param('create_stats')) {
      $self->store_previous_release_stats($species_name, $species_dir, $self->param('non_dbSNP_only'));
    }
    my $batches = $self->get_batches($species_name, $self->param('non_dbSNP_only'));
    push @post_processing_input, {
      species_name => $species_name,
      species_dir => $species_dir, 
    };
    foreach my $batch (@$batches) {
      push @input, {
        batch => $batch,
        species_name => $species_name,
        species_dir => $species_dir,
        ancestral_file => $run_pipeline->{$species_name}->{file},
      };  
    }
  }
  $self->param('input', \@input);
  $self->param('post_processing_input', {post_process_species => \@post_processing_input});
}

sub create_species_dir {
  my $self = shift;
  my $species_name = shift;
  my $pipeline_dir = $self->param('pipeline_dir');
 
  unless (-d "$pipeline_dir/$species_name") {
    my $err;
    make_path("$pipeline_dir/$species_name", {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }
  return "$pipeline_dir/$species_name"; 
}

sub store_previous_release_stats {
  my $self = shift;
  my $species_name = shift;
  my $species_dir = shift;
  my $non_dbSNP_only = shift;
  my $sql = qq/SELECT ancestral_allele, COUNT(*) FROM variation_feature GROUP BY ancestral_allele;/;
  if ($non_dbSNP_only) {
    my $dbSNP_source_id = $self->get_source_id('dbSNP');
    $sql = qq/SELECT ancestral_allele, COUNT(*) FROM variation_feature WHERE source_id != $dbSNP_source_id GROUP BY ancestral_allele;/;
  }

  if (! -e "$species_dir/previous_release_stats") {
    my $vdba = $self->param('registry')->get_DBAdaptor($species_name, 'variation');
    my $registry = $self->param('registry');
    my $dbc = $vdba->dbc;
    my $ancestral_allele_counts = $dbc->sql_helper()->execute( -SQL => $sql);
    my $fh = FileHandle->new("$species_dir/previous_release_stats", 'w');
    foreach (sort {$b->[1] <=> $a->[1]} @$ancestral_allele_counts) {
      my $allele = $_->[0] || 'NULL';
      my $count = $_->[1];
      print $fh "$allele $count\n"; 
    }
    $fh->close;
  } 
}

sub get_batches {
  my $self = shift;
  my $species_name = shift;
  my $non_dbSNP_only = shift;
  my $batch_size = $self->param('batch_size');
  my @batches = ();
  my ($vf_start, $vf_end);
  if ($non_dbSNP_only) {
    ($vf_start, $vf_end) = @{$self->get_variation_feature_start_end_for_non_dbSNP($species_name)};
  } else {
    $vf_start = 1;
    $vf_end = $self->get_variation_feature_count($species_name);
  }
  my $end = $vf_start + $batch_size; 
  my $id = 1;
  while ($end < $vf_end) {
    push @batches, {start => $vf_start, end => $end, id => $id};
    $vf_start = $end + 1;
    $end = $end + $batch_size;
    $id++;
  }
  push @batches, {start => $vf_start, end => $vf_end, id => $id};
  return \@batches;
}

sub get_variation_feature_start_end_for_non_dbSNP {
  my $self = shift;
  my $species_name = shift;
  my $vdba = $self->param('registry')->get_DBAdaptor($species_name, 'variation');
  my $registry = $self->param('registry');
  my $dbc = $vdba->dbc;
  my $start_end = $dbc->sql_helper()->execute( -SQL =>qq/select min(variation_feature_id), max(variation_feature_id) from variation_feature where source_id != 1;/);
  return $start_end->[0];
}

sub get_variation_feature_count {
  my $self = shift;
  my $species_name = shift;
  my $pipeline_dir = $self->param('pipeline_dir'); 
  my $vf_count;
  my $vf_count_file = "$pipeline_dir/$species_name/vf_count";
  if (! -e $vf_count_file) {
    my $vdba = $self->param('registry')->get_DBAdaptor($species_name, 'variation');
    my $registry = $self->param('registry');
    my $dbc = $vdba->dbc;
    $vf_count = $dbc->sql_helper()->execute_single_result( -SQL =>qq/select max(variation_feature_id) from variation_feature;/);
    my $fh = FileHandle->new($vf_count_file, 'w');
    print $fh "$vf_count\n";
    $fh->close();
  } else {
    my $fh = FileHandle->new($vf_count_file, 'r');
    while (<$fh>) {
      chomp;
      $vf_count = $_;  
    } 
    $fh->close;
  }
  return $vf_count;
}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id($self->param('input'), 2);
  $self->dataflow_output_id($self->param('post_processing_input'), 1);
}

sub get_assembly {
  my $self = shift;
  my $core_dba = $self->get_species_adaptor('core');
  my $dbc = $core_dba->dbc;
  my $current_db_name = $dbc->dbname();
  my $species = $self->param('species');
  my $species_id = $self->get_species_id($dbc, $current_db_name, $species);
  my $sth = $dbc->prepare("SELECT version FROM ".$current_db_name.".coord_system WHERE species_id = ".$species_id." ORDER BY rank LIMIT 1;");
  $sth->execute();
  my $assembly;
  $sth->bind_columns(\$assembly);
  $sth->execute();
  $sth->fetch();
  $sth->finish();
  return $assembly;
}

sub get_species_id {
  my ($self, $dbc, $current_db_name, $species) = @_;
  my $species_id = $dbc->sql_helper()->execute_simple( -SQL =>qq/select species_id from $current_db_name.meta where meta_key = 'species.production_name' and meta_value ='$species';/);
  return $species_id->[0];
}

1;

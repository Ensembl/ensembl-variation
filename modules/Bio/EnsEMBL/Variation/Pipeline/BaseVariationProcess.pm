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

package Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess;

use strict;
use warnings;

use POSIX;

use Bio::EnsEMBL::Variation::Pipeline::TranscriptFileAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Hive::AnalysisJob;
use Digest::MD5 qw(md5_hex);

use base qw(Bio::EnsEMBL::Hive::Process);

sub param {
    my $self = shift;
    
    unless ($self->input_job) {
        # if we don't have an input job, add a dummy one (used when we're not 
        # running as part of a pipeline proper)
        $self->input_job(Bio::EnsEMBL::Hive::AnalysisJob->new);
    }

    return $self->SUPER::param(@_);
}

sub required_param {
    my $self        = shift;
    my $param_name  = shift;
   
    my $param_value = $self->param($param_name, @_);

    die "$param_name is a required parameter" unless defined $param_value;
    
    return $param_value;
}

sub get_transcript_file_adaptor {
    my $self = shift;
    my $transcripts = shift;
    
    unless ($self->{tfa}) {
        $self->{tfa} = Bio::EnsEMBL::Variation::Pipeline::TranscriptFileAdaptor->new(
            fasta_file  => $self->param('fasta_file'),
            transcripts => $transcripts,
        );
    }

    return $self->{tfa};
}

sub get_species_adaptor {
    my ($self, $group) = @_;

    my $species = $self->required_param('species');
    
    return $self->get_adaptor($species, $group);
}

sub get_adaptor {
    my ($self, $species, $group) = @_;
   
    my $dba;

    eval {
        $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
    };

    unless (defined $dba) {
        $self->_load_registry();
        $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
    }
    
    unless (defined $dba) {
        die "Failed to a get DBA for $species and group $group";
    }
    return $dba;
}

sub _load_registry {
    my ($self) = @_;
    
    my $reg_file = $self->required_param('ensembl_registry');

    die("ERROR: Registry file $reg_file not found\n") unless -e $reg_file;
    
    Bio::EnsEMBL::Registry->load_all($reg_file, 0, 1);
    
    return;
}

## use same date format for all pipelines 
sub run_date{
    my ($self) = @_;

    return strftime("%Y-%m-%d", localtime);
}

sub get_files_dir {
  my $self = shift;
  my $id = shift;
  my $folder_name = shift;

  my $dir = $self->required_param('pipeline_dir').'/'.$folder_name.'_files/'.substr(md5_hex($id), 0, 2);

  unless(-d $dir) {
    mkdir($dir) or die "ERROR: Could not make directory $dir\n";
  }

  return $dir;
}

sub run_cmd {
  my $self = shift;
  my $cmd = shift;
  if (my $return_value = system($cmd)) {
    $return_value >>= 8;
    die "system($cmd) failed: $return_value";
  }
}

sub get_assembly {
  my $self = shift;
  my $core_dba = $self->get_species_adaptor('core');
  my $genome_container_adaptor = $core_dba->get_GenomeContainerAdaptor;
  return $genome_container_adaptor->get_version();
}

sub get_schema_version {
  my $self = shift;
  my $core_dba = $self->get_species_adaptor('core');
  my $meta_container_adaptor = $core_dba->get_MetaContainerAdaptor;
  return $meta_container_adaptor->get_schema_version();
}

sub get_species_id {
  my $self = shift;
  my $core_dba = $self->get_species_adaptor('core');
  my $meta_container_adaptor = $core_dba->get_MetaContainerAdaptor;
  return $meta_container_adaptor->species_id();
}

sub get_source_id {
  my $self = shift;
  my $source_name = shift;
  my $variation_dba = $self->get_species_adaptor('variation');
  my $source_adaptor = $variation_dba->get_SourceAdaptor; 
  my $source = $source_adaptor->fetch_by_name($source_name);
  if (!defined $source) {
    die "Could not fetch source for name: $source_name" 
  }
  return $source->dbID();
} 

1;

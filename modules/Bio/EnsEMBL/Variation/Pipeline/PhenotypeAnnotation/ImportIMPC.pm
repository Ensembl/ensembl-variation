=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

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


=head1 ImportIMPC

This module imports Mouse IMPC data. The module fetched the data from
http://www.mousephenotype.org.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportIMPC;

use warnings;
use strict;

use File::Path qw(make_path);
use File::stat;
use POSIX 'strftime';
use LWP::Simple;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::MouseBasePhenotypeAnnotation');

my %source_info;
my $workdir;
my ($logFH, $errFH);

my $core_dba;
my $variation_dba;

my $debug;

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $coord_file   = $self->required_param('coord_file');

  $core_dba    = $self->get_species_adaptor('core');
  $variation_dba  = $self->get_species_adaptor('variation');

  $debug        = $self->param('debug_mode');
  $self->SUPER::set_debug($self->param('debug_mode'));

  %source_info = (source_description => 'International Mouse Phenotyping Consortium',
                  source_url => 'http://www.mousephenotype.org/',
                  object_type => 'Gene',
                  source_status => 'germline',
                  source_name => 'IMPC',        #source name in the variation db
                  source_name_short => 'IMPC',  #source identifier in the pipeline
                  #source version is set based on the EBI fetched data hash response_date
                  );

  $workdir = $pipeline_dir."/".$source_info{source_name_short}."/".$species;
  make_path($workdir);

  open ($logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species);
  open ($errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species);
  $self->SUPER::set_logFH($logFH);
  $self->SUPER::set_errFH($errFH);

  my $url ='/mi/impc/solr/genotype-phenotype';
  my $file_impc = "IMPC_phenotypes.txt";
  print $logFH "Found file (".$workdir."/".$file_impc."), will skip new fetch\n" if -e $workdir."/".$file_impc;
  $file_impc = $self->get_mouse_phenotype_data($workdir, $source_info{source_name}, $url) unless -e $workdir."/".$file_impc;

  $self->param('impc_file', $file_impc);
}

sub run {
  my $self = shift;

  my $file_impc = $self->required_param('impc_file');
  my $coord_file = $self->required_param('coord_file');

  #get source ids
  my $mouse_phenotype_source_ids = $self->get_mouse_phenotype_data_source_ids($workdir."/".$file_impc, $source_info{source_name}, $variation_dba);
  print $logFH "$source_info{source_name} source_id is ", $mouse_phenotype_source_ids->{$source_info{source_name}}, "\n" if ($debug);

  $source_info{source_version} = $self->update_mouse_phenotype_data_version($mouse_phenotype_source_ids, $variation_dba); 
  $self->clear_mouse_phenotype_data_from_last_release($mouse_phenotype_source_ids, $variation_dba);
  my $marker_coords = $self->get_marker_coords($workdir."/".$file_impc, $coord_file, $core_dba);

  # get phenotype data
  my $results = $self->parse_mouse_phenotype_data($workdir."/".$file_impc, $marker_coords, $source_info{source_name}, $mouse_phenotype_source_ids, $variation_dba);
  print $logFH "Got ".(scalar @{$results->{'phenotypes'}})." new phenotypes \n" if $debug ;

  #get source id
  my $source_id = $self->get_or_add_source(\%source_info,$variation_dba);
  print $logFH "$source_info{source_name} source_id is $source_id\n" if ($debug);

  # save phenotypes
  $self->set_skip_synonyms(1);
  $self->save_phenotypes(\%source_info, $results, $core_dba, $variation_dba);

  my %param_source = (source_name => $source_info{source_name},
                      type => $source_info{object_type});
  $self->param('output_ids', { source => \%param_source,
                              species => $self->required_param('species'),
                              workdir => $workdir
                            });
  close($logFH);
  close($errFH);
}

sub write_output {
  my $self = shift;

  if ($self->param('debug_mode')) {
    open (my $logPipeFH, ">", $workdir."/".'log_import_debug_pipe');
    print $logPipeFH "Passing $source_info{source_name} import (".$self->required_param('species').") for checks (check_phenotypes)\n";
    close ($logPipeFH);
  }
  $self->dataflow_output_id($self->param('output_ids'), 1);

}

1;

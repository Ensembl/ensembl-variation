=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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


=head1 ImportMGI

This module imports Mouse MGI data. The module fetched the data from
https://www.mousephenotype.org.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportMGI;

use warnings;
use strict;

use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::MouseBasePhenotypeAnnotation');

my %source_info;

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $coord_file   = $self->required_param('coord_file');

  $self->debug($self->param('debug_mode'));

  %source_info = (source_description => 'Mouse Genome Informatics',
                  source_url => 'http://www.informatics.jax.org/',
                  object_type => 'Gene',
                  source_status => 'germline',
                  source_name => 'MGI',        #source name in the variation db
                  source_name_short => 'MGI',  #source identifier in the pipeline
                  data_types => 'phenotype_feature',
                  #source version is set based on EBI fetch data hash response_date
                  );
  my $workdir = $pipeline_dir."/".$source_info{source_name_short};
  make_path($workdir) or die "Failed to create $workdir $!\n";
  $self->workdir($workdir);

  open(my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  #fetch source version at runtime:
  $source_info{source_version} = $self->fetch_data_version;

  #get input file MGI:
  #Note: if file already exists, new fetch will be skipped, this can lead to
  #discrepancies between file and version if a new release is out between runs
  my $url = '/mi/impc/solr/mgi-phenotype';
  my $file_mgi = "MGI_phenotypes.txt";
  print $logFH "Found file (".$workdir."/".$file_mgi."), will skip new fetch\n" if -e $workdir."/".$file_mgi;
  $file_mgi = $self->get_input_file($workdir, $source_info{source_name}, $url) unless -e $workdir."/".$file_mgi;

  $self->param('mgi_file', $file_mgi);
}

sub run {
  my $self = shift;

  my $file_mgi = $self->required_param('mgi_file');
  my $coord_file = $self->required_param('coord_file');

  # dump and clean pre-existing phenotypes
  $self->dump_phenotypes($source_info{source_name}, 1);

  #get source id
  my $source_id = $self->get_or_add_source(\%source_info);
  $self->print_logFH("$source_info{source_name} source_id is $source_id\n") if ($self->debug);

  my $marker_coords = $self->get_marker_coords($self->workdir."/".$file_mgi, $coord_file);

  # get phenotype data
  my $results = $self->parse_input_file($self->workdir."/".$file_mgi, $marker_coords, $source_info{source_name}, $source_id);
  $self->print_logFH("Got ".(scalar @{$results->{'phenotypes'}})." new phenotypes \n") if ($self->debug);

  # save phenotypes
  $self->skip_synonyms(1);
  $self->save_phenotypes(\%source_info, $results);

  my %param_source = (source_name => $source_info{source_name},
                      type => $source_info{object_type});
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species'),
                               workdir => $self->workdir,
                               run_type => $self->required_param('run_type')
                             });
}

sub write_output {
  my $self = shift;

  $self->print_pipelogFH("Passing $source_info{source_name_short} import (".$self->required_param('species').") for checks (check_phenotypes)\n") if ($self->debug);
  close($self->logFH) if defined $self->logFH ;
  close($self->errFH) if defined $self->errFH ;
  close($self->pipelogFH) if defined $self->pipelogFH ;

  $self->dataflow_output_id($self->param('output_ids'), 2);

}

1;

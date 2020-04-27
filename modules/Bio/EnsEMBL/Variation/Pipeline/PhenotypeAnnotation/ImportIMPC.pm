=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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
https://www.mousephenotype.org.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportIMPC;

use warnings;
use strict;

use File::Path qw(make_path);
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::MouseBasePhenotypeAnnotation');

my %source_info;

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $coord_file   = $self->required_param('coord_file');

  $self->debug($self->param('debug_mode'));

  %source_info = (source_description => 'International Mouse Phenotyping Consortium',
                  source_url => 'https://www.mousephenotype.org/',
                  object_type => 'Gene',
                  source_status => 'germline',
                  source_name => 'IMPC',        #source name in the variation db
                  source_name_short => 'IMPC',  #source identifier in the pipeline
                  data_types => 'phenotype_feature',
                  #source version is set based on the EBI fetched data hash response_date
                  );

  my $workdir = $pipeline_dir."/".$source_info{source_name_short};
  unless (-d $workdir) {
    my $err;
    make_path($workdir, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }
  $self->workdir($workdir);

  open(my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  #fetch source version at runtime:
  $source_info{source_version} = $self->fetch_data_version;

  #get input file IMPC:
  #Note: if the file already exists, a new file fetch will be skipped, this can lead to
  #discrepancies between file and version if a new release is out between runs
  my $url ='/mi/impc/solr/genotype-phenotype';
  my $file_impc = "IMPC_phenotypes.txt";
  print $logFH "Found file (".$workdir."/".$file_impc."), will skip new fetch\n" if -e $workdir."/".$file_impc;
  $file_impc = $self->get_input_file($workdir, $source_info{source_name}, $url) unless -e $workdir."/".$file_impc;

  $self->param('impc_file', $file_impc);
}

sub run {
  my $self = shift;

  my $file_impc = $self->required_param('impc_file');
  my $coord_file = $self->required_param('coord_file');

  # dump and clean pre-existing phenotypes
  $self->dump_phenotypes($source_info{source_name}, 1);

  #get source id and update with latest source_info data
  my $source_id = $self->get_or_add_source(\%source_info);
  $self->print_logFH("$source_info{source_name} source_id is $source_id\n") if ($self->debug);

  my $marker_coords = $self->get_marker_coords($self->workdir."/".$file_impc, $coord_file);

  # get phenotype data
  my $results = $self->parse_input_file($self->workdir."/".$file_impc, $marker_coords, $source_info{source_name}, $source_id);
  $self->print_logFH("Got ".(scalar @{$results->{'phenotypes'}})." new phenotypes \n") if ($self->debug);

  # save phenotypes
  $self->skip_synonyms(1);
  $self->save_phenotypes(\%source_info, $results);

  my %param_source = (source_name => $source_info{source_name},
                      type => $source_info{object_type});
  $self->param('output_ids', { source => \%param_source,
                              species => $self->required_param('species'),
                              workdir => $self->workdir
                            });
}

sub write_output {
  my $self = shift;

  $self->print_pipelogFH("Passing $source_info{source_name_short} import (".$self->required_param('species').") for checks (check_phenotypes)\n") if ($self->debug);
  close($self->logFH) if defined $self->logFH ;
  close($self->errFH) if defined $self->errFH ;
  close($self->pipelogFH) if defined $self->pipelogFH ;

  $self->dataflow_output_id($self->param('output_ids'), 1);

}

1;

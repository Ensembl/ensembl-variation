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


=head1 ImportMGI

This module imports Mouse MGI data. The module fetched the data from
http://www.mousephenotype.org.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportMGI;

use warnings;
use strict;

use File::Path qw(make_path);
use File::stat;
use POSIX 'strftime';
use LWP::Simple;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::MouseBasePhenotypeAnnotation');

my %source_info;
my $workdir;

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $coord_file   = $self->required_param('coord_file');

  $self->debug($self->param('debug_mode'));
  $self->core_db_adaptor($self->get_species_adaptor('core'));
  $self->variation_db_adaptor($self->get_species_adaptor('variation'));

  %source_info = (source_description => 'Mouse Genome Informatics',
                  source_url => 'http://www.informatics.jax.org/',
                  object_type => 'Gene',
                  source_status => 'germline',
                  source_name => 'MGI',        #source name in the variation db
                  source_name_short => 'MGI',  #source identifier in the pipeline
                  #source version is set based on EBI fetch data hash response_date
                  );
  $workdir = $pipeline_dir."/".$source_info{source_name_short}."/".$species;
  make_path($workdir);

  open (my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open (my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open (my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  #get input file MGI:
  my $url = '/mi/impc/solr/mgi-phenotype';
  my $file_mgi = "MGI_phenotypes.txt";
  print $logFH "Found file (".$workdir."/".$file_mgi."), will skip new fetch\n" if -e $workdir."/".$file_mgi;
  $file_mgi = $self->get_mouse_phenotype_data($workdir, $source_info{source_name}, $url) unless -e $workdir."/".$file_mgi;

  $self->param('mgi_file', $file_mgi);
}

sub run {
  my $self = shift;

  my $file_mgi = $self->required_param('mgi_file');
  my $coord_file = $self->required_param('coord_file');

  #get source ids
  my $mouse_phenotype_source_ids = $self->get_mouse_phenotype_data_source_ids($workdir."/".$file_mgi, $source_info{source_name});
  $self->print_logFH($source_info{source_name}, " source_id is ", $mouse_phenotype_source_ids->{$source_info{source_name}}, "\n") if ($self->debug);

  $source_info{source_version} = $self->update_mouse_phenotype_data_version($mouse_phenotype_source_ids);
  $self->clear_mouse_phenotype_data_from_last_release($mouse_phenotype_source_ids);
  my $marker_coords = $self->get_marker_coords($workdir."/".$file_mgi, $coord_file);

  # get phenotype data
  my $results = $self->parse_mouse_phenotype_data($workdir."/".$file_mgi, $marker_coords, $source_info{source_name}, $mouse_phenotype_source_ids);
  $self->print_logFH("Got ".(scalar @{$results->{'phenotypes'}})." new phenotypes \n") if ($self->debug);

  #get source id
  my $source_id = $self->get_or_add_source(\%source_info);
  $self->print_logFH("$source_info{source_name} source_id is $source_id\n") if ($self->debug);

  # save phenotypes
  $self->skip_synonyms(1);
  $self->save_phenotypes(\%source_info, $results);

  my %param_source = (source_name => $source_info{source_name},
                      type => $source_info{object_type});
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species'),
                               workdir => $workdir
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

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

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportMGI;

use warnings;
use strict;

use File::Path qw(make_path);
use File::stat;
use POSIX 'strftime';
use LWP::Simple;
use Data::Dumper; #TODO: remove if not needed
use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation qw($variation_dba);

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::MouseBasePhenotypeAnnotation');

my %source_info;
my $workdir;
my $core_dba;
my $variation_dba;
my $phenotype_dba;

my $debug;

sub fetch_input {
    my $self = shift;

    my $pipeline_dir = $self->required_param('pipeline_dir');
    my $species      = $self->required_param('species');
    my $coord_file   = $self->required_param('coord_file');

    $core_dba    = $self->get_species_adaptor('core');
    $variation_dba  = $self->get_species_adaptor('variation'); #TODO: why did the init -> Base class -> this not work?
    $phenotype_dba  = $variation_dba->get_PhenotypeAdaptor; 

    $debug        = $self->param('debug_mode');

    #TODO: 
    # original call: perl ensembl-variation/scripts/import/import_phenotype_data.pl -host ${host} -user ${user} -pass ${pass}  -dbname ${dbname} \
    #-source mgi -infile ${datadir}/mgifile?.txt -verbose -version 20121031

    %source_info = (source_description => 'Mouse Genome Informatics',
                    source_url => 'http://www.informatics.jax.org/',
                    object_type => 'Gene',

                    source_name => 'MGI', #TODO: figure out where this is used
                    source_status => 'germline',
                    source => 'mgi',
                    #source version is set based on ebi fetch data hash response_date
                    );
    
    $workdir = $pipeline_dir."/".$source_info{source_name}."/".$species;
    make_path($workdir);

    #get input file MGI:
    my $url = '/mi/impc/solr/mgi-phenotype';
    my $file_mgi = "mgi_phenotypes.txt";
    print "Found file (".$workdir."/".$file_mgi."), will skip new fetch\n" if -e $workdir."/".$file_mgi;
    $file_mgi = $self->get_mouse_phenotype_data($workdir, $source_info{source}, $url) unless -e $workdir."/".$file_mgi;

    $self->param('mgi_file', $file_mgi);
}

sub run {
  my $self = shift;

  my $file_mgi = $self->required_param('mgi_file');
  my $coord_file = $self->required_param('coord_file');

  local (*STDOUT, *STDERR);
  open STDOUT, ">>", $workdir."/".'log_import_out_'.$file_mgi; #TODO: what is best error/out log naming convention?
  open STDERR, ">>", $workdir."/".'log_import_err_'.$file_mgi;

  #get source ids
  my $mouse_phenotype_source_ids = $self->get_mouse_phenotype_data_source_ids($workdir."/".$file_mgi, $source_info{source},$variation_dba);
  print STDOUT $source_info{source}, " source_id is ", $mouse_phenotype_source_ids->{$source_info{source_name}}, "\n" if ($debug);

  $source_info{source_version} = $self->update_mouse_phenotype_data_version($mouse_phenotype_source_ids,$variation_dba); 
  $self->clear_mouse_phenotype_data_from_last_release($mouse_phenotype_source_ids,$variation_dba);
  my $marker_coords = $self->get_marker_coords($workdir."/".$file_mgi, $coord_file,$core_dba);

  # get phenotype data
  my $results = $self->parse_mouse_phenotype_data($workdir."/".$file_mgi, $marker_coords, $source_info{source}, $mouse_phenotype_source_ids, $variation_dba);
  print "Got ".(scalar @{$results->{'phenotypes'}})." new phenotypes \n" if $debug ;

  #get source id
  my $source_id = $self->get_or_add_source(\%source_info,$variation_dba);
  print STDOUT "$source_info{source} source_id is $source_id\n" if ($debug);

  # save phenotypes
  $self->set_skip_synonyms(1); #TODO:check if it was set in BasePhenotypeAnnotation 
  $self->save_phenotypes(\%source_info, $results, $core_dba, $variation_dba);

  my %param_source = (source_name => $source_info{source_name},
                      type => $source_info{object_type});
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species'),
                               workdir => $workdir
                             });
}

sub write_output {
  my $self = shift;

  $self->dataflow_output_id($self->param('output_ids'), 1);
  print "Passing MGI import for checks\n" if $self->param('debug_mode');
}

1;

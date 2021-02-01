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


=head1 ImportOrphanet

This module imports Orphanet data. The module fetches the data from
www.orphadata.org.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportOrphanet;

use warnings;
use strict;

use File::Path qw(make_path);
use LWP::Simple;
use XML::LibXML;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $run_type     = $self->required_param('run_type');

  $self->debug($self->param('debug_mode'));

  my $orphanet_data_url = 'http://www.orphadata.org/data/xml/en_product6.xml';

  %source_info = (source_description => 'The portal for rare diseases and drugs',
                  source_url => 'https://www.orpha.net/',
                  object_type => 'Gene',
                  source_status => 'germline',
                  source_name => 'Orphanet',        #source name in the variation db
                  source_name_short => 'Orphanet',  #source identifier in the pipeline
                  data_types => 'phenotype_feature',
                #  source_version is set later in import based on record contained in input file
                  );

  #get input file Orphanet:
  my $workdir = $pipeline_dir."/".$source_info{source_name_short}."/".$species;
  make_path($workdir) or die "Failed to create $workdir $!\n";
  $self->workdir($workdir);

  open(my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  my $file_orphanet = "en_product6.xml";
  getstore($orphanet_data_url, $workdir."/".$file_orphanet) unless -e $workdir."/".$file_orphanet;

  print $logFH "Found files (".$workdir."/".$file_orphanet."), will skip new fetch\n" if -e $workdir."/".$file_orphanet;
  $self->param('orphanet_file', $file_orphanet);
}

sub run {
  my $self = shift;

  my $file_orphanet = $self->required_param('orphanet_file');

  # dump and clean pre-existing phenotype features
  $self->dump_phenotypes($source_info{source_name}, 1);

  # get phenotype data + save it (all in one method)
  my ($results,$source_date) = $self->parse_input_file($file_orphanet);
  $self->print_logFH("Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n") if ($self->debug);

  #save source_date
  $source_info{source_version} = $source_date;

  # save phenotypes
  $self->save_phenotypes(\%source_info, $results);

  my %param_source = (source_name => $source_info{source_name_short},
                      type => $source_info{object_type});
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species'),
                               run_type => $self->required_param('run_type'),
                             });
}

sub write_output {
  my $self = shift;

  $self->print_pipelogFH("Passing $source_info{source_name_short} import (".$self->required_param('species').") for checks (check_phenotypes)\n") if ($self->debug);
  close($self->logFH) if defined $self->logFH;
  close($self->errFH) if defined $self->errFH;
  close($self->pipelogFH) if defined $self->pipelogFH;

  $self->dataflow_output_id($self->param('output_ids'), 2);
}


=head2 parse_input_file

  Arg [1]    : string $infile
               The input file name.
  Example    : $results = $obj->parse_input_file($infile)
  Description: Parse phenotypes from Orphanet species input file, uses gene symbols lookup in core
  Returntype : hashref with results (key 'phenotypes')
  Exceptions : none

=cut

sub parse_input_file {
  my ($self, $infile) = @_;

  my $ga = $self->core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor\n") unless defined($ga);

  my $errFH1;
  open($errFH1, ">", $self->workdir."/".'log_import_err_'.$infile) || die ("Could not open file for writing: $!\n");

  my @phenotypes;

  my $xml_parser   = XML::LibXML->new();
  my $orphanet_doc = $xml_parser->parse_file($self->workdir."/".$infile);

  #get date:
  my $first_node = $orphanet_doc->findnodes('JDBOR')->get_node(1);
  $first_node->getAttribute('date') =~ /(\d+)-(\d{2})-(\d{2})/;
  my $date = $1.$2.$3;

  foreach my $disorder ($orphanet_doc->findnodes('JDBOR/DisorderList/Disorder')) {
    my ($orpha_number_node) = $disorder->findnodes('./OrphaCode');
    my $orpha_number = $orpha_number_node->to_literal;
    my ($name_node) = $disorder->findnodes('./Name');
    my $name = $name_node->to_literal;

    # Special encoding conversion
    utf8::encode($name);

    my @gene_nodes = $disorder->findnodes('./DisorderGeneAssociationList/DisorderGeneAssociation/Gene');

    foreach my $gene_node(@gene_nodes) {
      my $ref;
      my %ens_ids;

      #get the HGNC xref
      foreach my $external_reference_node ($gene_node->findnodes('./ExternalReferenceList/ExternalReference')) {
        my ($source_node) = $external_reference_node->findnodes('./Source');
        if ($source_node->to_literal =~ /HGNC/) {
          my ($ref_node) = $external_reference_node->findnodes('./Reference');
          $ref = $ref_node->to_literal;
        }
        if ($source_node->to_literal =~ /ensembl/i) {
          my ($ens_node) = $external_reference_node->findnodes('./Reference');
          $ens_ids{$ens_node->to_literal} = 1;
        }
      }
      if (defined($ref) || scalar keys %ens_ids) {
        my $genes = [];

        if(scalar keys %ens_ids) {
          foreach my $ens_id(keys %ens_ids) {
            my $g = $ga->fetch_by_stable_id($ens_id);
            push @$genes, $g if defined($g);
          }
        }

        if(scalar @$genes == 0 && defined($ref)) {
          my $genes = $ref =~ /^\d+$/ ? $ga->fetch_all_by_description('%HGNC:'.$ref.']%') : $ga->fetch_all_by_external_name($ref, 'HGNC');

          # we don't want any LRG genes
          @$genes = grep {$_->stable_id !~ /^LRG_/} @$genes;

          # and we don't want any duplicates if possible
          @$genes = grep {$ens_ids{$_->stable_id}} @$genes if scalar keys %ens_ids;
        }

        if ( scalar(@$genes) > 0) {
          if(scalar @$genes != 1) {
            print $errFH1 "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for HGNC ID $ref or Ens IDs ".join(",", keys %ens_ids)."\n";
          }

          next unless scalar @$genes;

          foreach my $gene(@$genes) {
            push @phenotypes, {
              'id' => $gene->stable_id,
              'description' => $name,
              'external_id' => $orpha_number,
              'accessions'  => [ 'Orphanet:' . $orpha_number],
               ontology_mapping_type =>'is', 
              'seq_region_id' => $gene->slice->get_seq_region_id,
              'seq_region_start' => $gene->seq_region_start,
              'seq_region_end' => $gene->seq_region_end,
              'seq_region_strand' => $gene->seq_region_strand
            };
          }
        }
      }
    }
  }
  close($errFH1);

  my %result = ('phenotypes' => \@phenotypes);
  return (\%result, $date);
}

1;

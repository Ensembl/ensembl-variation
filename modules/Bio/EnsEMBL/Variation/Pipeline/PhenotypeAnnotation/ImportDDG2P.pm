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


=head1 ImportDDG2P

This module imports DDG2P (Developmental Disorders Genotype-to-Phenotype) data. The module fetched the data from
EBI Gene2Phenotype project https://www.ebi.ac.uk/gene2phenotype/.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportDDG2P;

use warnings;
use strict;

use File::Path qw(make_path);
use POSIX qw(strftime);
use IO::Uncompress::Gunzip qw(gunzip);
use Text::CSV;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');

  $self->debug($self->param('debug_mode'));

  %source_info = (source_description => 'Developmental Disorders Genotype-to-Phenotype Database',
                  source_url => 'https://decipher.sanger.ac.uk/',
                  object_type => 'Gene',
                  source_status => 'germline',
                  source_name => 'DDG2P',       #source name in the variation db
                  source_name_short => 'DDG2P', #source identifier in the pipeline
                  source_version => strftime("%Y%m%d", localtime), # it is current month
                  data_types        => 'phenotype_feature',
                  );

  my $workdir = $pipeline_dir."/".$source_info{source_name}."/".$species;
  make_path($workdir) or die "Failed to create $workdir $!\n";
  $self->workdir($workdir);

  open(my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  #get input file DDG2P:
  my $ddg2p_url = 'https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz';
  my $dateStrURL = strftime("%d_%m_%Y", localtime);
  my $file_ddg2p_gz = "DDG2P_$dateStrURL.csv.gz";

  print $logFH "Found files (".$workdir."/".$file_ddg2p_gz."), will skip new fetch\n" if -e $workdir."/".$file_ddg2p_gz;
  my $resHTTPcode = qx{curl -L -w %{http_code} -X GET $ddg2p_url -o $workdir/$file_ddg2p_gz} unless -e $workdir."/".$file_ddg2p_gz ;
  print $errFH "WARNING: File cound not be retrieved (HTTP code: $resHTTPcode)" if defined($resHTTPcode) && $resHTTPcode != 200;

  my $file_ddg2p = "DDG2P.csv";
  gunzip $workdir."/".$file_ddg2p_gz => $workdir."/".$file_ddg2p;
  $self->param('ddg2p_file', $file_ddg2p);
}

sub run {
  my $self = shift;

  my $file_ddg2p = $self->required_param('ddg2p_file');

  # dump and clean pre-existing phenotype features
  $self->dump_phenotypes($source_info{source_name}, 1);

  #get source id
  my $source_id = $self->get_or_add_source(\%source_info);
  $self->print_logFH("$source_info{source_name} source_id is $source_id\n") if ($self->debug);

  # get phenotype data + save it (all in one method)
  my $results = $self->parse_input_file($file_ddg2p);
  $self->print_logFH("Got ".(scalar @{$results->{'phenotypes'}})." new phenotypes \n") if ($self->debug);

  # save phenotypes
  $self->save_phenotypes(\%source_info, $results);

  my %param_source = (source_name => $source_info{source_name},
                      type => $source_info{object_type});
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species'),
                               run_type => $self->required_param('run_type'),
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


=head2 parse_input_file

  Arg [1]    : string $infile
               The input file name.
  Example    : $results = $obj->parse_input_file($infile)
  Description: Specific parsing method for DDG2P data, uses gene symbols lookup in core
  Returntype : hashref with results (key 'phenotypes')
  Exceptions : none

=cut

sub parse_input_file {
  my ($self, $infile) = @_ ;

  my $ga = $self->core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor\n") unless defined($ga);

  my $errFH1;
  open($errFH1, ">", $self->workdir."/".'log_import_err_'.$infile) ;

  my @phenotypes;
  my $fh;

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open($fh, "zcat ".$self->workdir."/$infile |") || die ("Could not open $infile for reading: $!\n");
  }
  else {
    open($fh,'<',$self->workdir."/".$infile) || die ("Could not open $infile for reading: $!\n");
  }

  my %headers;

  my $csv = Text::CSV->new ( { binary => 1 } )  # should set binary attribute.
            || die ("Cannot use CSV: ".Text::CSV->error_diag ()."\n");

  # Get columns headers
  $csv->column_names ($csv->getline ($fh));

  # Read through the file and parse out the desired fields
  while (my $content = $csv->getline_hr ($fh)) {

    # get data from the line
    my $symbol  = $content->{'gene symbol'};
    my $allelic = $content->{'allelic requirement'};
    my $mode    = $content->{'mutation consequence'};
    my $phen    = $content->{'disease name'};
    my $id      = $content->{'disease mim'};
    my @accns   = split/\;/,$content->{'phenotypes'};

    if ($symbol && $phen) {
      $phen =~ s/\_/ /g;

      my $genes = $ga->fetch_all_by_external_name($symbol, 'HGNC');

      # we don't want any LRG genes
      @$genes = grep {$_->stable_id !~ /^LRG_/} @$genes;

      # try restricting by name
      if (scalar(@$genes) > 1) {
        my @tmp = grep {$_->external_name eq $symbol} @$genes;
        $genes = \@tmp if scalar @tmp;
      }

      if (scalar @$genes != 1) {
        print $errFH1 "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for HGNC ID $symbol\n";
      }

      next unless scalar @$genes;

      foreach my $gene(@$genes) {
        push @phenotypes, {
          'id' => $gene->stable_id,
          'description' => $phen,
          'external_id' => $id,
          'seq_region_id' => $gene->slice->get_seq_region_id,
          'seq_region_start' => $gene->seq_region_start,
          'seq_region_end' => $gene->seq_region_end,
          'seq_region_strand' => $gene->seq_region_strand,
          'mutation_consequence' => $mode,
          'inheritance_type' => $allelic,
          'accessions' => \@accns,
          ontology_mapping_type =>'involves' 
        };
      }
    }
  }
  close($fh);
  close($errFH1);

  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

1;

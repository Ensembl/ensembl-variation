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


=head1 ImportZFIN

This module imports ZFIN (The Zebrafish Information Network) data.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportZFIN;

use strict;
use warnings;

use File::Path qw(make_path);
use File::stat;
use POSIX qw(strftime);
use LWP::Simple;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;

sub fetch_input {
  #create output folder structure and fetches input files
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $run_type     = $self->required_param('run_type');

  $self->debug($self->param('debug_mode'));

  # import specific constants
  %source_info = (source_description => 'The zebrafish model organism database',
                  source_url => 'https://www.zfin.org/',
                  object_type => 'Gene',
                  #source_version  will be set based on the date of the fetched input file  (year/month/day-> yyyymmdd)

                  source_status => 'germline',

                  source_name => 'ZFIN',        #source name in the variation db
                  source_name_short => 'ZFIN',  #source identifier in the pipeline
                  data_types => 'phenotype_feature',
                  );
  my $inputFile = 'phenoGeneCleanData_fish.txt';
  my $zfin_url = 'https://zfin.org/downloads/'.$inputFile;

  #create workdir folder
  my $workdir = $pipeline_dir."/".$source_info{source_name_short}."/".$species;
  make_path($workdir) or die "Failed to create $workdir $!\n";
  $self->workdir($workdir);

  open(my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  getstore($zfin_url, $workdir."/".$inputFile) unless -e $workdir."/".$inputFile;
  print $logFH "Found files (".$workdir."/".$inputFile.") and will skip new fetch\n" if -e $workdir."/".$inputFile;
  $source_info{source_version} = strftime("%Y%m%d", localtime(stat($workdir."/".$inputFile)->mtime));

  $self->param('zfin_file', $inputFile);
}

sub run {
  my $self = shift;

  #Process QTLs file
  my $input_file = $self->required_param('zfin_file');   #Go through files and parse them in the correct format

  # dump and clean pre-existing phenotype features
  $self->dump_phenotypes($source_info{source_name}, 1);

  # get phenotype data
  my $results = $self->parse_input_file($input_file);
  $self->print_logFH("Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n") if ($self->debug);

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
  Description: Parse phenotypes from ZFIN input file, uses gene symbols lookup in core
  Returntype : hashref with results (key 'phenotypes')
  Exceptions : none

=cut

sub parse_input_file {
  my ($self, $infile) = @_;

  my $ga = $self->core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor\n") unless defined($ga);

  my $errFH1;
  open($errFH1, ">", $self->workdir."/"."log_import_err_".$infile) ;

  my @phenotypes;

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open(IN, "zcat ".$self->workdir."/$infile |") || die ("Could not open $infile for reading\n");
  }
  else {
    open(IN,'<',$self->workdir."/".$infile) || die ("Could not open $infile for reading\n");
  }

  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;

    my @data = split /\t/, $_;

    my $symbol  = $data[1];
    my $gene_id = $data[2];
    my $phen    = $data[4];

    for my $i(6, 8, 10, 12) {
      $phen .= ($phen ? ', ' : '').$data[$i] if $data[$i];
    }

    if($symbol && $phen) {
      my $genes = $ga->fetch_all_by_external_name($symbol);

      # try restricting by name
      if(scalar(@$genes) > 1) {
        my @tmp = grep {$_->external_name eq $symbol} @$genes;
        $genes = \@tmp if scalar @tmp;
      }

      if(scalar @$genes != 1) {
        print $errFH1 "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for gene ID $symbol\n";
      }

      next unless scalar @$genes;

      foreach my $gene(@$genes) {
        push @phenotypes, {
          'id' => $gene->stable_id,
          'description' => $phen,
          'external_id' => $gene_id,
          'seq_region_id' => $gene->slice->get_seq_region_id,
          'seq_region_start' => $gene->seq_region_start,
          'seq_region_end' => $gene->seq_region_end,
          'seq_region_strand' => $gene->seq_region_strand,
        };
      }
    }
  }
  close (IN);
  close ($errFH1);

  return {'phenotypes' => \@phenotypes};
}

1;


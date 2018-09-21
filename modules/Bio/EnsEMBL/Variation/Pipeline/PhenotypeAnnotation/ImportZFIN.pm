=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

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


package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportZFIN;

use strict;
use warnings;

use File::Copy;
use File::Path qw(make_path remove_tree);
use LWP::Simple;
use Data::Dumper;

use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;
my $workdir;
my $debug;
my $inputFile;

my $core_dba;
my $variation_dba;
my $phenotype_dba;

#-source zfin

sub fetch_input {
    #create output folder structure and fetches input files 
    my $self = shift;

    my $pipeline_dir = $self->required_param('pipeline_dir');
    my $species      = $self->required_param('species');

    $core_dba    = $self->get_species_adaptor('core');
    $variation_dba  = $self->get_species_adaptor('variation');
    $phenotype_dba  = $variation_dba->get_PhenotypeAdaptor;

    $debug        = $self->param('debug_mode');

    # import specific constants
    %source_info = (source_description => 'The zebrafish model organism database',
                    source_url => 'http://www.zfin.org/',
                    object_type => 'Gene',

                    source_name => 'ZFIN',
                    source_version => $self->required_param('zfin_version'),

                    source_status => 'germline',
                    source => 'ZFIN',
                    );
    my $zfin_url = 'http://zfin.org/downloads/phenoGeneCleanData_fish.txt';
    $inputFile = 'phenoGeneCleanData_fish.txt';

    #create workdir folder
    $workdir = $pipeline_dir."/".$source_info{source_name}."/".$species;
    make_path($workdir);
    local (*STDOUT, *STDERR);
    open STDOUT, ">", $workdir."/".'log_import_ZFIN_out'; #TODO: what is best error/out log naming convention?
    open STDERR, ">", $workdir."/".'log_import_ZFIN_err';

    getstore($zfin_url, $workdir."/".$inputFile) unless -e $workdir."/".$inputFile;
    print "Found files (".$workdir."/".$inputFile.") and will skip new fetch\n" if -e $workdir."/".$inputFile;
}

sub run {
    my $self = shift;
    
    local (*STDOUT, *STDERR);
    open STDOUT, ">", $workdir."/".'log_import_out_'.$inputFile; #TODO: what is best error/out log naming convention?
    open STDERR, ">", $workdir."/".'log_import_err_'.$inputFile;

    # get phenotype data
    my $results = parse_zfin($workdir."/".$inputFile);
    print "Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n" if $debug ;

    # save phenotypes
    $self->save_phenotypes(\%source_info, $results, $core_dba, $variation_dba);

    my %param_source = (source_name => $source_info{source_name},
                        type => ['Gene']);
    $self->param('output_ids', { source => \%param_source,
                                 species => $self->required_param('species')
                               });
    
}

sub write_output {
  my $self = shift;

  $self->dataflow_output_id($self->param('output_ids'), 1);
  print "Passing RGD import for checks\n" if $self->param('debug_mode');
}


# ZFIN specific phenotype parsing method for txt files
sub parse_zfin {
  my $infile = shift;

  my $ga = $core_dba->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor") unless defined($ga);

  my @phenotypes;

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
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
        print STDERR "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for gene ID $symbol\n";
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

  return {'phenotypes' => \@phenotypes};
}

1;


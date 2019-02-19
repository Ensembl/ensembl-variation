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

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportDDG2P;

use warnings;
use strict;

use File::Path qw(make_path);
use POSIX 'strftime';
use IO::Uncompress::Gunzip qw(gunzip);
use Text::CSV;
use Data::Dumper; #TODO: remove if not needed
use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation qw($variation_dba);

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;
my $workdir;
my $core_dba;
my $variation_dba;

my $debug;

sub fetch_input {
    my $self = shift;

    my $pipeline_dir = $self->required_param('pipeline_dir');
    my $species      = $self->required_param('species');

    $core_dba    = $self->get_species_adaptor('core');
    $variation_dba  = $self->get_species_adaptor('variation'); #TODO: why did the init -> Base class -> this not work?

    $debug        = $self->param('debug_mode');

    my $dateStr = strftime "%Y%m%d", localtime;

    %source_info = (source_description => 'Developmental Disorders Genotype-to-Phenotype Database',
                    source_url => 'http://decipher.sanger.ac.uk/',
                    object_type => 'Gene',

                    source_name => 'DDG2P', #TODO: figure out where this is used
                    source_version => $dateStr, # it is current month
                    source_status => 'germline', #TODO: check statuses
                    source => 'ddg2p',
                    );

    $workdir = $pipeline_dir."/".$source_info{source_name}."/".$species;
    make_path($workdir);

    #get input file DDG2P:
    my $ddg2p_url = 'http://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz';
    my $dateStrURL = strftime "%d_%m_%Y", localtime;
    my $file_ddg2p_gz = "DDG2P_$dateStrURL.csv.gz";

    print "Found files (".$workdir."/".$file_ddg2p_gz."), will skip new fetch\n" if -e $workdir."/".$file_ddg2p_gz;
    my $resHTTPcode = qx{curl -L -w %{http_code} -X GET $ddg2p_url -o $workdir/$file_ddg2p_gz} unless -e $workdir."/".$file_ddg2p_gz ;
    warn "WARNING: File cound not be retrieved (HTTP code: $resHTTPcode)" if defined($resHTTPcode) && $resHTTPcode != 200;
    
    my $file_ddg2p = "DDG2P.csv";
    gunzip $workdir."/".$file_ddg2p_gz => $workdir."/".$file_ddg2p; 
    $self->param('ddg2p_file', $file_ddg2p);
}

sub run {
  my $self = shift;

  my $file_ddg2p = $self->required_param('ddg2p_file');

  local (*STDOUT, *STDERR);
  open STDOUT, ">>", $workdir."/".'log_import_out_'.$file_ddg2p; #TODO: what is best error/out log naming convention?
  open STDERR, ">>", $workdir."/".'log_import_err_'.$file_ddg2p;

  #get source id
  my $source_id = $self->get_or_add_source(\%source_info,$variation_dba);
  print STDOUT "$source_info{source} source_id is $source_id\n" if ($debug);

  # get phenotype data + save it (all in one method)
  my $results = parse_ddg2p($workdir."/".$file_ddg2p, $core_dba);
  print "Got ".(scalar @{$results->{'phenotypes'}})." new phenotypes \n" if $debug ;

  # save phenotypes
  $self->save_phenotypes(\%source_info, $results, $core_dba, $variation_dba);

  my %param_source = (source_name => $source_info{source_name},
                      type => $source_info{object_type});
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species')
                             });
}

sub write_output {
  my $self = shift;

  $self->dataflow_output_id($self->param('output_ids'), 1);
  print "Passing DDG2P import for checks\n" if $self->param('debug_mode');
}

# DDG2P specific phenotype parsing method
sub parse_ddg2p {
  my $infile = shift;
  my $core_dba = shift;

  my $ga = $core_dba->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor") unless defined($ga);

  my @phenotypes;
  my $fh;

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open $fh, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open($fh,'<',$infile) or die ("Could not open $infile for reading");
  }

  my %headers;

  my $csv = Text::CSV->new ( { binary => 1 } )  # should set binary attribute.
            or die "Cannot use CSV: ".Text::CSV->error_diag ();

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
        $DB::single = 1;
        print STDERR "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for HGNC ID $symbol\n";
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

  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

1;

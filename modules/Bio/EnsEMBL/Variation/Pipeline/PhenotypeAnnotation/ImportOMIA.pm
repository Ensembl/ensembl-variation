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


=head1 ImportOMIA

This module imports OMIA (Online Mendelian Inheritance in Animals) data.

Imports gene symbols and the associated phenotypes and is not affected
by assembly/coordinate changes.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportOMIA;

use warnings;
use strict;

use File::Path qw(make_path);
use File::stat;
use POSIX 'strftime';
use LWP::Simple;
use HTTP::Tiny;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;
my $workdir;
my ($logFH, $errFH);

my $core_dba;
my $variation_dba;

my $debug;

my %species = (
  9685  => 'cat',
  9031  => 'chicken',
  9598  => 'chimpanzee',
  9913  => 'cow',
  9615  => 'dog',
  61853 => 'gibbon',
  9925  => 'goat',
  9796  => 'horse',
  9544  => 'macaque',
  10090 => 'mouse',
  13616 => 'opossum',
  9601  => 'orangutan',
  9825  => 'pig', # Sus scrofa domesticus
  9258  => 'platypus',
  10116 => 'rat',
  9940  => 'sheep',
  99883 => 'tetraodon',
  9103  => 'turkey',
  59729 => 'zebra_finch',
  7955  => 'zebrafish'
);

my %species_synonyms = (
  'cat'  => 'felis_catus',
  'chicken'  => 'gallus_gallus',
  'chimpanzee'  => 'pan_troglodytes',
  'chimp' => 'pan_troglodytes',
  'cow'  => 'bos_taurus',
  'dog'  => 'canis_familiaris',
  'gibbon' => 'nomascus_leucogenys',
  'goat'  => 'capra_hircus',
  'horse' => 'equus_caballus',
  'macaque'  => 'macaca_mulatta',
  'mouse' => 'mus_musculus',
  'opossum' => 'monodelphis_domestica',
  'orangutan' => 'pongo_abelii',
  'pig'  => 'sus_scrofa', # Sus scrofa domesticus
  'platypus' => 'ornithorhynchus_anatinus',
  'rat' => 'rattus_norvegicus',
  'sheep'  => 'ovis_aries',
  'tetraodon' => 'tetraodon_nigroviridis',
  'turkey'  => 'Meleagris_gallopavo',
  'zebra_finch' => 'taeniopygia_guttata',
  'zebrafish' => 'danio_rerio'
);

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');

  $core_dba       = $self->get_species_adaptor('core');
  $variation_dba  = $self->get_species_adaptor('variation');

  $debug        = $self->param('debug_mode');
  $self->SUPER::set_debug($self->param('debug_mode'));

  my $omia_url = 'https://omia.org/curate/causal_mutations/?format=gene_table';

  %source_info = (source_description => 'Online Mendelian Inheritance in Animals (OMIA) is a database of genes, inherited disorders and traits in more than 200 animal species (other than human and mouse)',
                  source_url => 'http://omia.angis.org.au/home/',
                  object_type => 'Gene',
                  #source_version  will be set based on the date of the fetched input file (year/month/day-> yyyymmdd)
                  source_name => 'OMIA',        #source name in the variation db
                  source_name_short => 'OMIA',  #source identifier in the pipeline
                  );

  my $workdir_fetch = $pipeline_dir."/".$source_info{source_name_short};
  make_path($workdir_fetch);
  my $file_omia = 'omia_gene_table.txt';

  $workdir = $pipeline_dir."/".$source_info{source_name_short}."/".$species;
  make_path($workdir);
  open ($logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species);
  open ($errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species);
  $self->SUPER::set_logFH($logFH);
  $self->SUPER::set_errFH($errFH);

  #get input files OMIA gene_table, this file contains multiple species
  print $logFH "Found file (".$workdir_fetch."/".$file_omia.") and will skip new fetch\n" if -e $workdir_fetch."/".$file_omia;
  getstore($omia_url, $workdir_fetch."/".$file_omia) unless -e $workdir_fetch."/".$file_omia;
  $source_info{source_version} = strftime "%Y%m%d", localtime(stat($workdir_fetch."/".$file_omia)->mtime);

  #get section specific for this species
  print $logFH "Found folder (".$workdir_fetch."/omia_split) and will skip new split\n" if -e $workdir_fetch."/omia_split";
  split_omia($workdir_fetch,$file_omia) unless -e $workdir_fetch."/"."omia_split";

  my $org_path = $pipeline_dir."/".$source_info{source_name_short}."/omia_split/omia_".$species.".txt";
  my $new_link = $workdir."/omia_".$species.".txt";
  symlink $org_path, $new_link unless -e $new_link;

  $self->param('omia_file', "omia_".$species.".txt");
}

sub run {
  my $self = shift;

  my $omia_file = $self->required_param('omia_file');

  # get phenotype data
  my $results = parse_omia($omia_file);
  print $logFH "Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n" if $debug ;

  # save phenotypes
  $self->save_phenotypes(\%source_info, $results, $core_dba, $variation_dba);

  my %param_source = (source_name => $source_info{source_name_short},
                      type => [$source_info{object_type}]);
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species')
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


# OMIA specific phenotype parsing methods
sub split_omia {
  my $workdir = shift;
  my $all_file = shift;

  my $http = HTTP::Tiny->new();
  my $server = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=';
  my $ext = "&retmode=xml";

  my $prefix = 'omia_';
  my $suffix = '.txt';

  my %data;

  open F, "< $workdir/$all_file" or die $!;
  while(<F>) {
    chomp ($_);
    next if ($_ =~ /^gene_symbol/);
    my $line = $_;
    my @line_content = split("\t",$line);
    my $taxo = $line_content[3];
    die("ERROR: Input file ($all_file) not in expected format") unless defined($taxo);

    if ($data{$taxo}) {
      push(@{$data{$taxo}}, $line);
    }
    else {
      $data{$taxo} = [$line];
    }
  }
  close(F);

  foreach my $taxo_id (sort(keys(%data))) {
    my $id = $taxo_id;

    if ($species{$taxo_id}) {
       $id = $species_synonyms{$species{$taxo_id}};
     }
     else {
      my $response = $http->get($server.$taxo_id.$ext);

      if ($response->{success} && length($response->{content})) {
        my $content = $response->{content};
        if ($content =~ /<GenbankCommonName>(.+)<\/GenbankCommonName>/) {
          $id = $1;
        }
      }
    }

    $id = lc($id);
    $id =~ s/ /_/g;
    $id =~ s/'//g;
    $id =~ s/^domestic_//g;

    make_path($workdir."/"."omia_split");
    open OUT, "> $workdir/omia_split/$prefix$id$suffix" or die $!;
    foreach my $line (@{$data{$taxo_id}}) {
      print OUT "$line\n";
    }
    close(OUT);
  }

}

sub parse_omia {
  my $infile = shift;

  my $errFH1;
  open ($errFH1, ">", $workdir."/".'log_import_err_'.$infile) ;

  my $ga = $core_dba->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor") unless defined($ga);

  my @phenotypes;

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $workdir."/".$infile |" or die ("Could not open $workdir."/".$infile for reading");
  }
  else {
    open(IN,'<',$workdir."/".$infile) or die ("Could not open $workdir."/".$infile for reading");
  }

  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;

    next if /^gene_symbol/;

    my @data = split /\t/, $_;

    my $genes = $ga->fetch_all_by_external_name($data[0]);

    if(scalar @$genes != 1) {
      print $errFH1 "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for gene name $data[0]\n";
    }

    next unless scalar @$genes;

    foreach my $gene(@$genes) {
      push @phenotypes, {
        'id' => $gene->stable_id,
        'description' => $data[5],
        'external_id' => 'OMIA'.$data[2],
        'seq_region_id' => $gene->slice->get_seq_region_id,
        'seq_region_start' => $gene->seq_region_start,
        'seq_region_end' => $gene->seq_region_end,
        'seq_region_strand' => $gene->seq_region_strand
      };
    }
  }
  close(IN);
  close ($errFH1);

  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

1;

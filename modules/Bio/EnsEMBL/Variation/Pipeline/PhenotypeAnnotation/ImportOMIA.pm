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

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportOMIA;

use warnings;
use strict;

use File::Path qw(make_path);
use LWP::Simple;
use HTTP::Tiny;

use Data::Dumper; #TODO: remove if not needed
use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(GWAS NONE);
#use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation qw($variation_dba);

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

#TODO: should contain everything about OMIA, including getting the file(s) and for all OMIA species

my %source_info;
my $workdir;
my $core_dba;
my $variation_dba;
my $phenotype_dba;

my $pubmed_prefix = 'PMID:';

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

    $core_dba    = $self->get_species_adaptor('core');
    $variation_dba  = $self->get_species_adaptor('variation'); #TODO: why did the init -> Base class -> this not work?
    $phenotype_dba  = $variation_dba->get_PhenotypeAdaptor; 

    $debug        = $self->param('debug_mode');

    #TODO: 
    # original call: perl ensembl-variation/scripts/import/import_phenotype_data.pl -host ${host} -user ${user} -pass ${pass}  -dbname ${dbname} \
    #-source omia -infile ${datadir}/gwascatalog.txt -verbose -version 20121031

    my $omia_url = 'http://omia.angis.org.au/curate/causal_mutations/?format=gene_table';

    %source_info = (source_description => 'Online Mendelian Inheritance in Animals (OMIA) is a database of genes, inherited disorders and traits in more than 200 animal species (other than human and mouse)',
                    source_url => 'http://omia.angis.org.au/home/',
                    object_type => 'Gene',

                    source_name => 'OMIA',
                    source_version => $self->required_param('omia_version'), #TODO: figure out how to get the original file_name and get : 'Date on the file name'

                    source => 'omia', #used in BasePhenotypeAnnotation multiple lines 83, ..., 972 TODO: check if really needed
                    );

    $workdir = $pipeline_dir."/".$source_info{source_name};
    my $file_omia = 'omia_gene_table.txt';

    #get input files OMIA gene_table, this file contains multiple species
    make_path($workdir);
    print "Found file (".$workdir."/".$file_omia.") and will skip new fetch\n" if -e $workdir."/".$file_omia;
    getstore($omia_url, $workdir."/".$file_omia) unless -e $workdir."/".$file_omia;

    #get section specific for this species
    print "Found folder (".$workdir."/omia_split) and will skip new split\n" if -e $workdir."/omia_split";
    split_omia($workdir,$file_omia) unless -e $workdir."/"."omia_split";

    $workdir = $pipeline_dir."/".$source_info{source_name}."/".$species;
    make_path($workdir);
    my $org_path = $pipeline_dir."/".$source_info{source_name}."/omia_split/omia_".$species.".txt";
    my $new_link = $workdir."/omia_".$species.".txt";
    symlink $org_path, $new_link unless -e $new_link;

    $self->param('omia_file', "omia_".$species.".txt");
}

sub run {
  my $self = shift;

  my $omia_file = $self->required_param('omia_file');

  local (*STDOUT, *STDERR);
  open STDOUT, ">", $workdir."/".'log_import_out_'.$omia_file; #TODO: what is best error/out log naming convention?
  open STDERR, ">", $workdir."/".'log_import_err_'.$omia_file;

  # get phenotype data
  my $results = parse_omia($workdir."/".$omia_file);
  print "Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n" if $debug ;

  # save phenotypes
  $self->save_phenotypes(\%source_info, $results, $core_dba, $variation_dba);

  my %param_source = (source_name => $source_info{source_name},
                      type => [$source_info{object_type}]);
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species')
                             });
}

sub write_output {
  my $self = shift;

  $self->dataflow_output_id($self->param('output_ids'), 1);
  print "Passing OMIA import (".$self->param('species').") for checks\n" if $self->param('debug_mode');
} #TODO: add to all multi-species runs, the log + species info as above


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
    
    next if /^gene_symbol/;
    
    my @data = split /\t/, $_;
    
    my $genes = $ga->fetch_all_by_external_name($data[0]);
    
    if(scalar @$genes != 1) {
      print STDERR "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for gene name $data[0]\n";
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
  
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

1;
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


=head1 ImportAnimalQTL

This module import AnimalQTL data. The module relies on the existance of the files to be imported.

# Currently human download is needed as a licence aggreement is needed
# steps:
# 1. follow link https://www.animalgenome.org/cgi-bin/QTLdb/index
# 2. select species e.g. Chicken QTL
# 3. Downloads -> 'All data by bp' gff format
# 4. get the *.txt.gz file and rename to species(_assembly)?.gff3.gz eg. gallus_gallus_GRCg6a.gff3.gz or gallus_gallus.gff3

# NOTE: make sure the Animal QTL db has the data mapped to the Ensembl assembly name
=cut


package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportAnimalQTL;

use strict;
use warnings;

use File::Copy;
use File::Path qw(make_path remove_tree);
use File::stat;
use POSIX 'strftime';
use PerlIO::gzip;

use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;
my $workdir;
my $debug;
my $inputFile;
my ($logFH, $errFH);

my $core_dba;
my $variation_dba;
my $ontology_db_adaptor;

#my $animalqtl_url = 'https://www.animalgenome.org/cgi-bin/QTLdb/index';
#my %animalQTL_species_url = (
#  gallus_gallus => 'https://www.animalgenome.org/cgi-bin/QTLdb/GG/download?file=gbpGG_5.0', #Gallus gallus
#  pig => 'https://www.animalgenome.org/cgi-bin/QTLdb/SS/download?file=gbpSS_11.1', #Sus scrofa
#  sheep => 'https://www.animalgenome.org/cgi-bin/QTLdb/OA/download?file=gbpOAR_3.1',  # Ovis aries
#  cow => 'https://www.animalgenome.org/cgi-bin/QTLdb/BT/download?file=gbpUMD_3.1', #Bos taurus
#  horse => 'https://www.animalgenome.org/cgi-bin/QTLdb/EC/download?file=gbpEC_2.0', #Equus caballus
#);

sub fetch_input {
  #create output folder structure and fetches input files
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $animalqtl_inputDir = $self->required_param('animalqtl_input_dir');

  $debug = $self->param('debug_mode');
  $self->SUPER::set_debug($self->param('debug_mode'));

  $core_dba    = $self->get_species_adaptor('core');
  $variation_dba  = $self->get_species_adaptor('variation');
  $ontology_db_adaptor = $self->get_adaptor('multi', 'ontology');

  # import specific constants
  %source_info = (source_description => 'The Animal Quantitative Trait Loci (QTL) database (Animal QTLdb) is designed to house all publicly available QTL and association data on livestock animal species',
                  source_url => 'http://www.animalgenome.org/cgi-bin/QTLdb/index',
                  object_type => 'QTL',
                  #source_version  will be set based on file date
                  source_status     => 'germline',
                  source_name       => 'Animal_QTLdb', #source name in the variation db
                  source_name_short => 'AnimalQTLdb',  #source identifier in the pipeline
                  threshold => $self->required_param('threshold_qtl'),
                  );

  #create workdir folder
  $workdir = $pipeline_dir."/".$source_info{source_name_short}."/".$species;
  make_path($workdir);
  open ($logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species);
  open ($errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species);
  $self->SUPER::set_logFH($logFH);
  $self->SUPER::set_errFH($errFH);

  # check if the species assembly specific file exists
  print $logFH "AnimalQTL import expects input folder with gff3 files: example format gallus_gallus.*.gff3  \n" if $debug;
  print $logFH "using input folder: $animalqtl_inputDir for species: $species \n" if $debug;

  opendir(INDIR, $animalqtl_inputDir);
  my @files = readdir(INDIR);
  closedir(INDIR);
  my $ok = 0;
  foreach my $file (@files){
    if ($file =~/^$species.*gff3$/ || $file =~/^$species.*gff3.gz$/){
      $inputFile = $file;
      $ok = 1;
    }
  }
  print $errFH "ERROR: Animal_QTLdb file not found for $species in inputDir ($animalqtl_inputDir)!\n" unless $ok;
  die "Animal_QTLdb file not found for $species in inputDir ($animalqtl_inputDir)!\n" unless $ok;

  $source_info{source_version} = strftime "%Y%m%d", localtime(stat($animalqtl_inputDir."/".$inputFile)->mtime);
  print $logFH "Found inputDir file: $inputFile \n";
  copy($animalqtl_inputDir."/".$inputFile, $workdir."/".$inputFile) unless -e $workdir."/".$inputFile;
  print $logFH "Found file (".$workdir."/".$inputFile."), will skip new copy of inputData\n" if -e $workdir."/".$inputFile;

}

sub run {
  my $self = shift;

  # get seq_region_ids
  my $seq_region_ids = $self->get_seq_region_ids($variation_dba);

  # get phenotype data
  my $results = $self->parse_animal_qtl($seq_region_ids, $workdir."/".$inputFile);
  print $logFH "Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n" if $debug ;

  # save phenotypes
  $self->save_phenotypes(\%source_info, $results, $core_dba, $variation_dba);

  my %param_source = (source_name => $source_info{source_name_short},
                      type => $source_info{object_type});
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
    print $logPipeFH "Passing AnimalQTL import (".$self->required_param('species').") for checks (check_phenotypes)\n";
    close ($logPipeFH);
  }
  $self->dataflow_output_id($self->param('output_ids'), 1);
}

# AnimalQTL specific phenotype parsing method for gff3 files
sub parse_animal_qtl {
  my ($self, $seq_region_ids, $infile) = @_;

  my $ontology_term_adaptor = $ontology_db_adaptor->get_OntologyTermAdaptor;

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

    next if /^(\#|\s)/ || !$_;

    my @data = split /\t/, $_;

    # fix chr
    $data[0] =~ s/^chr(om)?\.?//i;

    if(!defined($seq_region_ids->{$data[0]})) {
      print $errFH "WARNING: Could not find seq_region_id for chromosome name $data[0]\n";
      next;
    }

    # parse "extra" GFF fields
    my $extra = {};

    foreach my $chunk(split /\;/, $data[8]) {
      my ($key, $value) = split /\=/, $chunk;
      next unless defined($key) && defined($value);
      $value =~ s/\"|\'//g;
      $extra->{$key} = $value;
    }

    if ($extra->{'Map_Type'} eq 'Linkage' || $data[3] eq '' || $data[3] == 0) {
      print $errFH "WARNING: Could not find a precise location for the QTL ".$extra->{QTL_ID}."\n";
      next;
    }

    if ($data[4] !~ /^\d+$/) {
      print $errFH "WARNING: Could not find a numeric seq_region_end for the QTL ".$extra->{QTL_ID}."\n";
      next;
    }

    # create phenotype hash
    my $phenotype = {
      'id' => $extra->{QTL_ID},
      'description' => $extra->{Name},
      'seq_region_id' => $seq_region_ids->{$data[0]},
      'seq_region_start' => $data[3] || 1,
      'seq_region_end' => $data[4],
      'seq_region_strand' => 1
    };

    my @accessions = ();
    # CM ontology
    my $CMO_name = $extra->{CMO_name};
    if ($CMO_name) {
      my $terms = $ontology_term_adaptor->fetch_all_by_name($CMO_name, 'CMO');
      foreach my $term (@$terms) {
        push @accessions, $term->accession;
      }
    }

    # VT ontology
    my $VTO_name = $extra->{VTO_name};
    if ($VTO_name) {
      my $terms = $ontology_term_adaptor->fetch_all_by_name($VTO_name, 'VT');
      foreach my $term (@$terms) {
        push @accessions, $term->accession;
      }
    }
    # Ontology
    if (scalar(@accessions) > 0) {
      $phenotype->{accessions} = \@accessions;
      $phenotype->{ontology_mapping_type} = 'is';
    }

    if ($phenotype->{'seq_region_start'} > $phenotype->{'seq_region_end'}) {
      my $tmp_end = $phenotype->{'seq_region_end'};
      $phenotype->{'seq_region_start'} = $phenotype->{'seq_region_end'};
      $phenotype->{'seq_region_end'}   = $tmp_end;
    }

    # add additional fields if found
    $phenotype->{'study'} = $self->get_pubmed_prefix().$extra->{'PUBMED_ID'} if defined($extra->{'PUBMED_ID'} && $extra->{'PUBMED_ID'} =~ /^\d+$/);
    $phenotype->{'p_value'} = $extra->{'P-value'} if defined($extra->{'P-value'});
    $phenotype->{'f_stat'} = $extra->{'F-stat'} if defined($extra->{'F-stat'});
    $phenotype->{'lod_score'} = $extra->{'LOD-score'} if defined($extra->{'LOD-score'});
    $phenotype->{'variance'} = $extra->{'Variance'} if defined($extra->{'Variance'});
    $phenotype->{'associated_gene'} = $extra->{'Candidate_Gene_Symble'} if defined($extra->{'Candidate_Gene_Symble'});

    push @phenotypes, $phenotype;
  }

  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

1;


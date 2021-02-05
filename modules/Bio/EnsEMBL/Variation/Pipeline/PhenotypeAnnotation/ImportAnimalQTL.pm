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

use File::Path qw(make_path);
use File::stat;
use POSIX qw(strftime);
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;

# AnimalQTLdb URL set up for Ensembl: https://www.animalgenome.org/QTLdb/export/ENS83H19HZS/
#my $animalqtl_url = 'https://www.animalgenome.org/cgi-bin/QTLdb/index';
my $animalqtl_baseURL='https://www.animalgenome.org/QTLdb/export/ENS83H19HZS/';

my %animalQTL_species_url = (
  gallus_gallus => $animalqtl_baseURL.'QTL_GG_5.0.gff.txt.gz', #Gallus gallus
  sus_scrofa => $animalqtl_baseURL.'QTL_SS_11.1.gff.txt.gz', #Sus scrofa
  ovis_aries => 'https://www.animalgenome.org/QTLdb/tmp/QTL_OAR_3.1.gff.txt.gz',  # Ovis aries #TODO: replace with the one in export once it is there
  bos_taurus => $animalqtl_baseURL.'QTL_ARS_UCD1.gff.txt.gz', #Bos taurus
  equus_caballus => $animalqtl_baseURL.'QTL_EquCab2.0.gff.txt.gz', #Equus caballus
  ovis_aries_rambouillet => "",
);

my %animalQTL_species_fileNames = (
  gallus_gallus => 'QTL_gallus_gallus_gbp_6.0.gff3.gz', #Gallus gallus, remapped file
  sus_scrofa => 'QTL_sus_scrofa_gbp_11.1.gff3.gz', #Sus scrofa
  ovis_aries => 'QTL_ovis_aries_gbp_3.1.gff3.gz',  # Ovis aries, remapped file
  bos_taurus => 'QTL_bos_taurus_gbp_1.2.gff3.gz', #Bos taurus
  equus_caballus => 'QTL_equus_caballus_gbp_3.0.gff3.gz', #Equus caballus, remapped file
  ovis_aries_rambouillet => 'QTL_ovis_aries_rambouillet_gbp_1.0.gff3.gz', # remapped from ovis aries data
);

# use '0' if the data is not the same assembly and import should be skipped
my %animalQTL_species_ok = (
  gallus_gallus => 1, #Gallus gallus, remapped data, not same Ensembl assembly as AnimalQTL
  sus_scrofa => 1, #Sus scrofa
  ovis_aries => 1,  #Ovis aries: remapped data, not same Ensembl assembly as AnimalQTL
  bos_taurus => 1, #Bos taurus
  equus_caballus => 1, #Equus caballus: remapped data, not same Ensembl assembly as AnimalQTL
  ovis_aries_rambouillet => 1, #remapped data from Ovis aries (ovis_aries_variation)
);


sub fetch_input {
  #create output folder structure and fetches input files
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $threshold    = $self->param('threshold_qtl');

  # import specific constants
  %source_info = (source_description => 'The Animal Quantitative Trait Loci (QTL) database (Animal QTLdb) is designed to house all publicly available QTL and association data on livestock animal species',
                  source_url => 'https://www.animalgenome.org/cgi-bin/QTLdb/index',
                  object_type => 'QTL',
                  #source_version  will be set based on file date
                  source_status     => 'germline',
                  source_name       => 'Animal_QTLdb', #source name in the variation db
                  source_name_short => 'AnimalQTLdb',  #source identifier in the pipeline
                  data_types        => 'phenotype_feature,study',
                  threshold => $threshold,
                  );

  #create workdir folder
  my $workdir = $pipeline_dir."/".$source_info{source_name_short}."/".$species;
  unless (-d $workdir) {
    my $err;
    make_path($workdir, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }
  $self->workdir($workdir);

  # if assembly does not match, don't import data
  return unless $animalQTL_species_ok{$species};

  $self->debug($self->param('debug_mode'));

  open(my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  # check if the species assembly specific file exists
  my $animalqtl_inputDir = $pipeline_dir."/".$source_info{source_name_short}."/animalqtl_data";
  print $logFH "AnimalQTL import expects input folder with gff3 files: example format QTL_gallus_gallus.*.gff3  \n" if ($self->debug);
  print $logFH "using input folder: $animalqtl_inputDir for species: $species \n" if ($self->debug);

  my $inputFile = $animalqtl_inputDir."/".$animalQTL_species_fileNames{$species};
  my $url=$animalQTL_species_url{$species};

  # if the folder does not exist, try to fetch from ulr
  if (! -d $animalqtl_inputDir) {
    make_path($animalqtl_inputDir) or die "Failed to create $animalqtl_inputDir $!\n";
    my $fetch_cmd = "wget --content-disposition --no-check-certificate -O $inputFile \"$url\"";
    my $return_value = $self->run_cmd($fetch_cmd)
      unless -e $animalqtl_inputDir."/".$animalQTL_species_fileNames{$species};
    die ("File fetch failed code: $return_value!\n") unless defined($return_value) && $return_value == 0;
  } else {
    opendir(INDIR, $animalqtl_inputDir);
    my @files = readdir(INDIR);
    closedir(INDIR);
    my $ok = 0;
    foreach my $file (@files){
      if ($file =~/^QTL_$species.*gff3$/ || $file =~/^QTL_$species.*gff3.gz$/){
        $inputFile = $file;
        $ok = 1;
        last;
      }
    }
    # if directory exists and file not found try to fetch the file
    if (!$ok) {
      my $fetch_cmd = "wget --content-disposition --no-check-certificate -O $inputFile \"$url\"";
      my $return_value = $self->run_cmd($fetch_cmd)
        unless -e $animalqtl_inputDir."/".$animalQTL_species_fileNames{$species};
      die ("File fetch failed code: $return_value!\n") unless defined($return_value) && $return_value == 0;
      $ok = 1;
    }
    print $errFH "ERROR: Animal_QTLdb file not found for $species in inputDir ($animalqtl_inputDir)!\n" unless $ok;
    die ("Animal_QTLdb file not found for $species in inputDir ($animalqtl_inputDir)!\n") unless $ok;
  }

  #allow time between file download and file read for system to sync
  sleep(45);

  #fetch coreDB assembly, in future this should be tested against
  my $gc =  $self->core_db_adaptor->get_adaptor('GenomeContainer');
  $self->param('species_assembly', $gc->get_version);  #'GRCg6a' for gallus_gallus
  print $logFH 'INFO: Found core species_assembly:'. $self->param('species_assembly'). "\n" if ($self->debug);

  $source_info{source_version} = strftime("%Y%m%d", localtime(stat($animalqtl_inputDir."/".$inputFile)->mtime));
  print $logFH "Found inputDir file: $inputFile \n";
  if ( -e $workdir."/".$inputFile) {
    print $logFH "Found file (".$workdir."/".$inputFile."), will skip new copy of inputData\n";
  } else {
    my $cp_cmd = "cp -p $animalqtl_inputDir/$inputFile $workdir/$inputFile";
    my $return_value = $self->run_cmd($cp_cmd);
    die ("File copy failed code: $return_value!\n") unless defined($return_value) && $return_value == 0;
  }

  $self->param('qtl_file', $inputFile);
}

sub run {
  my $self = shift;

  return unless $animalQTL_species_ok{$self->required_param('species')};

  my $file_qtl = $self->required_param('qtl_file');

  # dump and clean pre-existing phenotypes
  $self->dump_phenotypes($source_info{source_name}, 1);

  # get seq_region_ids
  my $seq_region_ids = $self->get_seq_region_ids();

  # get phenotype data
  my $results = $self->parse_input_file($seq_region_ids, $file_qtl);
  $self->print_logFH("Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n") if ($self->debug);

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

  # default branch is also check source branch
  if ($animalQTL_species_ok{$self->required_param('species')}){
    $self->print_pipelogFH("Passing $source_info{source_name_short} import (".$self->required_param('species').") for checks (check_phenotypes)\n") if ($self->debug);
    close($self->logFH) if (defined $self && defined $self->logFH) ;
    close($self->errFH) if (defined $self && defined $self->errFH) ;
    close($self->pipelogFH) if (defined $self && defined $self->pipelogFH) ;

  } else {
    open(my $pipelogFH, ">", $self->workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$self->required_param('species')) || die ("Failed to open file: $!\n");
    print $pipelogFH "Ensembl species has different assembly than AnimalQTL, will exit!\n";
    close($pipelogFH);
  }

  $self->dataflow_output_id($self->param('output_ids'), 2);

}


=head2 parse_input_file

  Arg [1]    : arrayref $seq_region_ids
               The array of string seq_region_id(s)
  Arg [2]    : string $infile
               The input file name.
  Example    : $result = $self->parse_input_file($seq_region_ids, $infile)
  Description: Specific parsing method for AnimalQTL gff3 phenotype files.
  Returntype : hashref with results (key 'phenotypes')
  Exceptions : none

=cut

sub parse_input_file {
  my ($self, $seq_region_ids, $infile) = @_;

  my $errFH1;
  open($errFH1, ">", $self->workdir."/".'log_import_err_'.$infile) || die ("Could not open file (".$self->workdir."/log_import_err_".$infile.") for writing $!\n") ;

  my $ontology_term_adaptor = $self->ontology_db_adaptor->get_OntologyTermAdaptor;

  my @phenotypes;

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open(IN, "zcat ".$self->workdir."/$infile | ") || die ("Could not open $infile for reading\n");
  }
  else {
    open(IN,'<',$self->workdir."/".$infile) || die ("Could not open $infile for reading\n");
  }

  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;

    next if /^(\#|\s)/ || !$_;

    my @data = split /\t/, $_;

    # fix chr
    $data[0] =~ s/^chr(om)?\.?//i;

    if(!defined($seq_region_ids->{$data[0]})) {
      print $errFH1 "WARNING: Could not find seq_region_id for chromosome name $data[0]\n";
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
      print $errFH1 "WARNING: Could not find a precise location for the QTL ".$extra->{QTL_ID}."\n";
      next;
    }

    if ($data[4] !~ /^\d+$/) {
      print $errFH1 "WARNING: Could not find a numeric seq_region_end for the QTL ".$extra->{QTL_ID}."\n";
      next;
    }

    my $phenotype_description= $extra->{Name};
    if ($phenotype_description eq '' && defined ($extra->{trait})){
      $phenotype_description = $extra->{trait};
    }

    if ($phenotype_description eq '') {
      print $errFH1 "WARNING: Could not find a phenotype description for the QTL ".$extra->{QTL_ID}."\n";
      next;
    }

    # create phenotype hash
    my $phenotype = {
      'id' => $extra->{QTL_ID},
      'description' => $phenotype_description,
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
  close(IN);
  close ($errFH1);

  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

1;


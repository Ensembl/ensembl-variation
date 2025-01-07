=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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


=head1 ImportCGC

This module imports Cancer Gene Census data. The module fetched the data from
OpenTargets project https://www.targetvalidation.org/downloads/data

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportCGC;

use warnings;
use strict;

use File::Path qw(make_path);
use File::Basename;
use DateTime;
use JSON;

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;
my $debug;

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $repo_dir     = $self->required_param('repo_dir');

  $debug = $self->param('debug_mode');

  %source_info = (source_description => 'Catalog of genes of which mutations have been causally implicated in cancer',
                  source_url => 'https://cancer.sanger.ac.uk/census',
                  object_type => 'Gene',
                  #source_version  will be set based on the date in the file name (year/month-> yyyymm)
                  source_status => 'somatic',
                  source_name => 'Cancer Gene Census',     #source name in the variation db
                  source_name_short => 'CGC', #source identifier in the pipeline
                  data_types => 'phenotype_feature,study',
                  );

  my $workdir = $pipeline_dir."/".$source_info{source_name_short}."/".$species;
  unless (-d $workdir) {
    my $err;
    make_path($workdir, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }
  $self->workdir($workdir);

  open(my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  # get current date
  my $dt = DateTime->now;
  my $year = sprintf("%02d", $dt->year % 100);
  my $month = sprintf ("%02d",$dt->month);

  # get input file CGC via OpenTargets, get latest published file:
  my $script_dir = $repo_dir . "/ensembl-variation/scripts/python";
  my $file_opent = "cgc_input_latest.json";
  system("python3 $script_dir/download_cgc_file.py -d $workdir -r latest") unless -e $workdir."/".$file_opent;

  my $found = (-e $workdir."/".$file_opent ? 1: 0);

  if(!$found) {
    print $errFH "Input file ($file_opent) is missing\n";
    die "Input file is missing: check download script $script_dir/download_cgc_file.py\n";
  }

  print $logFH "INFO: Found input file ($file_opent)\n" ;
  my $valid_file = $self->validate_input_file($file_opent);
  if(!$valid_file) {
    print $errFH "Input file ($file_opent) has incorrect format\n";
    die ("Input file ($file_opent) has incorrect format or it's empty\n");
  }

  # Source version to be used to update the table source
  $source_info{source_version} = "20$year$month";

  my $file_cgc_ensembl = "cgc_ensembl_efo.txt";
  $self->print_logFH("Retrieving data from Ensembl and EFO based on CancerGeneCensus file\n");
  $self->print_logFH("INFO: Found file ($file_cgc_ensembl), will skip new fetch\n") if -e $workdir."/".$file_cgc_ensembl;
  $self->get_input_file($file_opent, $file_cgc_ensembl) unless -e  $workdir."/".$file_cgc_ensembl;
  $self->print_logFH("Done retrieving Ensembl core and ontology data, if file was already present it will be reused.\n");

  $self->param('cgc_file', $file_cgc_ensembl);
}

sub run {
  my $self = shift;

  my $file_cgc = $self->required_param('cgc_file');

  # dump and clean pre-existing phenotype features
  $self->dump_phenotypes($source_info{source_name}, 1);

  #get source id
  my $source_id = $self->get_or_add_source(\%source_info);
  $self->print_logFH("$source_info{source_name_short} source_id is $source_id\n") if ($debug);

  # get phenotype data + save it (all in one method)
  my $results = $self->parse_input_file($file_cgc);
  $self->print_logFH( "Parsed ".(scalar @{$results->{'phenotypes'}})." new phenotypes \n") if ($debug);

  # save phenotypes
  $self->save_phenotypes(\%source_info, $results) unless scalar(@{$results->{'phenotypes'}}) == 0;

  my %param_source = (source_name => $source_info{source_name_short},
                      type => $source_info{object_type});
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species'),
                               run_type => $self->required_param('run_type'),
                             });
  $self->clean_dir;
}

sub write_output {
  my $self = shift;

  $self->print_pipelogFH("Passing $source_info{source_name_short} import (".$self->required_param('species').") for checks (check_phenotypes)\n") if ($debug);
  close($self->logFH) if defined $self->logFH ;
  close($self->errFH) if defined $self->errFH ;
  close($self->pipelogFH) if defined $debug ;

  $self->dataflow_output_id($self->param('output_ids'), 2);

}

=head2 validate_input_file

  Arg [1]    : string $infile
               The input file name
  Example    : $obj->validate_input_file($infile)
  Description: Validates the format of the input file.
  Returntype : Boolean
  Exceptions : none

=cut

sub validate_input_file {
  my ($self, $input_file) = @_;

  my $valid = 1;

  # empty file
  if(-z $self->workdir."/".$input_file) {
    $valid = 0;
  }

  my $read;
  if($input_file =~ /gz$/) {
    open($read, "zcat ".$self->workdir."/$input_file |") || die ("Validate input file: could not open $input_file for reading: $!\n");
  }
  else {
    open($read,'<',$self->workdir."/".$input_file) || die ("Validate input file: could not open $input_file for reading: $!\n");
  }

  my $line = <$read>;
  chomp $line;
  my $json_hash = decode_json($line);

  if(!defined $json_hash->{'datasourceId'}) {
    $valid = 0;
  }
  elsif(!defined $json_hash->{'datatypeId'}) {
    $valid = 0;
  }
  elsif(!defined $json_hash->{'targetId'}) {
    $valid = 0;
  }

  close($read);

  return $valid;
}

=head2 get_input_file

  Arg [1]    : string $infile
               The input file name
  Arg [1]    : string $outfile
               The output file name.
  Example    : $obj->get_input_file($infile,$outfile)
  Description: Specific parsing method for OpenTarget json files:
               writes output file with only the relevant information from Cancer Gene Census.
  Returntype : none
  Exceptions : none

=cut

sub get_input_file {
  my ($self, $input_file, $output_file) = @_;

  my $ga = $self->core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor\n") unless defined($ga);
  my $ota = $self->ontology_db_adaptor->get_OntologyTermAdaptor;
  die("ERROR: Could not get ontology term adaptor\n") unless defined($ota);

  my $errFH1;
  open($errFH1, ">", $self->workdir."/".'log_import_err_'.$input_file) || die ("Failed to open file: $!\n");

  my %efos;
  my %data;

  if($input_file =~ /gz$/) {
    open(IN, "zcat ".$self->workdir."/$input_file |") || die ("Could not open $input_file for reading: $!\n");
  }
  else {
    open(IN,'<',$self->workdir."/".$input_file) || die ("Could not open $input_file for reading: $!\n");
  }

  while(<IN>) {
    chomp $_;
    my $json_hash = decode_json($_);

    my $source = $json_hash->{'datasourceId'};
    next unless ($source =~ /cancer_gene_census/);

    my $type = $json_hash->{'datatypeId'};
    my $gene_id = $json_hash->{'targetId'}; # gene ID (ENSG00000147889)
    my $pmids = parse_publications($json_hash->{'literature'}) if (defined $json_hash->{'literature'});
    my $phenotype_id = $json_hash->{'diseaseId'} if (defined $json_hash->{'diseaseId'});
    my $so_term;
    my $n_mutated_samples;
    my $n_samples_tested;
    my $n_samples_mutation;
    if (defined $json_hash->{'mutatedSamples'}) {
      my $mutated_samples_info = $json_hash->{'mutatedSamples'}->[0];
      $so_term = $mutated_samples_info->{'functionalConsequenceId'} if(defined $mutated_samples_info->{'functionalConsequenceId'});
      $n_mutated_samples = $mutated_samples_info->{'numberMutatedSamples'} if(defined $mutated_samples_info->{'numberMutatedSamples'});
      $n_samples_tested = $mutated_samples_info->{'numberSamplesTested'} if(defined $mutated_samples_info->{'numberSamplesTested'});
      $n_samples_mutation = $mutated_samples_info->{'numberSamplesWithMutationType'} if(defined $mutated_samples_info->{'numberSamplesWithMutationType'});
    }

    if (! defined $phenotype_id){
      print $errFH1 "phenotype_id (json_hash->{'diseaseId'}) not found for $gene_id!\n";
      next;
    }

    # Phenotype fetching and parsing from the Ontology db
    my $phenotype;
    if ($phenotype_id =~ /^EFO/) {
      $phenotype = ($efos{$phenotype_id}) ? $efos{$phenotype_id} : get_phenotype_desc($phenotype_id, $ota);
    }

    if (!$phenotype) {
      print $errFH1 "$gene_id: no phenotype desc found for $phenotype_id\n";
      $phenotype = 'ND';
    }
    else {
      $efos{$phenotype_id} = $phenotype;
    }

    my $gene = $ga->fetch_by_stable_id($gene_id, 'HGNC');		

    if(!$gene) {
      print $errFH1 "$gene_id could not be found in the Ensembl db!\n";
      next;
    }

    $data{$gene_id}{$phenotype}{'type'} = $type;
    $data{$gene_id}{$phenotype}{'source'} = $source;
    $data{$gene_id}{$phenotype}{'phenotype_id'} = $phenotype_id;
    foreach my $pmid (@$pmids) {
      $data{$gene_id}{$phenotype}{'pmids'}{$pmid} = 1; 
    }
    $data{$gene_id}{$phenotype}{'so_term'} = $so_term;
    $data{$gene_id}{$phenotype}{'n_mutated_samples'} = $n_mutated_samples;
    $data{$gene_id}{$phenotype}{'n_samples_tested'} = $n_samples_tested;
    $data{$gene_id}{$phenotype}{'n_samples_mutation'} = $n_samples_mutation;
  }
  close(IN);

  open(OUT, "> ".$self->workdir."/$output_file") || die ("Could not open file for writing: $!\n");
  foreach my $gene_id (sort(keys(%data))) {
      foreach my $phenotype (keys(%{$data{$gene_id}})) {
        my $type = $data{$gene_id}{$phenotype}{'type'};
        my $source = $data{$gene_id}{$phenotype}{'source'};
        my $phenotype_id = $data{$gene_id}{$phenotype}{'phenotype_id'};
        my $pmids = join(',', keys(%{$data{$gene_id}{$phenotype}{'pmids'}}));
        my $so_term = $data{$gene_id}{$phenotype}{'so_term'};
        my $n_mutated_samples = $data{$gene_id}{$phenotype}{'n_mutated_samples'};
        my $n_samples_tested = $data{$gene_id}{$phenotype}{'n_samples_tested'};
        my $n_samples_mutation = $data{$gene_id}{$phenotype}{'n_samples_mutation'};
        print OUT "$gene_id\t$type\t$source\t$phenotype\t$phenotype_id\t$pmids\t$so_term\t$n_mutated_samples\t$n_samples_tested\t$n_samples_mutation\n";
      }
  }

  close(OUT);

  close ($errFH1);
}


=head2 get_phenotype_desc

  Arg [1]    : string $phenotype_id
               Phenotype ontology accession.
  Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor $ota
               The ontology db term adaptor file.
  Example    : $obj->get_phenotype_desc($phenotype_id, $ota)
  Description: Fetch phenotype description based on phenotype id from Ontology database.
  Returntype : string phenotype description
  Exceptions : none

=cut

sub get_phenotype_desc {
  my $id = shift;
  my $ota = shift;

  $id =~ s/ //g;
  $id =~ s/_/:/g; #OpenTargets uses EFO_0000389 not EFO:0000389

  my $phenotype;
  my $term = $ota->fetch_by_accession($id);
  $phenotype = $term->name if ($term);

  return $phenotype;
}


=head2 parse_publications

  Arg [1]    : string $references_json
               Json references entry.
  Example    : $obj->parse_publications($references_json)
  Description: Parse the publication pmid from the OpenTargets json record.
  Returntype : arrayref
  Exceptions : none

=cut

sub parse_publications {
  my $pubs = shift;

  my @pmids = ();
  foreach my $pub_id (@$pubs) {
      my $pmid = "PMID:$pub_id";
      push @pmids,$pmid;
  }
  return \@pmids;
}


=head2 parse_input_file

  Arg [1]    : string $infile
               The input file name.
  Example    : $results = $obj->parse_input_file($infile)
  Description: Parse the specific Cancer Genome Consensus data (from OpenTargets) reformated in tabulated file.
  Returntype : hashref with results (key 'phenotypes')
  Exceptions : none

=cut

sub parse_input_file {
  my ($self, $infile) = @_;

  my $ga = $self->core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor\n") unless defined($ga);

  my $errFH1;
  open($errFH1, ">", $self->workdir."/".'log_import_err_'.$infile) || die ("Failed to open file".$self->workdir."/".'log_import_err_'.$infile.": $!\n");

  my @phenotypes;

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open(IN, "zcat ".$self->workdir."/$infile |") || die ("Could not open $infile for reading: $!\n");
  }
  else {
    open(IN,'<',$self->workdir."/".$infile) || die ("Could not open $infile for reading: $!\n");
  }

  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;

    my @row_data = split(/\t/, $_);

    # get data
    my $gene_id   = $row_data[0];
    my $phen      = $row_data[3];
    my $accession = $row_data[4];
    my $pmids     = $row_data[5];
    my $so_term   = $row_data[6];
    my $n_mutated_samples  = $row_data[7];
    my $n_samples_tested   = $row_data[8];
    my $n_samples_mutation = $row_data[9];

    ## change accession format
    $accession =~ s/\_/\:/;

    my $gene = $ga->fetch_by_stable_id($gene_id);

    if(!$gene) {
      print $errFH1 "WARNING: Ensembl gene $gene_id not found in the Ensembl Core database. Skipped!\n";
      next;
    }

    my %data = (
      'id' => $gene_id,
      'description' => $phen,
      'seq_region_id' => $gene->slice->get_seq_region_id,
      'seq_region_start' => $gene->seq_region_start,
      'seq_region_end' => $gene->seq_region_end,
      'seq_region_strand' => $gene->seq_region_strand,
      'SO_accession' => $so_term,
      'mutated_samples' => $n_mutated_samples,
      'samples_tested' => $n_samples_tested,
      'samples_mutation' => $n_samples_mutation,
      'accessions' => [$accession],
      'ontology_mapping_type' => 'is'
    );

    $data{'study'} = $pmids if ($pmids);

    push(@phenotypes,\%data);
  }
  close(IN);
  close($errFH1);

  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

1;

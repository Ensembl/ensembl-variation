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


=head1 ImportGWAS

This module imports NHGRI-EBI GWAS catalog data. The module fetches the data from
the API hosted at EBI: https://www.ebi.ac.uk/gwas/.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportGWAS;

use warnings;
use strict;

use File::Path qw(make_path);
use File::stat;
use File::Basename;
use POSIX qw(strftime);

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;

sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');
  my $run_type     = $self->required_param('run_type');

  $self->debug($self->param('debug_mode'));
  $self->skip_sets(0); #add variation set data
  $self->gwas(1); #run on gwas specific behaviour

  my $gwas_url = 'https://www.ebi.ac.uk/gwas/api/search/downloads/alternative';

  %source_info = (source_description => 'Variants associated with phenotype data from the NHGRI-EBI GWAS catalog',
                  source_url => 'https://www.ebi.ac.uk/gwas/',
                  object_type => 'Variation',
                  source_status => 'germline',
                  set => 'ph_nhgri',
                  source_name       => 'NHGRI-EBI GWAS catalog', #source name in the variation db
                  source_name_short => 'GWAS',                   #source identifier in the pipeline
                  data_types => 'phenotype_feature,study',
                # source_version is set based on file name
                  );

  my $workdir = $pipeline_dir."/".$source_info{source_name_short}."/".$species;
  make_path($workdir) or die "Failed to create $workdir $!\n";
  $self->workdir($workdir);

  open(my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open(my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);

  #get input files:
  my $file_gwas;
  # first look if one already was downloaded, if not then fetch latest
  my @files = glob($workdir. '/gwas_catalog_v*-associations_e*_r*\.tsv' );
  my ($latestFile, $latestFileTime);
  foreach my $f (@files){
    if (defined $latestFile) {
      my $fileTime = strftime("%Y%m%d", localtime(stat($f)->mtime)); #get file date
      if ($fileTime > $latestFileTime) {
        $latestFile = $f; $latestFileTime = $fileTime;
      }
    } else {
      $latestFile = $f;
      $latestFileTime = strftime("%Y%m%d", localtime(stat($latestFile)->mtime));
    }
  }
  if (defined $latestFile) {
    $file_gwas = basename($latestFile);
    print $logFH "Found file (".$workdir."/".$file_gwas.") and will skip new fetch\n";
  } else {
    chdir $workdir;
    my $redirInfo=`curl -w '%{redirect_url}' -JLO $gwas_url `;
    my @tmp = split(/'/, $redirInfo);
    $file_gwas = $tmp[1];
    print $logFH "Fetched new file : $file_gwas\n";
  }

  #parse date from file name
  if ($file_gwas  =~ m/gwas_catalog_v.*-associations_e.*_r(.*)\.tsv/){
    $source_info{source_version} = $1;
    $source_info{source_version} =~ s/-//g;
  } else {
    print $errFH "WARNING: could not parse version from file : $file_gwas\n";
  }
  $self->param('gwas_file', $file_gwas);
}

sub run {
  my $self = shift;

  my $gwas_file = $self->required_param('gwas_file');

  # dump and clean pre-existing phenotype features
  $self->dump_phenotypes($source_info{source_name}, 1);

  # get phenotype data
  my $results = $self->parse_input_file($gwas_file);
  $self->print_pipelogFH("Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n") if ($self->debug);

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
  close($self->logFH) if defined $self->logFH ;
  close($self->errFH) if defined $self->errFH ;
  close($self->pipelogFH) if defined $self->pipelogFH ;

  $self->dataflow_output_id($self->param('output_ids'), 2);
}


=head2 parse_input_file

  Arg [1]    : string $infile
               The input file name.
  Example    : $results = $obj->parse_input_file($infile)
  Description: Parse phenotypes from NHGRI-EBI GWAS input file
  Returntype : hashref with results (key 'phenotypes')
  Exceptions : none

=cut

sub parse_input_file {
  my ($self, $infile) = @_;

  my %headers;
  my @phenotypes;

  my $errFH1;
  open($errFH1, ">", $self->workdir."/".'log_import_err_'.$infile) || die ("Could not open file for writing: $!\n");

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

    my @row_data = split(/\t/,$_);

    # header
    if(/^DATE\s+ADDED\s+TO\s+CATALOG/) {
      $headers{uc($row_data[$_])} = $_ for 0..$#row_data;
    }
    else {
      die ("ERROR: Could not find header data\n") unless %headers;

      my %content;
      $content{$_} = $row_data[$headers{$_}] for keys %headers;

      my $pubmed_id      = $content{'PUBMEDID'};
      my $study          = $content{'STUDY'};
      my $phenotype      = $content{'DISEASE/TRAIT'};
      my $gene           = ($content{'REPORTED GENE(S)'} =~ /\?/) ? '' : $content{'REPORTED GENE(S)'};
      my $rs_risk_allele = ($content{'STRONGEST SNP-RISK ALLELE'} =~ /\?/) ? '' : $content{'STRONGEST SNP-RISK ALLELE'};
      my $rs_id          = $content{'SNPS'};
      my $pvalue         = ($content{'P-VALUE'} ne '') ? $content{'P-VALUE'} : '';
      my $ratio          = $content{'OR OR BETA'};
      my $ratio_info     = $content{'95% CI (TEXT)'};
      my @accessions     = split/\,/, $content{'MAPPED_TRAIT_URI'};

      print $errFH1 "WARNING: 'DISEASE/TRAIT' entry is empty for '$rs_id'\n" if ($phenotype eq '');
      next if ($phenotype eq '');

      my $risk_frequency = '';
      if ($rs_risk_allele =~ /^\s*$rs_id-+\s*(\w+)\s*$/i) {
        $rs_risk_allele = $1;
        $risk_frequency = $content{'RISK ALLELE FREQUENCY'};
      }

      $gene =~ s/\s+//g;
      $gene =~ s/–/-/g;
      $gene =~ s/[^\x00-\x7F]//g; # Remove non ASCII characters
      $gene = '' if $gene eq '-' or $gene eq 'NR'; #Skip uninformative entries, missing data in original curation see GWAS catalog curation

      my %data = (
        'study_type' => 'GWAS',
        'description' => $phenotype,
        'associated_gene' => $gene,
        'risk_allele' => $rs_risk_allele,
        'risk_allele_freq_in_controls' => $risk_frequency,
        'p_value' => $pvalue,
        'study_description' => $study,
        'accessions'   => \@accessions,
        'ontology_mapping_type' =>'is'
      );

      # Post process the ratio data
      if (defined($ratio)) {
        if ($ratio =~ /(\d+)?(\.\d+)$/) {
          my $pre  = $1;
          my $post = $2;
          $ratio = (defined($pre)) ? "$pre$post" : "0$post";
          $ratio = 0 if ($ratio eq '0.00');
        } else {
          $ratio = undef;
        }
      }

      # Add ratio/coef
      if (defined($ratio)) {
        # Parse the ratio info column to extract the unit information (we are not interested in the confidence interval)
        if ($ratio_info =~ /^\s*(\[.+\])?\s*(.+)$/) {
          my $unit = $2;
             $unit =~ s/\(//g;
             $unit =~ s/\)//g;
             $unit =~ s/µ/micro/g;
          if ($unit =~ /^\s+$/ || $ratio >= 1) {
            $data{'odds_ratio'} = $ratio;
          }
          else {
            $data{'beta_coef'} = "$ratio $unit";
          }
        }
        else {
          $data{'odds_ratio'} = $ratio;
        }
      }

      # Parse the ids
      my @ids;
      $rs_id ||= "";
      while ($rs_id =~ m/(rs[0-9]+)/g) {
        push(@ids,$1);
      }
      $data{'variation_names'} = join(',',@ids);
      $data{'study'} = $self->get_pubmed_prefix . $pubmed_id if (defined($pubmed_id));

      # If we did not get any rsIds, skip this row (this will also get rid of the header)
      print $errFH1 "WARNING: Could not parse any rsIds from string '$rs_id'\n" if (!scalar(@ids));
      next if (!scalar(@ids));

      map {
        my %t_data = %{\%data};
        $t_data{'id'} = $_;
        push(@phenotypes,\%t_data)
      } @ids;
    }
  }
  close(IN);
  close($errFH1);

  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

1;

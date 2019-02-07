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

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportGWAS;

use warnings;
use strict;

use File::Path qw(make_path);
use LWP::Simple;
use Data::Dumper; #TODO: remove if not needed
use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(GWAS NONE);
#use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation qw($variation_dba);

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

#TODO: should contain everything about RGD, including getting the file(s) and for all RGD species

my %source_info;
my $workdir;
my $core_dba;
my $variation_dba;
my $phenotype_dba;

my $pubmed_prefix = 'PMID:';

my $debug;

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
    #-source nhgri -infile ${datadir}/gwascatalog.txt -verbose -version 20121031

    my $gwas_url = 'http://www.ebi.ac.uk/gwas/api/search/downloads/alternative';

    %source_info = (source_description => 'Variants associated with phenotype data from the NHGRI-EBI GWAS catalog',
                    source_url => 'http://www.ebi.ac.uk/gwas/',
                    object_type => 'Variation',
                    set => 'ph_nhgri',

                    source_name_use => 'GWAS',
                    source_name => 'NHGRI-EBI GWAS catalog', #TODO: figure out where this is used
                    source_version => $self->required_param('nhgri_version'), #TODO: figure out how to get the original file_name and get : 'Date on the file name'

                    source => 'nhgri',

                    #  object_type => 'QTL', #default type, import code will switch to Gene for the Gene-type ones
                  #  source_mapped_attrib_type => 'Rat Genome Database', #for ontology mapping (attr_type_id 509) entry in phenotype_ontology_accession (attr_id 588)
                  #  source_status => 'germline',
                  #  threshold => $self->required_param('threshold_qtl'),
                    );

    $workdir = $pipeline_dir."/".$source_info{source_name_use}."/".$species;
    my $file_gwas = 'gwas_catalog_v1.0.2-associations.txt';

    #get input files RGD eQTL, GENE:
    make_path($workdir);
    getstore($gwas_url, $workdir."/".$file_gwas) unless -e $workdir."/".$file_gwas;
    print "Found files (".$workdir."/".$file_gwas.") and will skip new fetch\n" if -e $workdir."/".$file_gwas;
    $self->param('gwas_file', $file_gwas);
}

sub run {
  my $self = shift;

  my $gwas_file = $self->required_param('gwas_file');

  local (*STDOUT, *STDERR);
  open STDOUT, ">", $workdir."/".'log_import_out_'.$gwas_file; #TODO: what is best error/out log naming convention?
  open STDERR, ">", $workdir."/".'log_import_err_'.$gwas_file;

  # get phenotype data
  my $results = parse_nhgri($workdir."/".$gwas_file);
  print "Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n" if $debug ;

  # save phenotypes
  $self->save_phenotypes(\%source_info, $results, $core_dba, $variation_dba);

  my %param_source = (source_name => $source_info{source_name_use},
                      type => $source_info{object_type};
  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species')
                             });
}

sub write_output {
  my $self = shift;

  $self->dataflow_output_id($self->param('output_ids'), 1);
  print "Passing NHGRI-EBI GWAS import for checks\n" if $self->param('debug_mode');
}


# NHGRI-EBI GWAS specific phenotype parsing method
sub parse_nhgri {
  my $infile = shift;

  my %headers;
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

    my @row_data = split(/\t/,$_);

    # header
    if(/^DATE\s+ADDED\s+TO\s+CATALOG/) {
      $headers{uc($row_data[$_])} = $_ for 0..$#row_data;
    }
    else {
      die "ERROR: Couldn't find header data\n" unless %headers;

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

      my $risk_frequency = '';
      if ($rs_risk_allele =~ /^\s*$rs_id-+\s*(\w+)\s*$/i) {
        $rs_risk_allele = $1;
        $risk_frequency = $content{'RISK ALLELE FREQUENCY'};
      }

      $gene =~ s/\s+//g;
      $gene =~ s/–/-/g;
      $gene =~ s/[^\x00-\x7F]//g; # Remove non ASCII characters

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
      $data{'study'} = $pubmed_prefix . $pubmed_id if (defined($pubmed_id));

      # If we didn't get any rsIds, skip this row (this will also get rid of the header)
      warn("Could not parse any rsIds from string '$rs_id'\n") if (!scalar(@ids));
      next if (!scalar(@ids));

      map {
        my %t_data = %{\%data};
        $t_data{'id'} = $_;
        push(@phenotypes,\%t_data)
      } @ids;
    }
  }
  close(IN);

  my %result = ('phenotypes' => \@phenotypes);

  return \%result;
}

1;

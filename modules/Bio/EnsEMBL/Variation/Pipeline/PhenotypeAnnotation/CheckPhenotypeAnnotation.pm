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


=head1 CheckPhenotypeAnnotation

This module runs checks and produces reports after the phenotype import stage

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation);
use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(species);

my $source;
my $workdir;
my $report;

sub fetch_input {
  my $self = shift;
  $source = $self->required_param('source');
  $workdir = $self->param('workdir');
  $workdir ||= $self->required_param('pipeline_dir')."/".$source->{source_name}."/".$self->required_param('species');
}

sub run {
  my $self = shift;

  my $dbh = $self->get_species_adaptor('variation')->dbc;

  open $report, ">$workdir/REPORT_QC.txt"||die "Failed to open report file for summary info :$!\n";

  print $report "\nRunning checks on phenotype, phenotype_feature, phenotype_feature_attrib and phenotype_ontology_accession data\n";
  check_phenotype_description($dbh);
  print $report "\n";
  check_fk($dbh);

  $self->param('output_ids', { source => $self->required_param('source'),
                      species => $self->required_param('species'),
                      workdir => $workdir,
                    });
  close $report;
}

sub write_output {
  my $self = shift;

  # if species is an 'ontology' term species then move to dataflow 2
  my %import_species = &species;
  my %ontology_species = map { $_ => 1 } @{$import_species{'ontology'}};
  if (defined ($ontology_species{$self->param('species')})) {
    if ($self->param('debug_mode')) {
      open (my $logPipeFH, ">>", $workdir."/".'log_import_debug_pipe');
      print $logPipeFH "Passing $source->{source_name} import (".$self->param('species').") for adding ontology accessions (ontology_mapping)\n";
      close ($logPipeFH);
    }
    $self->dataflow_output_id($self->param('output_ids'), 2);
  } else {
    if ($self->param('debug_mode')) {
      open (my $logPipeFH, ">>", $workdir."/".'log_import_debug_pipe');
      print $logPipeFH "Passing $source->{source_name} import (".$self->param('species').") for summary counts (finish_phenotype_annotation)\n";
      close ($logPipeFH);
    }
    $self->dataflow_output_id($self->param('output_ids'), 3);
  }
}


## check for unsupported characters in phenotype names 
sub check_phenotype_description{
  my $dbh = shift;

  my $pheno_ext_sth = $dbh->prepare(qq[ select phenotype_id, description from phenotype ]);
  $pheno_ext_sth->execute()||die;

  my $ph =  $pheno_ext_sth->fetchall_arrayref();
  foreach my $l (@{$ph}){
    print $report "WARNING: Phenotype id:$l->[0] has no description!\n" unless defined $l->[1];
    next unless defined $l->[1];
    print $report "WARNING: Phenotype id:$l->[0] has empty description!\n" if $l->[1] eq '';

    my $full = $l->[1];
    $l->[1] =~ s/\w+|\-|\,|\(|\)|\s+|\/|\.|\;|\+|\'|\:|\@|\*|\%//g;
    print $report "WARNING: Phenotype : $full (id:$l->[0]) looks suspect!\n" if(length($l->[1]) >0);

    # check for characters which will be interpreted a new lines
    $l->[1] =~ /.*\n.*/;
    print $report "WARNING: Phenotype : $full (id:$l->[0]) contains a newline \n" if(length($l->[1]) >0);

    # check for phenotype descriptions suggesting no phenotype
    print $report "WARNING: Phenotype : $full (id:$l->[0]) is not useful \n" if !checkNonTerms( $l->[1] );

    # check for unsupported individual character
    print $report "WARNING: Phenotype : $full (id:$l->[0]) has suspect start or unsupported characters \n" if !checkUnsupportedChar( $l->[1] );

  }

}

sub checkNonTerms {
  my $desc = shift;

  my $is_ok = 1;
  my @junk = ("None", "Not provided", "not specified", "Not in OMIM", "Variant of unknown significance", "not_provided", "?",".", "ClinVar: phenotype not specified");

  for my $check (@junk){
    if (index($desc, $check) != -1) {
      $is_ok  = 0;
      return $is_ok;
    }
  }

  return $is_ok;
}

sub checkUnsupportedChar {
  my $desc = shift;

  my $is_ok = 1;
  my $i = 0;
  for my $c (split //, $desc) {
    # get ascii code
    my $ascii_val = ord($c);

    # check code in supported range
    if($ascii_val < 32 || $ascii_val  > 126 || $ascii_val == 60 || $ascii_val == 62 ){
      $is_ok = 0;
    }

    # also check first character makes sense
    if($i == 0 && ( $ascii_val < 48 ||
      ($ascii_val  > 57 && $ascii_val < 65) ||
      ($ascii_val  > 90 && $ascii_val < 97) ||
      $ascii_val  > 122)){
      $is_ok = 0;
      $i++;
    }
  }
  return $is_ok;
}

# basic report on imported data
sub check_fk{
  my $dbh = shift;

  # check FKs: that all phenotypes have a phenotype_feature, all pf have a pfa, and vice versa
  my $featless_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype
                                                  where phenotype_id not in (select phenotype_id from phenotype_feature) ]);
  my $attribless_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype_feature
                                                 where phenotype_feature_id not in (select phenotype_feature_id from phenotype_feature_attrib) ]);
  my $featless_attrib_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype_feature_attrib
                                                  where phenotype_feature_id not in (select phenotype_feature_id from phenotype_feature) ]);
  my $phenoless_feat_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype_feature
                                                  where phenotype_id not in (select phenotype_id from phenotype) ]);
  my $phenoless_acc_count_ext_sth = $dbh->prepare(qq[ select count(*) from phenotype_ontology_accession
                                                where phenotype_id not in (select phenotype_id from phenotype) ]);

  $featless_count_ext_sth->execute()||die;
  my $featless_count = $featless_count_ext_sth->fetchall_arrayref();
  print $report "$featless_count->[0]->[0] phenotype entries with no phenotype_feature entry (expected: 0)\n";

  $attribless_count_ext_sth->execute()||die;
  my $attrless_count = $attribless_count_ext_sth->fetchall_arrayref();
  print $report "$attrless_count->[0]->[0] phenotype_feature entries with no phenotype_feature_attrib entry (can be valid cases)\n";

  $featless_attrib_count_ext_sth->execute()||die;
  my $featless_attrib_count = $featless_attrib_count_ext_sth->fetchall_arrayref();
  print $report "$featless_attrib_count->[0]->[0] phenotype_feature_attrib entries with missing phenotype_feature entry (expected: 0)\n";

  $phenoless_feat_count_ext_sth->execute()||die;
  my $phenoless_count = $phenoless_feat_count_ext_sth->fetchall_arrayref();
  print $report "$phenoless_count->[0]->[0] phenotype_feature entries with missing phenotype entry (expected: 0)\n";

  $phenoless_acc_count_ext_sth->execute()||die;
  my $phenoless_acc_count = $phenoless_acc_count_ext_sth->fetchall_arrayref();
  print $report "$phenoless_acc_count->[0]->[0] phenotype_ontology_accession rows with missing phenotype entry (expected: 0)\n";

}

1;


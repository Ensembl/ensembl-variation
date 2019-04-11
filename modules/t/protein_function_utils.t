# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


use strict;
use warnings;
use Data::Dumper;
use Test::More;
use Test::Deep;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::Utils::DbNSFPProteinFunctionAnnotation;
use Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation;
our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);

my $dir = $multi->curr_dir();
my $dbNSFP = Bio::EnsEMBL::Variation::Utils::DbNSFPProteinFunctionAnnotation->new(
  -working_dir => $dir,
  -species => 'Homo_sapiens',
  -annotation_file => $dir . '/testdata/dbNSFP3.5a_grch37.txt.gz',
  -assembly => 'GRCh37',
  -annotation_file_version => '3.5a',
  -pipeline_mode => 1,
  -write_mode => 0,
  -debug_mode => 1,
);

ok($dbNSFP->working_dir eq $dir, 'working_dir');

$dbNSFP->run('4d08f77b4cb14259684ce086ba089565', {'4d08f77b4cb14259684ce086ba089565' => 'ENSP00000435699'});
my @results = ();
my $pfpma = $vdb->get_ProteinFunctionPredictionMatrixAdaptor;
my $pred_matrices = $dbNSFP->{pred_matrices};
foreach my $analysis (keys %$pred_matrices) {
  my $pred_matrix = $pred_matrices->{$analysis};
  if ($dbNSFP->{results_available}->{$analysis}) {
    my $matrix = $pfpma->fetch_by_analysis_translation_md5($analysis, '4d08f77b4cb14259684ce086ba089565');
    my $debug_data = $dbNSFP->{debug_data};
    my $i = 588;
    my $aa = 'S';
    foreach my $prediction (keys %{$debug_data->{$analysis}->{$i}->{$aa}}) {
      my ($new_pred, $new_score) = $matrix->get_prediction($i, $aa);
      push @results, {position => $i, aa => $aa, new_pred => $new_pred, new_score => $new_score, analysis => $analysis}; 
    }
  }
}
my @sorted_results = sort {$a->{new_score} <=> $b->{new_score}} @results;
my $expected_results = [
          {
            'new_score' => '0.033',
            'new_pred' => 'tolerated',
            'position' => 588,
            'aa' => 'S',
            'analysis' => 'dbnsfp_meta_lr'
          },
          {
            'new_score' => '0.301',
            'new_pred' => 'likely benign',
            'position' => 588,
            'aa' => 'S',
            'analysis' => 'dbnsfp_revel'
          },
          {
            'new_score' => '0.552',
            'new_pred' => 'medium',
            'position' => 588,
            'aa' => 'S',
            'analysis' => 'dbnsfp_mutation_assessor'
          },
        ];
 
cmp_deeply(\@sorted_results, $expected_results, "dbNSFP - Retrieve scores and predictions from stored matrices.");

my $amino_acids = join('', @{$dbNSFP->amino_acids}[0..9]);
ok($amino_acids eq 'MPIGSKERPT', 'Compare first 10 amino acids');

my $all_triplets = $dbNSFP->get_triplets('ENSP00000435699'); # on forward strand
my $first_triplet = $all_triplets->[0];
my $last_triplet = $all_triplets->[$#$all_triplets];
my $expected_first_triplet = {
          'triplet_seq' => 'ATG',
          'coords' => [
                        [
                          32890598,
                          32890600
                        ]
                      ],
          'aa_position' => 1,
          'chrom' => '13',
          'new_triplets' => {
                              'ATG' => {
                                         '1' => {
                                                  'A' => 'AAG',
                                                  'T' => 'ATG',
                                                  'C' => 'ACG',
                                                  'G' => 'AGG'
                                                },
                                         '0' => {
                                                  'A' => 'ATG',
                                                  'T' => 'TTG',
                                                  'C' => 'CTG',
                                                  'G' => 'GTG'
                                                },
                                         '2' => {
                                                  'A' => 'ATA',
                                                  'T' => 'ATT',
                                                  'C' => 'ATC',
                                                  'G' => 'ATG'
                                                }
                                       }
                            },
          'aa' => 'M'
        };

my $expected_last_triplet =  {
          'triplet_seq' => 'AAA',
          'coords' => [
                        [
                          32907425,
                          32907427
                        ]
                      ],
          'aa_position' => 602,
          'chrom' => '13',
          'new_triplets' => {
                              'AAA' => {
                                         '1' => {
                                                  'A' => 'AAA',
                                                  'T' => 'ATA',
                                                  'C' => 'ACA',
                                                  'G' => 'AGA'
                                                },
                                         '0' => {
                                                  'A' => 'AAA',
                                                  'T' => 'TAA',
                                                  'C' => 'CAA',
                                                  'G' => 'GAA'
                                                },
                                         '2' => {
                                                  'A' => 'AAA',
                                                  'T' => 'AAT',
                                                  'C' => 'AAC',
                                                  'G' => 'AAG'
                                                }
                                       }
                            },
          'aa' => 'K'
        };

cmp_deeply($last_triplet, $expected_last_triplet, "last triplet structure");
cmp_deeply($first_triplet, $expected_first_triplet, "first triplet structure");
my @rows = ();
my $iter = $dbNSFP->get_tabix_iterator(13, 32907425, 32907425);
while (my $line = $iter->next) {
  push @rows, $dbNSFP->get_dbNSFP_row($line);
}
my $expected_rows =  [
          {
            'alt' => 'C',
            'chr' => '13',
            'ref' => 'A',
            'aaref' => 'K',
            'mutation_assessor_pred' => 'M',
            'mutation_assessor_score' => '0.53678',
            'refcodon' => 'AAA',
            'meta_lr_pred' => 'T',
            'pos' => '32907425',
            'meta_lr_score' => '0.0099',
            'aaalt' => 'Q',
            'revel_score' => '0.121'
          },
          {
            'alt' => 'G',
            'chr' => '13',
            'ref' => 'A',
            'aaref' => 'K',
            'mutation_assessor_pred' => 'L',
            'mutation_assessor_score' => '0.32339',
            'refcodon' => 'AAA',
            'meta_lr_pred' => 'T',
            'pos' => '32907425',
            'meta_lr_score' => '0.0086',
            'aaalt' => 'E',
            'revel_score' => '0.166'
          },
          {
            'alt' => 'T',
            'chr' => '13',
            'ref' => 'A',
            'aaref' => 'K',
            'mutation_assessor_pred' => '.',
            'mutation_assessor_score' => '.',
            'refcodon' => 'AAA',
            'meta_lr_pred' => '.',
            'pos' => '32907425',
            'meta_lr_score' => '.',
            'aaalt' => 'X',
            'revel_score' => '.'
          }
        ];

cmp_deeply(\@rows, $expected_rows, "dbNSFP rows for location 13:32907425-32907425");

my $cadd = Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation->new(
  -species => 'Homo_sapiens',
  -annotation_file => $dir . '/testdata/cadd_v1.3_grch37.txt.gz',
  -assembly => 'GRCh37',
  -annotation_file_version => 'v1.3',
  -pipeline_mode => 1,
  -write_mode => 0,
  -debug_mode => 1,
);
$cadd->run('4d08f77b4cb14259684ce086ba089565', {'4d08f77b4cb14259684ce086ba089565' => 'ENSP00000435699'});
@results = ();
$pred_matrices = $cadd->{pred_matrices};
foreach my $analysis (keys %$pred_matrices) {
  my $pred_matrix = $pred_matrices->{$analysis};
  if ($cadd->{results_available}->{$analysis}) {
    my $matrix = $pfpma->fetch_by_analysis_translation_md5($analysis, '4d08f77b4cb14259684ce086ba089565');
    my $debug_data = $cadd->{debug_data};
    my $i = 588;
    my $aa = 'S';
    foreach my $prediction (keys %{$debug_data->{$analysis}->{$i}->{$aa}}) {
      my ($new_pred, $new_score) = $matrix->get_prediction($i, $aa);
      push @results, {position => $i, aa => $aa, new_pred => $new_pred, new_score => $new_score, analysis => $analysis}; 
    }
  }
}

$expected_results =  [
          {
            'new_score' => 28,
            'new_pred' => 'likely benign',
            'position' => 588,
            'aa' => 'S',
            'analysis' => 'cadd'
          }
        ];

cmp_deeply(\@results, $expected_results, "CADD - Retrieve scores and predictions from stored matrices.");


@rows = ();
$iter = $cadd->get_tabix_iterator(13, 32907425, 32907425);
while (my $line = $iter->next) {
  push @rows, $cadd->get_CADD_row($line);
}

$expected_rows =  [
          {
            'RawScore' => '1.819199',
            '#Chrom' => '13',
            'Ref' => 'A',
            'PHRED' => '15.10',
            'Pos' => '32907425',
            'Alt' => 'C'
          },
          {
            'RawScore' => '1.149219',
            '#Chrom' => '13',
            'Ref' => 'A',
            'PHRED' => '11.48',
            'Pos' => '32907425',
            'Alt' => 'G'
          },
          {
            'RawScore' => '9.837337',
            '#Chrom' => '13',
            'Ref' => 'A',
            'PHRED' => '36',
            'Pos' => '32907425',
            'Alt' => 'T'
          }
        ];

cmp_deeply(\@rows, $expected_rows, "CADD rows for location 13:32907425-32907425");

done_testing();


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
my $pfpma = $vdb->get_ProteinFunctionPredictionMatrixAdaptor;

$vdb->dnadb($db);

my $tr_ad  = $db->get_TranscriptAdaptor;
my $transcript_stable_id = 'ENST00000238738';
my $transcript = $tr_ad->fetch_by_stable_id($transcript_stable_id); #get the Transcript object
my $translation = $transcript->translation;
my $dir = $multi->curr_dir();

# Test reading from dbNSFP file for different dbNSFP file versions and assemblies:
# The tests use a codon from ENST00000238738 and all the mutated versions of that codon
# The results are the predicitions from dbNSFP for each of the resulting amino acid changes
my $coords = {
  'GRCh37' => [[ 46770246, 46770248 ]],
  'GRCh38' => [[ 46543107, 46543109 ]],
};
my $triplet = {
  'triplet_seq' => 'GGC',
  'aa_position' => 21,
  'chrom' => '2',
  'new_triplets' => {
    'GGC' => {
      '0' => {
        'A' => 'AGC',
        'T' => 'TGC',
        'C' => 'CGC',
        'G' => 'GGC'
      },
      '1' => {
        'A' => 'GAC',
        'T' => 'GTC',
        'C' => 'GCC',
        'G' => 'GGC'
      },
      '2' => {
        'A' => 'GGA',
        'T' => 'GGT',
        'C' => 'GGC',
        'G' => 'GGG'
      }
    }
  },
  'aa' => 'G'
};

my @file_versions = ('3.5a', '4.0a', '4.1a', '4.2a', '4.3a', '4.4a', '4.5c', '4.6c', '4.7c', '4.8c', '4.9c', '4.9a');
my @assemblies = ('GRCh37', 'GRCh38');
my $expected_results = {
  'GRCh37' => {

    '3.5a' => { 'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99261' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},

    '4.0a' => {'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},

    '4.1a' => {'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},

    '4.2a' => {'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},

    '4.3a' => {'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},

    '4.4a' => {'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},
  
    '4.5c' => { 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
    '4.6c' => { 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
    '4.7c' => { 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
    '4.8c' => { 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
    '4.9c' => { 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
    '4.9a' => { 'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
  },

  'GRCh38' => {

    '3.5a' => { 'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99261' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},

    '4.0a' => {'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},

    '4.1a' => {'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},

    '4.2a' => {'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},

    '4.3a' => {'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},

    '4.4a' => {'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}},'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}}},

    '4.5c' => { 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
    '4.6c' => { 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
    '4.7c' => { 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
    '4.8c' => { 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
    '4.9c' => { 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
    '4.9a' => { 'dbnsfp_revel' => { '21' => { 'A' => { 'likely disease causing' => '0.940' }, 'S' => { 'likely disease causing' => '0.943' }, 'D' => { 'likely disease causing' => '0.908' }, 'C' => { 'likely disease causing' => '0.901' }, 'R' => { 'likely disease causing' => '0.944' }, 'V' => { 'likely disease causing' => '0.965' }}}, 'dbnsfp_mutation_assessor' => { '21' => { 'A' => { 'high' => '0.99996' }, 'S' => { 'high' => '0.99244' }, 'D' => { 'high' => '0.99996' }, 'C' => { 'high' => '0.99996' }, 'R' => { 'high' => '0.99996' }, 'V' => { 'high' => '0.99996' }}}, 'dbnsfp_meta_lr' => { '21' => { 'A' => { 'damaging' => '0.9814' },'S' => {'damaging' => '0.9798'},'D' => {'damaging' => '0.9836'},'C' => {'damaging' => '0.9836'},'R' => {'damaging' => '0.9814'},'V' => {'damaging' => '0.9814'}}},},
  }
};

foreach my $file_version (@file_versions) {
  foreach my $assembly (@assemblies) {
    my %dbNSFP_params = (
      -working_dir => $dir,
      -species => 'Homo_sapiens',
      -annotation_file => "$dir/testdata/test_data_dbNSFP$file_version\_$assembly.txt.gz",
      -assembly => $assembly,
      -annotation_file_version => $file_version,
      -pipeline_mode => 1,
      -write_mode => 0,
      -debug_mode => 1,
    );

    my $dbNSFP = Bio::EnsEMBL::Variation::Utils::DbNSFPProteinFunctionAnnotation->new(%dbNSFP_params);
    $dbNSFP->init_protein_matrix($translation, 'md5');
    $dbNSFP->init_header();
    $triplet->{'coords'} = $coords->{$assembly};
    $dbNSFP->load_predictions_for_triplets([$triplet]);
    cmp_deeply($dbNSFP->{debug_data}, $expected_results->{$assembly}->{$file_version}, "dbNSFP - Read from file and store predictions: $assembly, $file_version");
  }
}

my %dbNSFP_params = (
  -working_dir => $dir,
  -species => 'Homo_sapiens',
  -annotation_file => $dir . '/testdata/test_data_dbNSFP3.5a_GRCh37.txt.gz',
  -assembly => 'GRCh37',
  -annotation_file_version => '3.5a',
  -pipeline_mode => 1,
  -write_mode => 0,
  -debug_mode => 1,
);

my $dbNSFP = Bio::EnsEMBL::Variation::Utils::DbNSFPProteinFunctionAnnotation->new(%dbNSFP_params);
ok($dbNSFP->working_dir eq $dir, 'working_dir');
$dbNSFP->run('19c20411fb8a4a65deae8f1492bae6d4', {'19c20411fb8a4a65deae8f1492bae6d4' => 'ENSP00000436292'});
my $debug_data = get_debug_data($dbNSFP, '19c20411fb8a4a65deae8f1492bae6d4', 22, 'M');
$expected_results = [
  {
    'new_score' => '0.215',
    'new_pred' => 'tolerated',
    'position' => 22,
    'aa' => 'M',
    'analysis' => 'dbnsfp_meta_lr'
  },
  {
    'new_score' => '0.542',
    'new_pred' => 'likely disease causing',
    'position' => 22,
    'aa' => 'M',
    'analysis' => 'dbnsfp_revel'
  },
  {
    'new_score' => '0.648',
    'new_pred' => 'medium',
    'position' => 22,
    'aa' => 'M',
    'analysis' => 'dbnsfp_mutation_assessor'
  }
];
cmp_deeply($debug_data, $expected_results, "dbNSFP - Retrieve scores and predictions from stored matrices. Protein on reverse strand.");
my $amino_acids = join('', @{$dbNSFP->amino_acids}[0..9]);
ok($amino_acids eq 'MCAYPGCNKR', 'Compare first 10 amino acids ' . $amino_acids);

$dbNSFP = Bio::EnsEMBL::Variation::Utils::DbNSFPProteinFunctionAnnotation->new(%dbNSFP_params);
$dbNSFP->run('4d08f77b4cb14259684ce086ba089565', {'4d08f77b4cb14259684ce086ba089565' => 'ENSP00000435699'});
$debug_data = get_debug_data($dbNSFP, '4d08f77b4cb14259684ce086ba089565', 588, 'S');
$expected_results = [
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
 
cmp_deeply($debug_data, $expected_results, "dbNSFP - Retrieve scores and predictions from stored matrices.");
$amino_acids = join('', @{$dbNSFP->amino_acids}[0..9]);
ok($amino_acids eq 'MPIGSKERPT', 'Compare first 10 amino acids ' . $amino_acids);

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


$all_triplets = $dbNSFP->get_triplets('ENSP00000436292'); # on reverse strand
$first_triplet = $all_triplets->[0];
$last_triplet = $all_triplets->[$#$all_triplets];

$expected_first_triplet = {
  'triplet_seq' => 'ATG',
  'coords' => [
                [
                  32417876,
                  32417878
                ]
              ],
  'aa_position' => 1,
  'chrom' => '11',
  'new_triplets' => {
                      'ATG' => {
                                 '1' => {
                                          'A' => 'ATG',
                                          'T' => 'AAG',
                                          'C' => 'AGG',
                                          'G' => 'ACG'
                                        },
                                 '0' => {
                                          'A' => 'TTG',
                                          'T' => 'ATG',
                                          'C' => 'GTG',
                                          'G' => 'CTG'
                                        },
                                 '2' => {
                                          'A' => 'ATT',
                                          'T' => 'ATA',
                                          'C' => 'ATG',
                                          'G' => 'ATC'
                                        }
                               }
                    },
  'aa' => 'M'
};
$expected_last_triplet = {
  'triplet_seq' => 'ACT',
  'coords' => [
                [
                  32417804,
                  32417806
                ]
              ],
  'aa_position' => 25,
  'chrom' => '11',
  'new_triplets' => {
                      'ACT' => {
                                 '1' => {
                                          'A' => 'ATT',
                                          'T' => 'AAT',
                                          'C' => 'AGT',
                                          'G' => 'ACT'
                                        },
                                 '0' => {
                                          'A' => 'TCT',
                                          'T' => 'ACT',
                                          'C' => 'GCT',
                                          'G' => 'CCT'
                                        },
                                 '2' => {
                                          'A' => 'ACT',
                                          'T' => 'ACA',
                                          'C' => 'ACG',
                                          'G' => 'ACC'
                                        }
                               }
                    },
  'aa' => 'T'
};

cmp_deeply($first_triplet, $expected_first_triplet, "first triplet structure on reverse strand");
cmp_deeply($last_triplet, $expected_last_triplet, "last triplet structure on reverse strand");

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


my %cadd_params = ( 
  -species => 'Homo_sapiens',
  -annotation_file => $dir . '/testdata/test_data_cadd_v1.3_grch37.txt.gz',
  -assembly => 'GRCh37',
  -annotation_file_version => 'v1.3',
  -pipeline_mode => 1,
  -write_mode => 0,
  -debug_mode => 1,
);

my $cadd = Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation->new(%cadd_params);
$cadd->run('4d08f77b4cb14259684ce086ba089565', {'4d08f77b4cb14259684ce086ba089565' => 'ENSP00000435699'});
$debug_data = get_debug_data($cadd, '4d08f77b4cb14259684ce086ba089565', 588, 'S');
$expected_results =  [
          {
            'new_score' => 28,
            'new_pred' => 'likely benign',
            'position' => 588,
            'aa' => 'S',
            'analysis' => 'cadd'
          }
        ];
cmp_deeply($debug_data, $expected_results, "CADD - Retrieve scores and predictions from stored matrices.");

$cadd = Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation->new(%cadd_params);
$cadd->run('19c20411fb8a4a65deae8f1492bae6d4', {'19c20411fb8a4a65deae8f1492bae6d4' => 'ENSP00000436292'});
$debug_data = get_debug_data($cadd, '19c20411fb8a4a65deae8f1492bae6d4', 22, 'M');
$expected_results = [
  {
    'new_score' => 34,
    'new_pred' => 'likely deleterious',
    'position' => 22,
    'aa' => 'M',
    'analysis' => 'cadd'
  }
];
cmp_deeply($debug_data, $expected_results, "CADD - Retrieve scores and predictions from stored matrices. Protein on reverse strand.");

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

sub get_debug_data {
  my $annotation_source = shift;
  my $translation_md5 = shift;
  my $i = shift;
  my $aa = shift;
  my @results = ();
  my $pred_matrices = $annotation_source->{pred_matrices};
  foreach my $analysis (keys %$pred_matrices) {
    my $pred_matrix = $pred_matrices->{$analysis};
    if ($annotation_source->{results_available}->{$analysis}) {
      my $matrix = $pfpma->fetch_by_analysis_translation_md5($analysis, $translation_md5);
      my $debug_data = $annotation_source->{debug_data};
      if (defined $i && defined $aa) {
        foreach my $prediction (keys %{$debug_data->{$analysis}->{$i}->{$aa}}) {
          my ($new_pred, $new_score) = $matrix->get_prediction($i, $aa);
          push @results, {position => $i, aa => $aa, new_pred => $new_pred, new_score => $new_score, analysis => $analysis}; 
        }
      } else {
        foreach my $i (sort keys %{$debug_data->{$analysis}}) {
          foreach my $aa (keys %{$debug_data->{$analysis}->{$i}}) {
            foreach my $prediction (keys %{$debug_data->{$analysis}->{$i}->{$aa}}) {
              my ($new_pred, $new_score) = $matrix->get_prediction($i, $aa);
              push @results, {position => $i, aa => $aa, new_pred => $new_pred, new_score => $new_score, analysis => $analysis}; 
            }
          }
        }
      }
    }
  }
  my @sorted_results = sort {$a->{new_score} <=> $b->{new_score}} @results;
  return \@sorted_results;
}


done_testing();

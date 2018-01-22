#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

#### Check genomic, coding & non-coding HGVS strings give variation features which return the same HGVS strings
#### This is gene-annotation dependant
#### Exceptions should be reported when HGVS protein nomenclature cannot be reliably converted to genomic 

#### HGVS requires sequence accessions and versions throughout 
####   - we fail over to seq name where no acc/ver available; tested by X data here

use Test::More;

use Data::Dumper;

use FindBin qw($Bin);


use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdba = $multi->get_DBAdaptor('variation');
my $cdba = $multi->get_DBAdaptor('core');

use_ok('Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor');

my $DEBUG = 0;

## 3 pairs of hashes of equivalent names (some correct, some not) and output annotation
## Annotation format:  key =>  [hgvs_genomic,  variant_allele, hgvs_[non]coding,  variant_allele, hgvs_protein,  test_description, ref allele]

my %test_output = (
     1 => ["NC_000002.11:g.46746465G>A",
           "A",
           "ENST00000522587.1:c.-101-6514C>T", 
           "T",
           "", 
           "substitution, coding intron - downstream"
          ], 

     2 => ["NC_000002.11:g.46739156C>T", 
           "T",
           "ENST00000306448.4:c.*14G>A",
           "A",
           "",
           "substitution 3' UTR"
          ],

     3 => ["NC_000002.11:g.46746256G>A", 
           "A",
           "ENST00000306448.4:c.-274C>T", 
           "T",
           "",
           "substitution, 5' UTR" 
          ], 

     4 => ["NC_000002.11:g.46746507delinsACAA",
           "ACAA", 
           "ENST00000524249.1:n.775+16445delinsTTGT", 
           "TTGT",
           "", 
           "delins, non-coding ",
           "C"
          ], 
     5 => ["NC_000002.11:g.46739212C>G", 
           "G",
           "ENST00000522587.1:c.639G>C", 
           "C",
           "ENSP00000428141.1:p.Met213Ile", 
           "substitution, non_syn"
          ],
     6 => ["NC_000002.11:g.46739488G>A", 
           "A",
           "ENST00000306448.4:c.363C>T", 
           "T",
           "ENSP00000304891.4:p.Leu121=",
           "substitution, synonymous"
          ],
     7 => ["NC_000002.11:g.46731836A>G",   
           "G",
           "",                         
           "",
           "",
           "substitution, downstream" 
          ], 
     8 => ["NC_000002.11:g.46747460C>G", 
           "G",
           "",
           "",
           "",
           "substitution, upstream" 
          ], 
     9 => ["NC_000002.11:g.46732522G>A", 
           "A",
           "ENST00000524249.1:n.776-14552C>T",
           "T",
           "",
           "substitution, noncoding intron"
          ],
     10 => ["NC_000002.11:g.98275102C>T",
            "T",
            "ENST00000289228.5:c.445G>A",
            "A",
            "ENSP00000289228.5:p.Ala149Thr",
            "parseable protein change [-1]" 
           ], 
     11 => ["NC_000003.11:g.10191482_10191483insTTT",
            "TTT",
            "ENST00000345392.2:c.352_353insTTT",
            "TTT",
            "ENSP00000344757.2:p.Lys118delinsIleTer",
            "del ins, stop_gained",
           ],
 
     12 => ["NC_000004.11:g.41993003G>A", 
            "A", 
            "ENST00000264451.6:c.109+226G>A",   
            "A",
            "", 
            "substitution, coding intron upstream"
           ],

     13 => ["NC_000004.11:g.130032945A>G",
            "G",
            "ENST00000281146.4:c.599A>G",
            "G",
            "ENSP00000281146.4:p.Ter200TrpextTer2",
            "substitution, stop lost"
           ], 
     14 => ["NC_000004.11:g.130032948A>C",
            "C",
            "ENST00000281146.4:c.*2A>C",
            "C",
            "",
            "substitution, after stop "
            ], 
     15 => ["NC_000005.9:g.96232565_96232566insCC",
            "CC",
            "ENST00000508077.1:c.488_489insCC",
            "CC",
            "",
            "insertion, partial codon"
           ],
     16 => ["NC_000006.11:g.6649978_6649980dup",
           "AGG",
           "ENST00000230568.3:c.405+68_405+70dup",
           "AGG",
            "",
            "duplication, intronic - long",
            "AGG",
           ],
     17 =>  ["NC_000006.11:g.31997361del",
            "-",
            "ENST00000435363.2:c.3695del", 
            "-",
            "ENSP00000415941.2:p.Ser1232Ter",
            "deletion, stop gained",
            "C"
           ], 
     18 => ["NC_000007.13:g.143557504del", 
            "-",
            "ENST00000355951.2:c.1964del", 
            "-",
            "ENSP00000348220.2:p.Gln655ArgfsTer17", 
            "deletion, frameshift",
            "A"
           ],
     19 => ["NC_000007.13:g.7680048A>G",
            "G",
            "ENST00000223129.4:c.2T>C",
            "C",
            "ENSP00000223129.4:p.Met1?",
            "substitution,  start loss"
           ],

     20 => ["NC_000007.13:g.143175210G>A",
            "A",
            "ENST00000408916.1:c.245G>A",
            "A",
            "ENSP00000386201.1:p.Arg82Gln",
            "parseable protein change"
           ], 

     21 => ["NC_000012.11:g.102056227G>A",
            "A",
            "ENST00000360610.2:c.2049G>A",
            "A",
            "ENSP00000353822.2:p.Lys683=",
            "substitution synonymous"
           ], 
 
     22 => ["NC_000013.10:g.51519667dup",   #rs17857128
            "A",
            "ENST00000336617.2:c.615dup",
            "A",
            "ENSP00000337623.2:p.Glu206ArgfsTer13",
            "duplication, frameshift",
            "A"
           ],
 
     23 =>  ["NC_000017.10:g.7123233_7123234insCAGGACGTGGGCGTG",
            "CAGGACGTGGGCGTG",
            "ENST00000356839.4:c.68_69insCAGGACGTGGGCGTG",
            "CAGGACGTGGGCGTG",  
            "ENSP00000349297.4:p.Pro23_Gly24insArgThrTrpAlaTer",
            "insertion,  stop gained"
            ],

     24 => ["NC_000017.10:g.48452979_48452980insAGC",             ##rs67225428
            "AGC",
            "ENST00000393271.1:c.410_411insAGC",
            "AGC",
            "ENSP00000376952.1:p.Lys137_Pro138insAla", 
            "insertion,  codon gained"
            ], 

     25 =>  ["NC_000019.9:g.7706085T>C",                       ## rs144546645
            "C",
            "ENST00000320400.4:c.1186T>C",
            "C",
            "ENSP00000318233.4:p.Ter396GlnextTer?",
            "substitution, stop loss, no alt stop"
           ],

     26 => ["NC_000022.10:g.20920895_20920939dup",
            "CCACAGCCTCCGCCCTCCCAGGCTCTGCCCCAGCAGCTGCAGCAG",
            "ENST00000292733.7:c.832_876dup",
            "CCACAGCCTCCGCCCTCCCAGGCTCTGCCCCAGCAGCTGCAGCAG",
            "ENSP00000292733.7:p.Pro278_Gln292dup", 
            "insertion, peptide duplication",
            "CCACAGCCTCCGCCCTCCCAGGCTCTGCCCCAGCAGCTGCAGCAG",
           ],

     27 => ["MT:m.6721T>C",    ## rs199476127
            "C",
            "",
            "C",
            "", 
            "mitochondrial"
           ],

     35 => ["NC_000003.11:g.10191479_10191483delinsTTTTT",
            "TTTTT",
            "ENST00000345392.2:c.349_353delinsTTTTT",
            "TTTTT",
            "ENSP00000344757.2:p.Leu117_Lys118delinsPheLeu",
            "del ins if > 1 aa not a subs",
            "CTGAA"
           ],
     36 => ["NC_000019.9:g.7706085T>A",
            "A",
            "ENST00000320400.4:c.1186T>A",
            "A",
            "",
            "convert single base inv to subs"
           ],


);

my %test_output_shifted = (   
     28 => ["X:g.131215401dup",
            "A",
            "ENST00000298542.4:c.905+998dup", 
            "T",
            "",
            "duplication, intronic rc transcript",
            "T",
           ], 

     29 => ["NC_000011.9:g.32417913_32417914insCCTACGAGTACTACC",
            "ACCCCTACGAGTACT", 
            "ENST00000530998.1:c.454_455insAGTACTCGTAGGGGT",
            "AGTACTCGTAGGGGT", 
            "ENSP00000435307.1:p.Arg151_Ser152insTer",
            "insertion, stop gained [-1]",
            ],

     30 => ["NC_000013.10:g.51519669dup",
           "G",
           "ENST00000336617.2:c.616+1dup",
           "G",
           "",
           "insertion, frameshift lost on 3'shift",
           "G"
           ],

     31 => ["NC_000006.11:g.30558478dup",
            "A",
            "ENST00000396515.3:c.717dup",
            "A",
            "",
            "insertion, stop retained if not shifted",
            "A"
           ],

     32 => ["NC_000001.10:g.154140414_154140416del", 
            "-", 
            "ENST00000368530.2:c.857_*1del",
            "-",
            "",
            "deletion, stop loss unless shifted",
            "AAT"
           ],

     33 => ["NC_000012.11:g.102061071dup",
            "T",
            "ENST00000360610.2:c.2336-439dup",
            "T",
            "",
            "insertion, coding intron downstream",
            "NC_000012.11:g.102061071dup",
            "T"
           ],
     34 => ["NC_000019.9:g.48836480_48836482del", 
            "-",
            "ENST00000293261.2:c.1376_1378del",
            "-",
            "ENSP00000293261.2:p.Ser459del",
            "deletion, inframe codon loss",
            "CCT"
           ],

);


my %test_input = ( 
     1 => ["NC_000002.11:g.46746465G>A", 
           "ENST00000522587.1:c.-101-6514C>T", 
          ], 

     2 => ["NC_000002.11:g.46739156C>T",    
           "ENST00000306448.4:c.*14G>A",       
          ],
 
     3 => ["NC_000002.11:g.46746256G>A",    
           "ENST00000306448.4:c.-274C>T",      
          ],  

     4 => ["NC_000002.11:g.46746507delinsACAA",
           "NC_000002.11:g.46746507delCinsACAA",
           "ENST00000524249.1:n.775+16445delinsTTGT",
           "ENST00000524249.1:n.775+16445delGinsTTGT" 
          ], 
     5 => ["NC_000002.11:g.46739212C>G", 
           "ENST00000522587.1:c.639G>C",  
           "ENSP00000428141.1:p.Met213Ile",      
          ],
     6 => ["NC_000002.11:g.46739488G>A",    
           "ENST00000306448.4:c.363C>T",  
           "ENSP00000304891.4:p.Leu121=", 
           "ENST00000306448.4:c.363C>T(p.=)",    
          ],
     7 => ["NC_000002.11:g.46731836A>G",   
          ],
     8 => ["NC_000002.11:g.46747460C>G", 
          ],
     9 => ["NC_000002.11:g.46732522G>A", 
           "ENST00000524249.1:n.776-14552C>T",
          ],
     10 => ["NC_000002.11:g.98275102C>T",
            "ENST00000289228.5:c.445G>A",
            "ENSP00000289228.5:p.Ala149Thr",
           ],
     11 => ["NC_000003.11:g.10191482_10191483insTTT",
            "ENST00000345392.2:c.352_353insTTT",
            "ENSP00000344757.2:p.Lys118delinsIleTer",
           ],
 
     12 => ["NC_000004.11:g.41993003G>A", 
            "ENST00000264451.6:c.109+226G>A",   
           ],

     13 => ["NC_000004.11:g.130032945A>G",
            "ENST00000281146.4:c.599A>G",
            "ENSP00000281146.4:p.Ter200TrpextTer2",
           ], 
     14 => ["NC_000004.11:g.130032948A>C",
            "ENST00000281146.4:c.*2A>C",
           ], 
     15 => ["NC_000005.9:g.96232565_96232566insCC",
            "ENST00000508077.1:c.488_489insCC",
           ],
     16 => ["NC_000006.11:g.6649978_6649980dup",
            "ENST00000230568.3:c.405+68_405+70dup",
           ],
     17 => ["NC_000006.11:g.31997361del",
            "NC_000006.11:g.31997361delC",
            "ENST00000435363.2:c.3695del",           
            "ENST00000435363.2:c.3695delC",
            "ENSP00000415941.2:p.Ser1232Ter",
           ],    
     18 => ["NC_000007.13:g.143557504del",
            "NC_000007.13:g.143557504delT",   
            "NC_000007.13:g.143557504insTdelTT",
            "ENST00000355951.2:c.1964del", 
            "ENST00000355951.2:c.1964delA",  
            "ENSP00000348220.2:p.Gln655ArgfsTer17",  
           ],                     
     19 => ["NC_000007.13:g.7680048A>G",
            "ENST00000223129.4:c.2T>C",
            "ENSP00000223129.4:p.Met1?",
           ],

     20 => ["NC_000007.13:g.143175210G>A",
            "ENST00000408916.1:c.245G>A",
            "ENSP00000386201.1:p.Arg82Gln",
           ],  

     21 => ["NC_000012.11:g.102056227G>A",
            "ENST00000360610.2:c.2049G>A",
            "ENST00000360610.2:c.2049G>A(p.=)",
            "ENSP00000353822.2:p.Lys683=",
           ],          
       
     22 => ["NC_000013.10:g.51519667dup",  
            "NC_000013.10:g.51519667dupA",
            "ENST00000336617.2:c.615dup",
            "ENST00000336617.2:c.615dupA",
            "ENSP00000337623.2:p.Glu206ArgfsTer13",
           ],
               
     23 => ["NC_000017.10:g.7123233_7123234insCAGGACGTGGGCGTG",
            "ENST00000356839.4:c.68_69insCAGGACGTGGGCGTG",
            "ENSP00000349297.4:p.Pro23_Gly24insArgThrTrpAlaTer",
           ],
     24 => ["NC_000017.10:g.48452979_48452980insAGC",
            "ENST00000393271.1:c.410_411insAGC",
            "ENSP00000376952.1:p.Lys137_Pro138insAla", 
           ],          
     25 => ["NC_000019.9:g.7706085T>C",
           ],

     35 => ["NC_000003.11:g.10191479_10191483delinsTTTTT",
            "NC_000003.11:g.10191478_10191483delinsTTTTTT",
            "ENST00000345392.2:c.349_353delinsTTTTT",
            "ENSP00000344757.2:p.Leu117_Lys118delinsPheLeu",
            "ENSP00000344757.2:p.LeuLys117PheLeu"
            ],
     36 => ["NC_000019.9:g.7706085invT",
            "NC_000019.9:g.7706085T>A",
            "ENST00000320400.4:c.1186T>A",
           ],

);

my %test_input_shifted = (
     28 => ["X:g.131215392_131215393insA",
            "ENST00000298542.4:c.905+998dup", 
            "X:g.131215401dup",
           ],    

     29 => ["NC_000011.9:g.32417910_32417911insACCCCTACGAGTACT",
            "ENST00000530998.1:c.454_455insAGTACTCGTAGGGGT",
            "ENSP00000435307.1:p.Arg151_Ser152insTer",
           ],

     30 => ["NC_000013.10:g.51519667_51519668insG",
           "ENST00000336617.2:c.616+1dup",
           "NC_000013.10:g.51519669dup"
           ],

     31 => ["NC_000006.11:g.30558478dup",
            "NC_000006.11:g.30558478dupA",
            "ENST00000396515.3:c.717dup",
            "ENST00000396515.3:c.717dupA"
           ],

     32 => ["NC_000001.10:g.154140412_154140414del", 
            "ENST00000368530.2:c.857_*1del",
            "NC_000001.10:g.154140414_154140416del"
           ],

     33 => ["NC_000012.11:g.102061070_102061071insT",
            "NC_000012.11:g.102061071dupT",
            "ENST00000360610.2:c.2336-439dup",
            "ENST00000360610.2:c.2336-439dupT",
            "NC_000012.11:g.102061071dup"
           ],

     34 => ["NC_000019.9:g.48836480_48836482del", 
            "NC_000019.9:g.48836480_48836482delGAG",
            "ENST00000293261.2:c.1376_1378del",
            "ENST00000293261.2:c.1376_1378delCCT",
            "ENSP00000293261.2:p.Ser459del",
            "NC_000019.9:g.48836480_48836482del"
           ]
);

## results which change on left-shifting - not shifted
my %test_output_no_shift = ( 
     1 => ["X:g.131215393dup",
           "A",
           "ENST00000298542.4:c.905+997dup", 
           "T",
           "",
           "duplication, intronic rc transcript"
          ],    
     2 => ["NC_000011.9:g.32417913_32417914insCCTACGAGTACTACC",
           "CCTACGAGTACTACC", 
           "ENST00000530998.1:c.451_452insGGTAGTACTCGTAGG",
           "GGTAGTACTCGTAGG", 
           "ENSP00000435307.1:p.Arg151_Ser152insTer",
            "insertion, stop gained [-1]"
           ],
     3 =>  ["NC_000013.10:g.51519667_51519668insG",    ##rs17857128
           "G",
           "ENST00000336617.2:c.615_616insG",
           "G",
           "ENSP00000337623.2:p.Glu206GlyfsTer13",
           "insertion, frameshift",
           ],
     4 =>  ["NC_000006.11:g.30558477_30558478insA",
           "A",
           "ENST00000396515.3:c.716_717insA",
           "A", 
           "ENSP00000379772.3:p.Ter239=",
           "insertion, stop retained"
          ],
     5 => ["NC_000001.10:g.154140413_154140415del",
           "-", 
           "ENST00000368530.2:c.856_858del",
           "-",
           "ENSP00000357516.2:p.Ter286delextTer56",
            "deletion, stop loss",
            "TAA"
          ],
     6 => ["NC_000012.11:g.102061070_102061071insT",
           "T",
           "ENST00000360610.2:c.2336-440_2336-439insT", 
           "T",
           "",
           "insertion, coding intron downstream"
          ],
     7 => ["NC_000019.9:g.48836478_48836480del",  ## rs149734771
           "-",
           "ENST00000293261.2:c.1376_1378del",
           "-",
           "ENSP00000293261.2:p.Ser459del",
           "deletion, inframe codon loss",
           "CCT"
          ],
    );



my %test_input_no_shift = ( 
     1 => ["X:g.131215393dup",
           "ENST00000298542.4:c.905+997dup", 
          ],    
     2 => ["NC_000011.9:g.32417913_32417914insCCTACGAGTACTACC",
           "ENST00000530998.1:c.451_452insGGTAGTACTCGTAGG",
           "ENSP00000435307.1:p.Arg151_Ser152insTer",
          ],
     3 => ["NC_000013.10:g.51519667_51519668insG", 
           "ENST00000336617.2:c.615_616insG",
           "ENSP00000337623.2:p.Glu206GlyfsTer13",
           ],
     4 => ["NC_000006.11:g.30558477_30558478insA",
           "ENST00000396515.3:c.716_717insA",
           "ENST00000396515.3:c.716_717insA(p.=)",
           "ENSP00000379772.3:p.Ter239=",
          ],
     5 => ["NC_000001.10:g.154140413_154140415del",
           "ENST00000368530.2:c.856_858del",
           "ENSP00000357516.2:p.Ter286delextTer56",
          ],
     6 => ["NC_000012.11:g.102061070_102061071insT",
           "ENST00000360610.2:c.2336-440_2336-439insT", 
          ],
     7 => ["NC_000019.9:g.48836478_48836480del",
           "NC_000019.9:g.48836478_48836480delAGG",
           "ENST00000293261.2:c.1376_1378delCCT",
           "ENSP00000293261.2:p.Ser459del",
           "ENST00000293261.2:c.1376_1378del"
          ]
);


## default 3'shifted mode
foreach my $num (keys %test_input){
  foreach my $desc ( @{$test_input{$num}}){
    get_results($desc, \@{$test_output{$num}}, 1);
  }
}
foreach my $num (keys %test_input_shifted){
  foreach my $desc ( @{$test_input_shifted{$num}}){
    get_results($desc, \@{$test_output_shifted{$num}}, 1);
  }
}

print "\n\nTesting no  shift\n\n";
## non- 3' shifted mode
foreach my $num (keys %test_input){
  foreach my $desc ( @{$test_input{$num}}){
    get_results($desc, \@{$test_output{$num}}, 1);
  }
}
foreach my $num (keys %test_input_no_shift){
  foreach my $desc ( @{$test_input_no_shift{$num}}){
    get_results($desc, \@{$test_output_no_shift{$num}}, 0);
  }
}

# test an input looking up ref
is(
  $vdba->get_variationFeatureAdaptor->fetch_by_hgvs_notation(-hgvs => 'NC_000002.11:g.46746465N>A', -replace_ref => 1)->allele_string,
  'G/A',
  'replace_ref'
);


sub get_results{

  my $input    = shift;
  my $output   = shift;
  my $shift_it = shift;

  my $variationfeature_adaptor    = $vdba->get_variationFeatureAdaptor;
  my $transcript_adaptor          = $cdba->get_transcriptAdaptor;
  my $transcript_variation_adaptor= $vdba->get_transcriptVariationAdaptor;

  ## set flag if shifting not required
  $transcript_variation_adaptor->db->shift_hgvs_variants_3prime( $shift_it ) ;
 
 
  if($DEBUG==1){     print "\n\n\nStarting $input is_shifted: $shift_it\n";}
 
  ## create variation feature from hgvs string
  my $variation_feature ;
        
  eval{
    $variation_feature = $variationfeature_adaptor->fetch_by_hgvs_notation( $input );
  };
  unless($@ eq ""){
    print "fetch by hgvs error : $@\n" if $DEBUG ==1;
    return;
  }


  my $transcript ;
  if( $output->[2] =~/\w+/){
    my $name    = (split/\./, $output->[2])[0];
    $transcript = $transcript_adaptor->fetch_by_stable_id( $name );
  }

  test_output($variation_feature, $input, $output, $transcript );

}

sub test_output{

  my $variation_feature  = shift;
  my $input              = shift;
  my $output             = shift;
  my $transcript         = shift;

  ## may need to flip to transcript strand

  my $allele = $output->[3]; 
  $allele = $output->[1] if  $input =~ /g\./ ; 

  ## genomic level
  my $hgvs_genomic      = $variation_feature->get_all_hgvs_notations("", "g");

  ok(  $hgvs_genomic->{$allele} eq $output->[0], "$input genomic level $output->[0], $output->[5]" );
  if($DEBUG==1){print "TEMP: $input => gen; expected $output->[0]\t returns: $hgvs_genomic->{$allele} for allele:$allele\n";} 


  ## transcript level - transcript to be supplied as ref feature  - alt allele may be complimented wrt genomic reference 
  return unless defined $transcript;

  my $hgvs_coding  = $variation_feature->get_all_hgvs_notations($transcript, "c");

  ok( $hgvs_coding->{$allele} eq $output->[2], "$input ->  transcript level from VF $output->[2], $output->[5]");
  if($DEBUG==1){  print "TEMP: $input => trans; expected $output->[2]\t returns: $hgvs_coding->{$allele} for allele:$allele\n\n";}

  ## protein level - transcript to be supplied as ref feature - alt allele may be complimented
  return unless $output->[4]; ## only look for protein level annotation if expected

  my $transcript_variations = $variation_feature->get_all_TranscriptVariations();
  foreach my $transcript_variation (@{$transcript_variations}){

    next unless $transcript_variation->transcript->stable_id() eq $transcript->stable_id() ;

    my $hgvs_protein  = $variation_feature->get_all_hgvs_notations($transcript,  "p",  $transcript->stable_id,  "",  $transcript_variation  );

    ##get from TVA too
    my $tva = $transcript_variation->get_all_alternate_BaseVariationFeatureOverlapAlleles();
    ok( $tva->[0]->hgvs_transcript() eq $output->[2], "$input ->  transcript level from TVA $output->[2], $output->[5]");

    if (defined $output->[6]){ ## check reference sequence as used in HGVS is as expected
      ok( $tva->[0]->hgvs_transcript_reference() eq $output->[6], "$input ->  transcript level correct ref,$output->[6] ");
    }


    ok( $hgvs_protein->{$allele} eq $output->[4], "$input -> protein level - $output->[4]  $output->[5]");
    if($DEBUG==1){   print "TEMP: $input => prot; expected $output->[4]\t returns: $hgvs_protein->{$allele} for allele:$allele\n";}
  }
}

done_testing();


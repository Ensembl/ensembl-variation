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

## HGVS no longer requires the reference allele to be supplied for del, indel, dup, inv
## We still support the parsing of the old format.

use strict;
use warnings;

use Test::More;
use Test::Exception;
use FindBin qw($Bin);

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_hgvs_alleles);

my %alleles = ("ENST00000306448.4:c.*14G>A"                => ["G", "A"],
               "ENST00000524249.1:n.775+16445delGinsTTGT"  => ["G", "TTGT"],
               "ENST00000368530.2:c.857_*1delAAT"          => ["AAT", "-"], 
               "ENST00000360610.2:c.2336-439dupT"          => ["-", "T"],
               "ENST00000368530.2:c.857_*1invAAT"          => ["AAT", "ATT"],
               "ENST00000396515.3:c.716_717insA"           => ["-", "A"],
               "ENST00000396515.3:c.716_717[3]A"           => ["A", "AAA"],
               "ENST00000522587.1:c.-310+750[13]A"         => ["A", "AAAAAAAAAAAAA"],
               "ENST00000306448.4:c.*14G="                 => ["G", "G"]
);


foreach my $hgvs( keys %alleles){

  my ($ref, $alt) = get_hgvs_alleles($hgvs);

  ok($ref eq $alleles{$hgvs}->[0], "$hgvs => ref");
  ok($alt eq $alleles{$hgvs}->[1], "$hgvs => alt");
}

throws_ok {get_hgvs_alleles("ENST00000306448.4:c.*14G)A") ; } qr/is unknown or could not be correctly recognized/, 'Throw on bad hgvs';

done_testing();


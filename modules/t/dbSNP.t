#!/usr/bin/env perl

use strict;
use warnings;


use Test::More;

use Data::Dumper;

use FindBin qw($Bin);

use Bio::EnsEMBL::Variation::Utils::dbSNP qw(get_alleles_from_pattern);

my @test_patterns = ( ["(AGAC)25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/", 
                      "(AGAC)25/(AGAC)26/(AGAC)27/(AGAC)28/(AGAC)29/(AGAC)30/(AGAC)31/(AGAC)32/(AGAC)33/(AGAC)34/(AGAC)35/(AGAC)36/(AGAC)37/(AGAC)38/(AGAC)39/(AGAC)40/(AGAC)41/(AGAC)42/(AGAC)43/(AGAC)44"],
                      ["(AT)7/8", "(AT)7/(AT)8"],
		      ["T/C", "T/C"], 
		      ["-/A", "-/A"],
		      ["(FB66D06.X1)",  "(FB66D06.X1)"],
                      ["(TG)20/21/A/G", "(TG)20/(TG)21/A/G"],
                      ["(TTA)1/6/7/-",  "TTA/(TTA)6/(TTA)7/-"]
		      );


foreach my $t (@test_patterns){

  my $als = get_alleles_from_pattern($t->[0]);

  my $ens_string = join "/", @{$als}; 
  warn "\n$ens_string\n";

  ok( $t->[1] eq $ens_string , "[ testing: $t->[0] -> $t->[1] ]");


}

done_testing();

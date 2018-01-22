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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use strict;
use warnings;

my @files = @ARGV;

print qq{
<html>
  <title>VEP documentation</title>
  <body>
};

foreach my $file(@files) {
  open IN, $file or die("ERROR: Could not read from file $file\n");

  my $in_body = 0;

  while(<IN>) {
    if(/\<title\>/) {
      # s/title/h1"/g;
      # print;
      next;
    }

    elsif(/\<\/?body\>/) {
      $in_body = 1 - $in_body;
    }

    elsif($in_body) {
      print;
    }
  }

  close IN;

  print '<div style="clear:both"><h1 style="page-break-after: always">&nbsp</h1></div>';
}

print "</body></html>\n";

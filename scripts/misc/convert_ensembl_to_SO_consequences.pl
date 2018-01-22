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


use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);

if(@ARGV and ($ARGV[0] =~ /^(-)+h/)) {
    print
qq{
Use this script to convert any text file containing Ensembl consequence types to SO terms.

perl convert_ensembl_to_SO_consequences.pl input.txt > converted.txt
};
    exit(0);
}

my @cons = grep {$_->{display_term} =~ /\w+/} sort {$b->rank <=> $a->rank} values %OVERLAP_CONSEQUENCES;

while(<>) {
    foreach my $con(@cons) {
        my $ens = $con->{display_term};
        my $so  = $con->{SO_term};
        
        # special case ESSENTIAL_SPLICE_SITE
        if($ens eq 'ESSENTIAL_SPLICE_SITE') {
            $so = 'splice_region_variant';
        }
        
        s/(\s|,|\&|^|\"|\')$ens/$1$so/g;
    }
    print;
}

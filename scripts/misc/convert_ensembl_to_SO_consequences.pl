#!/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

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
        
        s/(\s|,|\&|^)$ens/$1$so/g;
    }
    print;
}
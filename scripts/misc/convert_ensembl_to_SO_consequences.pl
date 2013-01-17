#!/usr/bin/env perl

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
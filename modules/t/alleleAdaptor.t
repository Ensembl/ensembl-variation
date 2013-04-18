use strict;
use warnings;

use Test::More;
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Population;

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all("$Bin/test.ensembl.registry.72");

my $vdba = $reg->get_DBAdaptor('human', 'variation');
my $cdba = $reg->get_DBAdaptor('human', 'core');

my $name = 'rs144235347';

my $va = $vdba->get_VariationAdaptor();
my $variation = $va->fetch_by_name($name);
my $aa = $vdba->get_AlleleAdaptor();

my $alleles = $aa->fetch_all_by_Variation($variation);

is(scalar @$alleles, 24, 'alleles for variation');
my $hash;
foreach (@$alleles) {
    my $allele = $_->allele();
    my $population = $_->population();
    my $freq = $_->frequency();
    $freq ||= 'No frequency';
    if (defined $population) {
        $hash->{$population->name}->{$allele} = $freq;
        #print $population->name, ' ', $allele, ' ', $freq, "\n";
    }
}
is($hash->{'1000GENOMES:phase_1_LWK'}->{'T'}, 0.0103092783505155, 'allele freq for population');
is($hash->{'1000GENOMES:phase_1_GBR'}->{'C'}, 1, 'allele freq for population');

done_testing();


use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 16;
}


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::TranscriptVariation;
use Bio::EnsEMBL::Transcript;


# test constructor

my $tr = Bio::EnsEMBL::Transcript->new();

my $v = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs2421',
                                                -source => 'dbSNP');

my $vf = Bio::EnsEMBL::Variation::VariationFeature->new
  (-start => 100,
   -end   => 100,
   -strand => 1,
   -variation_name => 'rs2421',
   -map_weight => 1,
   -allele_string => 'A/T',
   -variation => $v);



my $pep_allele = 'K/N';

my $cdna_start = 1127;
my $cdna_end   = 1127;
my $tl_start   = 318;
my $tl_end     = 318;
my $type       = 'NON_SYNONYMOUS_CODING';


my $trvar = Bio::EnsEMBL::Variation::TranscriptVariation->new
  (-variation_feature => $vf,
   -transcript        => $tr,
   -pep_allele_string => $pep_allele,
   -cdna_start        => $cdna_start,
   -cdna_end          => $cdna_end,
   -translation_start => $tl_start,
   -translation_end   => $tl_end,
   -type              => $type);

ok($trvar->variation_feature() == $vf);
ok($trvar->transcript() == $tr);
ok($trvar->pep_allele_string() eq $pep_allele);
ok($trvar->cdna_start() == $cdna_start);
ok($trvar->cdna_end() == $cdna_end);
ok($trvar->translation_start() == $tl_start);
ok($trvar->translation_end() == $tl_end);
ok($trvar->type() eq $type);


# test getter/setters
my $tr_new = Bio::EnsEMBL::Transcript->new();
ok(test_getter_setter($trvar, 'transcript', $tr_new));


my $vf_new = Bio::EnsEMBL::Variation::VariationFeature->new();
ok(test_getter_setter($trvar, 'variation_feature', $vf_new));


ok(test_getter_setter($trvar, 'pep_allele_string', $pep_allele));

ok(test_getter_setter($trvar, 'cdna_start', 1));
ok(test_getter_setter($trvar, 'cdna_end', 12));

ok(test_getter_setter($trvar, 'translation_start', 4));
ok(test_getter_setter($trvar, 'translation_end', 10));

ok(test_getter_setter($trvar, 'type', 'INTRONIC'));



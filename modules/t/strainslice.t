use lib 't';
use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 23;
}

use lib '/nfs/acari/yuan/ensembl/src/ensembl-variation_normalized-alleles/modules/';
use lib '/nfs/acari/yuan/ensembl/src/ensembl_normalized-alleles/modules/';

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');
my $core = $multi->get_DBAdaptor('core');
$vdb->dnadb($core);

my $slice_adaptor = $core->get_SliceAdaptor();
my $slice = $slice_adaptor->fetch_by_region("chromosome", '20',30263615,30263640); #two mixed snps

my $sliceStrain1 = $slice->get_by_strain('FVB'); #a insertion '-'->"GG" and a SNP 'G'->'A';
my $sliceStrain2 = $slice->get_by_strain('SPRET/Ei'); #a insertion '-'->'G' and a deletion 'G'->'-';
my $sliceStrain3 = $slice->get_by_strain('BALB/cJ'); #a no change and a mixed 'G'->'AA';

if ($verbose) {
  print STDERR "strain 1: "," "x14,$sliceStrain1->seq," start: ",$sliceStrain1->start," end: ",$sliceStrain1->end,"\n";
  print STDERR "strain 2: "," "x14,$sliceStrain2->seq," start: ",$sliceStrain2->start," end: ",$sliceStrain2->end,"\n";
  print STDERR "strain 3: "," "x14,$sliceStrain3->seq," start: ",$sliceStrain3->start," end: ",$sliceStrain3->end,"\n";
}

ok($sliceStrain1->seq eq "CGTAGAGGCTTTACTGAAAATTTTTGTG");
ok($sliceStrain1->start == 30263615 and $sliceStrain1->end ==30263640);

my $diff_strain_1 = $sliceStrain1->get_all_differences_Slice();
my $diff_strain_2 = $sliceStrain2->get_all_differences_Slice();
my $diff_strain_3 = $sliceStrain3->get_all_differences_Slice();

my ($diff1, $diff2)= @{$diff_strain_1};
ok($diff1->start==7 and $diff1->end==6);
ok($diff1->allele_string eq "GG");
ok($diff2->start==17 and $diff2->end==17);
ok($diff2->allele_string eq "A");

my ($diff3, $diff4)= @{$diff_strain_2};
ok($diff3->start==7 and $diff3->end==6);
ok($diff3->allele_string eq "G");
ok($diff4->start==17 and $diff2->end==17);
ok($diff4->allele_string eq "-");

my $differences12 = $sliceStrain1->get_all_differences_StrainSlice($sliceStrain2);

my ($difference12, $difference12a) = @{$differences12};
ok($difference12->start == 7 and $difference12->end ==8);
ok($difference12->allele_string eq "GG");
ok($difference12a->start == 19 and $difference12a->end ==19);
ok($difference12a->allele_string eq "A");

my $differences21 = $sliceStrain2->get_all_differences_StrainSlice($sliceStrain1);

my ($difference21, $difference21a) = @{$differences21};
ok($difference21->start == 7 and $difference21->end ==7);
ok($difference21->allele_string eq "G");
ok($difference21a->start == 18 and $difference21a->end ==17);
ok($difference21a->allele_string eq "-");

my $subSlice = $slice->sub_Slice(5,10,1);
my $subSlice_strain = $sliceStrain1->sub_Slice(5,10,1);

ok($subSlice->start == 30263619 and $subSlice->end == 30263624);
ok($subSlice->seq eq "GACTTT");
ok($subSlice_strain->start == 30263619 and $subSlice_strain->end == 30263622);
ok($subSlice_strain->seq eq "GAGGCT");
ok($sliceStrain1->subseq(5,10,1) eq "GAGGCT");





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
use Test::More;
use Data::Dumper;
use Test::Exception;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');


$vdb->dnadb($cdb);

my $sa   = $cdb->get_SliceAdaptor();
my $svfa = $vdb->get_StructuralVariationFeatureAdaptor();
my $sva  = $vdb->get_StructuralVariationAdaptor();
my $sta  = $vdb->get_StudyAdaptor();
my $vsa  = $vdb->get_VariationSetAdaptor();
my $srca = $vdb->get_SourceAdaptor();

my $dbID = 4509635;
my $outer_start = 7823440;
my $inner_start = 7823450;
my $inner_end = 8819360;
my $outer_end = 8819373;
my $var_name = 'esv93078';
my $source_name = 'DGVa';
my $study_name  = 'estd59';
my $chr = '8';

ok($svfa && $svfa->isa('Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor'));


my $svf = $svfa->fetch_by_dbID($dbID);
my $source = $svf->source;
my $slice  = $svf->slice;

ok($svf->outer_start() == $outer_start, "svf_id -> outer_start");
ok($svf->start()       == $outer_start, "svf_id -> start");
ok($svf->inner_start() == $inner_start, "svf_id -> inner_start");
ok($svf->inner_end()   == $inner_end,   "svf_id -> inner_end");
ok($svf->end()         == $outer_end,   "svf_id -> end") ;
ok($svf->outer_end()   == $outer_end,   "svf_id -> outer_end");
ok($svf->strand()      == 1,            "svf_id -> strand");
ok($svf->variation_name() eq $var_name, "svf_id -> varname" );
ok($source->name() eq $source_name,     "svf_id -> source" );
ok($svf->study->name() eq $study_name,  "svf_id -> study" );
ok($slice->seq_region_name() eq $chr,   "svf_id -> slice name");



# test fetch_all_by_StructuralVariation
print "\n# Test - fetch_all_by_StructuralVariation\n";
my $sv = $sva->fetch_by_dbID(3506221);
my $svfs = $svfa->fetch_all_by_StructuralVariation($sv);
ok(@$svfs == 1,                         "sv -> vf count ");
throws_ok { $svfa->fetch_all_by_StructuralVariation('sv'); } qr/SupportingStructuralVariation arg expected/, 'Throw on wrong SV argument.';
my $sv2 = Bio::EnsEMBL::Variation::StructuralVariation->new(-variation_name => 'test');
throws_ok { $svfa->fetch_all_by_StructuralVariation($sv2); } qr/StructuralVariation arg must have defined dbID/, 'Throw on missing dbID.';

$svf = $svfs->[0];
$source = $svf->source;
$slice  = $svf->slice;

ok($svf->dbID() == $dbID,               "sv -> svf id");
ok($slice->seq_region_name() eq $chr,   "sv -> slice name ");
ok($svf->outer_start() == $outer_start, "sv -> outer_start");
ok($svf->start()       == $outer_start, "sv -> start");
ok($svf->inner_start() == $inner_start, "sv -> inner_start");
ok($svf->inner_end()   == $inner_end,   "sv -> inner_end");
ok($svf->end()         == $outer_end,   "sv -> end") ;
ok($svf->outer_end()   == $outer_end,   "sv -> outer_end");
ok($svf->strand()      == 1,            "sv -> strand");
ok($svf->variation_name() eq $var_name, "sv -> name");
ok($source->name() eq $source_name,     "sv -> source" );
ok($svf->study->name() eq $study_name,  "sv -> study" );


my $study = $sta->fetch_by_name('estd55');

## Slices ##
my @sv_names = ('esv2751608','esv2421345','esv2758415');
my $sv_somatic = 'esv2221103';

my $slice_test = $sa->fetch_by_region('chromosome','18');
my $slice_soma = $sa->fetch_by_region('chromosome','1');
my $slice_set  = $sa->fetch_by_region('chromosome','8');

# test fetch all by Slice
print "\n# Test - fetch_all_by_Slice\n";
my $svfs1 = $svfa->fetch_all_by_Slice($slice_test);
my $sv_count = 0;
foreach my $svf1 (@$svfs1) {
  my $name = $sv_names[$sv_count];
  $sv_count ++;
  ok($svf1->variation_name() eq $name, "svby slice - $sv_count");
}

# test fetch all by Slice constraint
print "\n# Test - fetch_all_by_Slice_constraint\n";
my $constraint_2 = "svf.variation_name='".$sv_names[0]."'";
my $svfs2 = $svfa->fetch_all_by_Slice_constraint($slice_test,$constraint_2);
ok($svfs2->[0]->variation_name() eq $sv_names[0], "sv by slice constraint");

# test fetch all by Slice size range
print "\n# Test - fetch_all_by_Slice_size_range\n";
my $size_min = 40_000;
my $size_max = 200_000;
my $svfs2a = $svfa->fetch_all_by_Slice_size_range($slice_test, $size_min, $size_max);
ok(scalar @$svfs2a == 2, "fetch all by Slice size range, size min and size max");
$svfs2a = $svfa->fetch_all_by_Slice_size_range($slice_test, $size_min);
ok(scalar @$svfs2a == 3, "fetch all by Slice size range, size min");
#4169389 4288923 120_000
#15976013 16025663 50_000
#22294124 22772309 480_000
throws_ok { $svfa->fetch_all_by_Slice_size_range('slice', $size_min, $size_max); } qr/Slice arg expected/, 'Throw on wrong slice argument.';

print "\n# Test - fetch_all_somatic_by_Slice_size_range\n";
$size_min = 2;
$size_max = 50;
my $svfs2b = $svfa->fetch_all_somatic_by_Slice_size_range($slice_soma, $size_min, $size_max);
ok(scalar @$svfs2b == 1, "fetch all somatic by Slice size range, size min and size max");
$svfs2b = $svfa->fetch_all_somatic_by_Slice_size_range($slice_soma, $size_min);
ok(scalar @$svfs2b == 1, "fetch all somatic by Slice size range, size min");
throws_ok { $svfa->fetch_all_somatic_by_Slice_size_range('slice', $size_min, $size_max); } qr/Slice arg expected/, 'Throw on wrong slice argument.';

# test fetch all somatic by Slice
print "\n# Test - fetch_all_somatic_by_Slice\n";
my $svfs3 = $svfa->fetch_all_somatic_by_Slice($slice_soma);
ok($svfs3->[0]->variation_name() eq $sv_somatic, "somatic sv by slice ");

# test fetch all somatic by Slice constraint
print "\n# Test - fetch_all_somatic_by_Slice_constraint\n";
my $constraint_4 = "svf.variation_name='$sv_somatic'";
my $svfs4 = $svfa->fetch_all_somatic_by_Slice($slice_soma);
ok($svfs4->[0]->variation_name() eq $sv_somatic, "somatic sv by slice constraint");

# test fetch all somatic by Slice Source
print "\n# Test - fetch_all_somatic_by_Slice_Source\n";
my $svfs4a = $svfa->fetch_all_somatic_by_Slice_Source($slice_soma, $source, 0);
ok($svfs4a->[0]->variation_name() eq $sv_somatic, "somatic sv by slice source");
throws_ok { $svfa->fetch_all_somatic_by_Slice_Source('slice', $source, 0); } qr/Slice arg expected/, 'Throw on wrong slice argument.';
throws_ok { $svfa->fetch_all_somatic_by_Slice_Source($slice_soma, 'source', 0); } qr/Source arg expected/, 'Throw on wrong source argument.';

# test fetch all by Slice SO term
print "\n# Test - fetch_all_by_Slice_SO_term\n";
my $SO_term_6 = 'inversion';
my $svfs6 = $svfa->fetch_all_by_Slice_SO_term($slice_test,$SO_term_6);
ok($svfs6->[0]->variation_name() eq $sv_names[1], "sv by slice and SO term");
throws_ok { $svfa->fetch_all_by_Slice_SO_term('slice', $SO_term_6); } qr/Slice arg expected/, 'Throw on wrong slice argument.';
warns_like { $svfa->fetch_all_by_Slice_SO_term($slice_test, ['SO_term']); } qr/The SO term/, 'Warn on SO term has not been found.';

# test fetch all cnv probe by Slice
print "\n# Test - fetch_all_cnv_probe_by_Slice\n";
my $svfs7 = $svfa->fetch_all_cnv_probe_by_Slice($slice_test);
ok($svfs7->[0]->variation_name() eq 'CN_674347', "cnv probe by slice");

# test fetch all by Slice Study
print "\n# Test - fetch_all_by_Slice_Study\n";
my $svfs8 = $svfa->fetch_all_by_Slice_Study($slice_test, $study);
ok($svfs8->[0]->variation_name() eq $sv_names[0], "sv by slice and study");
throws_ok { $svfa->fetch_all_by_Slice_Study('slice', $study); } qr/Slice arg expected/, 'Throw on wrong slice argument.';
throws_ok { $svfa->fetch_all_by_Slice_Study($slice_test, 'study'); } qr/Study arg expected/, 'Throw on wrong study argument.';

# test fetch all by Slice Source
print "\n# Test - fetch_all_by_Slice_Source\n";
my $svfs8a = $svfa->fetch_all_by_Slice_Source($slice_test, $source, 0);
ok(scalar @$svfs8a == 3, "svfs by slice and source");
throws_ok { $svfa->fetch_all_by_Slice_Source('slice', $source); } qr/Slice arg expected/, 'Throw on wrong slice argument.';
throws_ok { $svfa->fetch_all_by_Slice_Source($slice_test, 'source'); } qr/Source arg expected/, 'Throw on wrong source argument.';

# test fetch all by Slice VariationSet
print "\n# Test - fetch_all_by_Slice_VariationSet\n";
my $set = $vsa->fetch_by_name('1000 Genomes - High coverage - Trios');
my @sv_sets = ('esv89107','esv93078');
my $svfs9 = $svfa->fetch_all_by_Slice_VariationSet($slice_set, $set);
ok($svfs9->[0]->variation_name() eq $sv_sets[0], "sv by slice and variation set - 1");
ok($svfs9->[1]->variation_name() eq $sv_sets[1], "sv by slice and variation set - 2");
throws_ok { $svfa->fetch_all_by_Slice_VariationSet('slice', $set); } qr/Slice arg expected/, 'Throw on wrong slice argument.';
throws_ok { $svfa->fetch_all_by_Slice_VariationSet($slice_set, 'set'); } qr/VariationSet arg expected/, 'Throw on wrong variation set argument.';

## Other ##

# test list_dbIDs
print "\n# Test - list_dbIDs\n";
my $svfs10 = $svfa->list_dbIDs();
ok($svfs10->[0] == 1850296, "sv id by list of dbIDs");

# test fetch all by Study
print "\n# Test - fetch_all_by_Study\n";
my $svfs11 = $svfa->fetch_all_by_Study($study);
ok($svfs11->[0]->variation_name() eq $sv_names[0], "sv by study");
throws_ok { $svfa->fetch_all_by_Study('study'); } qr/Study arg expected/, 'Throw on wrong study argument.';
my $study11 = Bio::EnsEMBL::Variation::Study->new(-name => 'new_study');
warns_like { $svfa->fetch_all_by_Study($study11); } qr/Study does not have dbID/, 'Warn on missing dbID.';

# test fetch all by Source
print "\n# Test - fetch_all_by_Source\n";
my $svfs12 = $svfa->fetch_all_by_Source($source);
ok($svfs12->[0]->variation_name() eq 'esv93078', "sv by source");
throws_ok { $svfa->fetch_all_by_Source('source'); } qr/Source arg expected/, 'Throw on wrong source argument.';
my $source12 = Bio::EnsEMBL::Variation::Source->new(-name => 'new_source');
warns_like { $svfa->fetch_all_by_Source($source12); } qr/Source does not have dbID/, 'Warn on missing dbID.';

## store ##

# test store
print "\n# Test - store\n";
delete $svf->{$_} for qw(dbID seq_region_start variation_name);
my $new_seq_region_start = 1000;
my $new_var_name = 'test';
$svf->start($new_seq_region_start);
$svf->variation_name($new_var_name);

ok($svfa->store($svf), "store");

my $svfs_store = $svfa->fetch_all_by_Slice_constraint($slice_set, "svf.seq_region_start=$new_seq_region_start");
ok($svfs_store && $svfs_store->[0]->seq_region_start == $new_seq_region_start && $svfs_store->[0]->variation_name eq $new_var_name, "fetch stored");

done_testing();



use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 17;
}


use Bio::EnsEMBL::Test::TestUtils;


use Bio::EnsEMBL::Test::MultiTestDB;
use Data::Dumper;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);

my $vgfa = $vdb->get_VariationGroupFeatureAdaptor();
my $vga  = $vdb->get_VariationGroupAdaptor();

ok($vgfa && $vgfa->isa('Bio::EnsEMBL::Variation::DBSQL::VariationGroupFeatureAdaptor'));

my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_region('chromosome', '21');



my $vgfs = $vgfa->fetch_all_by_Slice($slice);
print_feats($vgfs);

ok(@$vgfs == 8);


my $vgf = $vgfa->fetch_by_dbID(1);

ok($vgf->start() == 28622831);
ok($vgf->end()   == 28622831);
ok($vgf->strand() == 1);
ok($vgf->variation_group()->name() eq 'PERLEGEN:B000077');
ok($vgf->display_id() eq 'PERLEGEN:B000077');
ok($vgf->variation_group_name() eq 'PERLEGEN:B000077');
ok($vgf->slice()->name() eq $slice->name());


# test fetch_all_by_VariationGroup

my $vg = $vga->fetch_by_dbID(78);
$vgfs = $vgfa->fetch_all_by_VariationGroup($vg);
ok(@$vgfs == 1);
$vgf = $vgfs->[0];

ok($vgf->dbID() == 1);
ok($vgf->slice->name() eq $slice->name());
ok($vgf->start == 28622831);
ok($vgf->end() == 28622831);
ok($vgf->strand() == 1);
ok($vgf->variation_group_name() eq 'PERLEGEN:B000077');
ok($vgf->slice()->name() eq $slice->name());



sub print_feats {
  my $feats = shift;
  return if(!$verbose);

  foreach my $f (@$feats) {
    print $f->seq_region_name(), ':', $f->start(), '-', $f->end(), ' ',
          $f->display_id(), ' ', "\n";

    foreach my $v (@{$f->variation_group()->get_all_Variations()}) {
      print "  ", $v->name(), "\n";
    }
  }
}


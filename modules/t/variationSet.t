use strict;
use warnings;


use Test::More;


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::VariationSet;
use Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor;


my $name          = 'Phenotype-associated variations';
my $description   = 'Variations that have been associated with a phenotype';
my $short_name    = 'ph_variants';


# test constructor
my $variation_set = Bio::EnsEMBL::Variation::VariationSet->new
      ( -dbID        => 12,
        -name        => $name,
        -description => $description,
        -short_name  => $short_name

);



ok($variation_set->name() eq $name, "name");
ok($variation_set->short_name() eq $short_name, "short_name");
ok($variation_set->description() eq $description, "description");



# test getter/setters


ok(test_getter_setter($variation_set, 'name', 'new name', "get/set new name"));
ok(test_getter_setter($variation_set, 'description', 'new description', "get/set new description"));


print "Expecting error message:\n";
# test constructor with a value which is too high
my $variation_set2= Bio::EnsEMBL::Variation::VariationSet->new
      ( -dbID        => 102,
        -name        => $name,
        -description => $description,
        -short_name  => $short_name
      );

done_testing();


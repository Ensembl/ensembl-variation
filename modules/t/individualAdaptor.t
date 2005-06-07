
use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 10;
}


use Bio::EnsEMBL::Test::TestUtils;


use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');

my $pa = $vdb->get_PopulationAdaptor();
my $ia = $vdb->get_IndividualAdaptor();

ok($ia && $ia->isa('Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor'));

# test fetch by dbID

my $individuals = $ia->fetch_individual_by_synonym(7);
my $ind = shift @{$individuals};
ok($ind->name() eq 'CEPH104.11');
ok($ind->dbID() == 781);
ok($ind->description() eq 'CEPH/VENEZUELAN PEDIGREE 104');
ok($ind->gender eq 'Female');


# test fetch all by name

($ind) = @{$ia->fetch_all_by_name('CL63')};
ok($ind->name() eq 'CL63');
ok($ind->dbID() == 2781);


# test fetch_all_by_Population
my $pop = $pa->fetch_by_name('TSC-CSHL:CEL_caucasian_CEPH');
my $inds = $ia->fetch_all_by_Population($pop);
ok(@$inds == 661);

# test fetch_all_by_parent_Individual

($ind) = @{$ia->fetch_all_by_name('CEPH104.01')}; # father
$inds = $ia->fetch_all_by_parent_Individual($ind);
ok(@$inds == 10);

($ind) = @{$ia->fetch_all_by_name('CEPH104.02')}; # mother
$inds = $ia->fetch_all_by_parent_Individual($ind);
ok(@$inds == 10);

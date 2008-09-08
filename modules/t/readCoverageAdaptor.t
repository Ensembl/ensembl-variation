
use strict;
use warnings;

BEGIN { $| = 1;
        use Test;
        plan tests => 1;
}


use Bio::EnsEMBL::Test::TestUtils;


use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 1;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);
my $rca = $vdb->get_ReadCoverageAdaptor();
my $sa  = $db->get_SliceAdaptor();
my $pa  = $vdb->get_PopulationAdaptor();
my $ia = $vdb->get_IndividualAdaptor();

ok($rca && $rca->isa('Bio::EnsEMBL::Variation::DBSQL::ReadCoverageAdaptor'));

# test get read coverage in a region for a certain population in in a certain level
my $slice = $sa->fetch_by_region('chromosome', '1', 1, 200_000);

use Data::Dumper;
#print STDERR Dumper($slice);

my $samples = ["P109"];
foreach my $sample_name (@{$samples}){
  my $sample = shift @{$ia->fetch_all_by_name($sample_name)};
  #print STDERR $sample;
  my $coverage = $rca->fetch_all_by_Slice_Sample_depth($slice,$sample,1);
  #print STDERR "\nsample_name is $sample_name\n";
  #print STDERR Dumper($coverage);

}
foreach my $rc (@{$rca->fetch_all_regions_covered($slice,$samples)}) {
  #print STDERR Dumper($rc);
  print STDERR "range is ", $rc->[0], '-', $rc->[1], "\n";
}



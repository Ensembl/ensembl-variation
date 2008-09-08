
use strict;
use warnings;

BEGIN { $| = 1;
        use Test;
        plan tests => 3;
}


use Bio::EnsEMBL::Test::TestUtils;


use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 1;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);

my $va  = $vdb->get_VariationAdaptor();
my $iga = $vdb->get_IndividualGenotypeAdaptor();

ok($va && $va->isa('Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor'));
ok($iga && $iga->isa('Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor'));

my $var = $va->fetch_by_dbID(1199);
ok($var->name() eq 'rs1205');

my $igs = $iga->fetch_all_by_Variation($var);

#use Data::Dumper;
#print STDERR Dumper($igs);

print_feats($igs);

sub print_feats {
  my $feats = shift;
  return if(!$verbose);

  foreach my $f (@$feats) {
    print STDERR $f->start(), '-', $f->end(), '-', $f->allele1, '-', $f->allele2, ' ', $f->strand(), "\n";
  }	
	
}	

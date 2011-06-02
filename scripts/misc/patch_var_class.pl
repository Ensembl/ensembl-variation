use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Pipeline::SetVariationClass;

my $species;
my $reg_file;
my $transcripts;
my $variations;
my $disamb = 0;
my $verbose;
my $help;

GetOptions(
    'variations=s'  => \$variations,
    'species=s'     => \$species,
    'registry=s'    => \$reg_file,
    'verbose|v'     => \$verbose,
    'help|h'        => \$help,
);

unless ($species && $reg_file && $variations) {
    warn "Missing required argument\n" unless $help;
    $help = 1;
}

if ($help) {
    print "Usage: $0 --species <species> --registry <file> --variations <file> --help --verbose\n";
    exit(0);
}

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_all($reg_file);

my $vdba = $reg->get_DBAdaptor($species, 'variation') or die "Failed to get variation DBA\n";

my $va  = $reg->get_adaptor($species, 'variation', 'variation');

open VARIATIONS, "<$variations" or die "Failed to open $variations";

my $updated = 0;

while (<VARIATIONS>) {
    next if /^\s*$/;
    chomp;

    my $v = $va->fetch_by_name($_);
    
    unless ($v) {
        warn "No variation found for name $_\n";
        next;
    }
    
    my $svc = Bio::EnsEMBL::Variation::Pipeline::SetVariationClass->new;

    $svc->param('species', $species);
    $svc->param('ensembl_registry', $reg_file);
    $svc->param('variation_id_start', $v->dbID);
    $svc->param('variation_id_stop', $v->dbID);

    $svc->run;
    $updated++;
}

close VARIATIONS;

print "Updated $updated variations\n" if $verbose;


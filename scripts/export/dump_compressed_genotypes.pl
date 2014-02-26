
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

my $config = {};

GetOptions(
  $config,
  'help|h',                  # displays help message
  
  'host|h=s',                  # DB options
  'user|u=s',
  'password|p=s',
  'port|P=i',
  'database|d=s',
  
  'registry|r=s',
  'species|s=s',
) or die "ERROR: Failed to parse command-line flags\n";

die("ERROR: No species defined (--species)") unless defined($config->{species});

$config->{reg} = 'Bio::EnsEMBL::Registry';

if(defined($config->{registry})) {
  $config->{reg}->load_all($config->{registry});
}

else {
  for(qw(host user port)) {
    die("ERROR: --$_ not defined") unless defined $config->{$_};
  }
  
  $config->{reg}->load_registry_from_db(
    -host       => $config->{host},
    -user       => $config->{user},
    -pass       => $config->{password},
    -port       => $config->{port},
  );
}

# get adaptors
my $ga = $config->{reg}->get_adaptor($config->{species}, "variation", "individualgenotype");
die("ERROR: Failed to get genotype adaptor") unless defined($ga);

# first get seq_regions
my $sth = $ga->dbc->prepare(qq{
  SELECT DISTINCT(seq_region_id)
  FROM variation_feature
});
$sth->execute();

my ($sr_id, @sr_ids);
$sth->bind_columns(\$sr_id);
push @sr_ids, $sr_id while $sth->fetch();
$sth->finish();

$sth = $ga->dbc->prepare(qq{
  SELECT vf.variation_name, vf.seq_region_id, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, g.*
  FROM compressed_genotype_var g, variation_feature vf
  WHERE vf.seq_region_id = ?
  AND vf.variation_id = g.variation_id
}, {'mysql_use_result' => 1});

foreach my $sr_id(@sr_ids) {
  $sth->execute($sr_id);
  
  my ($n, $i, $s, $e, $r, $v, $ss, $g);
  $sth->bind_columns(\$n, \$i, \$s, \$e, \$r, \$v, \$ss, \$g);
  
	my %done;
	while($sth->fetch) {
		my @genotypes = unpack("(ww)*", $g);
		
		while(@genotypes) {
			my $individual_id = shift @genotypes;
			my $gt_code = shift @genotypes;
      
      $ss = 0 unless defined $ss;
      next if $done{$individual_id}{$gt_code}{$ss};
      $done{$individual_id}{$gt_code}{$ss} = 1;
			
      print join("\t", ($n, $i, $s, $e, $r, $v, $ss, $individual_id, $gt_code));
      print "\n";
      
			#my $igtype  = Bio::EnsEMBL::Variation::IndividualGenotype->new_fast({
			#	_variation_id => $v,
			#	subsnp        => $ss,
			#	adaptor       => $ga,
			#});
		}
	}
  
  $sth->finish();
}

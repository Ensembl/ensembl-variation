#!/usr/bin/env perl
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



=head1 CONTACT

  Use this script to dump genotypes from compressed_genotype_var
  Re-import to compressed_genotype_region using
  ensembl-variation/scripts/import/compress_genotypes_by_region_from_file.pl

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut
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
) or die "ERROR: Failed to parse command-line flags\n";

foreach(qw(host database)) {
  die("ERROR: No $_ defined (--$_)") unless defined($config->{$_});
}

$config->{port} ||= 3306;
$config->{user} ||= 'ensro';

Bio::EnsEMBL::Registry->load_registry_from_db(
  -host       => $config->{host},
  -user       => $config->{user},
  -pass       => $config->{password},
  -port       => $config->{port},
);

# connect to DB
my $connection_string = sprintf(
  "DBI:mysql(RaiseError=>1):host=%s;port=%s;db=%s",
  $config->{host},
  $config->{port},
  $config->{database}
);
my $dbc = DBI->connect(
  $connection_string, $config->{user}, $config->{password}
);

# first get seq_regions
my $sth = $dbc->prepare(qq{
  SELECT DISTINCT(seq_region_id)
  FROM variation_feature
});
$sth->execute();

my ($sr_id, @sr_ids);
$sth->bind_columns(\$sr_id);
push @sr_ids, $sr_id while $sth->fetch();
$sth->finish();

$sth = $dbc->prepare(qq{
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
			
      print join("\t", ($n, $i, $s, $e, $r, $v, $ss, $individual_id, $gt_code));
      print "\n";
		}
	}
  
  $sth->finish();
}

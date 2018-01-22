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

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use strict;
use warnings;

use Getopt::Long;
use ImportUtils qw(load);
use Bio::EnsEMBL::Registry;
use FindBin qw( $Bin );

my ($TMP_DIR, $TMP_FILE, $species, $selected_seq_region, $registry_file);

GetOptions(
	'tmpdir=s'  => \$TMP_DIR,
	'tmpfile=s' => \$TMP_FILE,
	'species=s' => \$species,
	'registry_file=s' => \$registry_file,
);

warn("Make sure you have an updated ensembl.registry file!\n");

usage('-TMP_DIR argument is required') if(!$TMP_DIR);
usage('-TMP_FILE argument is required') if(!$TMP_FILE);
usage('-species argument is required') if(!$species);

$registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );

Bio::EnsEMBL::Registry->set_disconnect_when_inactive();

my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbVar = $vdba->dbc;

$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;

# create mapping table if not exists
$dbVar->do("CREATE TABLE IF NOT EXISTS subsnp_map (
  variation_id int(11) unsigned NOT NULL,
  subsnp_id int(11) unsigned DEFAULT NULL
)");
$dbVar->do("TRUNCATE subsnp_map;");

# first dump out the variation_id and subsnp_id from allele
my $sth = $dbVar->prepare("SELECT variation_id, subsnp_id FROM allele WHERE subsnp_id IS NOT NULL", {mysql_use_result => 1});

my $file = "$TMP_DIR\/$TMP_FILE";

open OUT, ">$file" or die "Could not write to temp file $TMP_FILE\/$TMP_DIR\n";

my ($variation_id, $subsnp_id);

$sth->execute;
$sth->bind_columns(\$variation_id, \$subsnp_id);

print OUT "$variation_id\t$subsnp_id\n" while $sth->fetch;

$sth->finish();

# get the size of the file we dumped
my @stats = stat(OUT);

close OUT;

my $mem = int ($stats[7] /(1024 * 1024 * 1024)) * 2;

my $cmd = 'bsub -o '.$TMP_DIR.'/sort'.$$.'.out -J sort'.$$; 

if($mem < 2) {
	$cmd .= ' -q normal ';
}
else {
	$cmd .= ' -q hugemem -M'.$mem.'000000 -R"select[mem>'.$mem.'000] rusage[mem='.$mem.'000]" ';
}

$cmd = "echo 'sort -u $file > $file\_sorted; mv $file\_sorted $file' | ".$cmd;

print "$cmd\n";

# now submit a hugemem farm job to unique sort it, wait for it to run
system($cmd);
system("bsub -K -w 'done(sort".$$.")' -J waiting sleep 1");

print "Loading data\n";

# now reimport
load($dbVar, qw(subsnp_map variation_id subsnp_id));

print "Creating index\n";

$dbVar->do("CREATE INDEX variation_idx ON subsnp_map(variation_id);");

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


use Getopt::Long;
use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);
use Bio::EnsEMBL::Registry;

my $config = {};

GetOptions(
    $config,
	'host|h=s',
	'user|u=s',
	'password|p=s',
	'port|P=i',
	'version|v=i',
) or die "ERROR: Could not parse command line options\n";

$config->{host} ||= 'ens-staging';
$config->{user} ||= 'ensro';
$config->{port} ||= 3306;

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => $config->{host}, -user => $config->{user}, -port => $config->{port}, -db_version => $config->{version});

my $oa = $reg->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );

die "ERROR: Could not get ontology term adaptor\n" unless defined($oa);

my %term_list = ();

foreach my $con(values %OVERLAP_CONSEQUENCES) {
	my $obj = $oa->fetch_by_accession($con->SO_accession);
  
    die("ERROR: Failed to fetch DB object for ".$con->SO_term." (".$con->SO_accession.")\n") unless defined($obj);
	
	get_parents($obj, \%term_list);
}

print "$_\n" for sort keys %term_list;

sub get_parents {
	my $obj = shift;
	my $term_list = shift;
	
	get_parents($_, $term_list) for @{$obj->parents};
	
	$term_list->{$obj->name} = 1;
}
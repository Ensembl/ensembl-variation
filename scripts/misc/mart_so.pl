#!/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

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
	
	get_parents($obj, \%term_list);
}

print "$_\n" for sort keys %term_list;

sub get_parents {
	my $obj = shift;
	my $term_list = shift;
	
	get_parents($_, $term_list) for @{$obj->parents};
	
	$term_list->{$obj->name} = 1;
}
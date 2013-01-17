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


#ÊScript to dump out a table of variation sets that can be used in the documentation

use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my $species    = shift;
my $host       = shift;
my $db_version = shift;

die ("Species, db_host and db_version must be specified") unless ($species && $host && $db_version);

# Filters
my @filters = ('fail_');

# Load the registry from db
$registry->load_registry_from_db(
    -host => $host,
    -user => 'ensro',
		-db_version => $db_version
);

my $vdb = $registry->get_DBAdaptor($species,'variation');
my $dbVar = $vdb->dbc->db_handle;

my $table_header = qq{
	<tr>
		<th style="width:200px">Name</th>
		<th>Size</th>
		<th>Description</th>
	</tr>
};



my %pops = ('1000 Genomes' => { 'term' => '1000GENOMES:phase_1_%', 'constraint' => 'size is not null'},
            'HapMap'       => { 'term' => 'CSHL-HAPMAP:%'},
           );


foreach my $project (sort(keys(%pops))) {
  
	my $bg = '';
  
	my $term = $pops{$project}{term};
	my $constraint  = ($pops{$project}{constraint}) ? $pops{$project}{constraint}.' AND ' : '';
	
  print "<h4>Populations from the $project Project (Human)</h4>\n";
	$project =~ s/ /_/g;
	$project = lc $project;
  print "<table id=\"$project\" class=\"ss\">\n";
  print "$table_header\n";


  my $stmt = qq{ SELECT sample_id, name, size, description FROM sample WHERE $constraint name like ? };
  my $sth = $dbVar->prepare($stmt);
  $sth->execute($term);
	
  while(my @data = $sth->fetchrow_array) {
	  $data[2] = '-' if (!$data[2]);
		
		my @desc = split(/\.,/, $data[3]);
		$data[3] = "$desc[0]. $desc[1]." if scalar(@desc > 1);
		
		print "\t<tr$bg>\n";
	  print "\t\t<td>$data[1]</td>\n";
	  print "\t\t<td>$data[2]</td>\n";
	  print "\t\t<td>$data[3]</td>\n";
	  print "\t</tr>\n";
		
		if ($bg eq '') { $bg = ' class="bg2"'; }	
	  else { $bg = ''; }
	}
	
	print "</table>\n\n";
	$sth->finish;
}

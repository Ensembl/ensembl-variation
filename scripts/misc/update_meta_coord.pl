#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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


use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Getopt::Long;
use Bio::EnsEMBL::Utils::CliHelper;
use FileHandle;

my @table_names = qw(
  compressed_genotype_region
  phenotype_feature
  read_coverage
  variation_feature
  structural_variation_feature
);

my $script_opts = [{args => ['host', 'dbhost', 'h'], type => '=s'},
        {args => ['port', 'dbport', 'P'], type => ':i'},
        {args => ['user', 'dbuser', 'u'], type => '=s'},
        {args => ['pass', 'dbpass', 'p'], type => ':s'},
        {args => ['dbname',  'D'],         type => ':s'},
        {args => ['type', 't'],    type => ':s'},
        {args => ['version'],    type => ':i'},];

if ( scalar(@ARGV) == 0 ) {
  usage();
  exit 0;
}

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();
my $optsd = get_opts();
my $opts = $cli_helper->process_args( $optsd, \&usage );

my ($variation_dbas, $core_dbas);
use Data::Dumper;
print Dumper($opts);
if (defined $opts->{dbname}) {
  $variation_dbas = get_dbas($opts);
  die "Variation database required" unless ($opts->{dbname} =~ m/variation/ || $opts->{type} eq "private");
  $opts->{dbname} =~ s/variation/core/;
  $core_dbas = get_dbas($opts);
} else {
  my $version = $opts->{version};
  die "Version argument is required" unless ($version);
  $opts->{pattern} = 'variation_'.$version.'_\d+$';
  $variation_dbas = get_dbas($opts);
  $opts->{pattern} = 'core_'.$version.'_\d+$';
  $core_dbas = get_dbas($opts);
}

foreach my $core_dbname (keys %$core_dbas) {
  my $variation_dbname = $core_dbname;
  $variation_dbname =~ s/core/variation/;
  my $vdba = $variation_dbas->{$variation_dbname};
  next unless ($vdba);
  my $file = $vdba->dbc()->dbname() . "_meta_coord.backup";
  my $fh = FileHandle->new($file, 'w');
  my $dbh = $vdba->dbc->db_handle;
  my $sth = $dbh->prepare(q/SELECT * FROM meta_coord/);
  $sth->execute;
  while ( my @row = $sth->fetchrow_array() ) {
    @row = map {defined() ? $_ : '\N'} @row;    
    print $fh join(' ', @row), "\n";
  } 
  $sth->finish;
  $fh->close;
  $vdba->dbc()->sql_helper()->execute_update(-SQL => "DELETE FROM meta_coord", -PARAMS => []);
  foreach my $table_name (@table_names) {
    print("Updating $table_name table entries... ");
# insertions where end < start cause an error. Because both start and end are unsigned, negative is not allowed as unsigned integer. If you cast the values to signed before doing the difference it works or add a 1 at the beginning
    my $sql = "INSERT INTO meta_coord "
        . "SELECT '$table_name', s.coord_system_id, "
        . "MAX( 1 + cast( t.seq_region_end as signed) - cast( t.seq_region_start as signed) ) "
        . "FROM $table_name t, $core_dbname.seq_region s, $core_dbname.coord_system c "
        . "WHERE t.seq_region_id = s.seq_region_id AND c.coord_system_id=s.coord_system_id AND c.species_id=?"
        . "GROUP BY s.coord_system_id";
    $vdba->dbc()->sql_helper()->execute_update(
      -SQL => $sql,
      -PARAMS => [ $vdba->species_id() ] );
      print("done\n");
  }
}

sub get_dbas { 
  my $opts = shift;
  my $dbas = {};
  for my $db_args ( @{ $cli_helper->get_dba_args_for_opts( $opts, 0 ) } ) {
    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(%$db_args);
    my $dbname = $dba->dbc()->dbname();
    $dbas->{$dbname} = $dba;
  }
  return $dbas;
}

sub get_opts {
  my @dba_opts = map {
    my $opt = join '|', @{$_->{args}};
    $opt . $_->{type};
  } @{$script_opts};
  return \@dba_opts;
}

sub usage {
  print <<USAGE_END;
USAGE:
  $0 --dbhost=XhostX [--dbport=3306] \\
  \t--dbuser=XuserX --dbpass=XXX \\
  $0 \$(XserverX details script_db) --dbname XdbnameX
  $0 \$(XserverX details script_db) --version Xrelease-versionX
  $0 --help
  --help        Displays this help text.
This script will dump the current meta_coord table to a backup file in
the current directory. Update meta_coord table of a single database by
providing the dbname. Or update all databases for a given release by
providing the version. Then it will update the meta_coord table for the
data in the following tables:
USAGE_END

  print( "\t", join( "\n\t", @table_names ), "\n" );

}

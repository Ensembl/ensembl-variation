#!/usr/bin/env perl
# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use DBI;
use Getopt::Long;

my $config = {};

GetOptions(
  $config,
  'host|h=s',
  'user|u=s',
  'port|P=i',
  'password|p=s',
  'db|d=s',
  'help=s'
) or die "ERROR: Could not parse command line options\n";

if($config->{help} || !@ARGV) {
  usage();
  exit(0);
}

$config->{port} ||= 3306;

for(qw(host user port password db)) {
  die "ERROR: --$_ not defined\n" unless defined($config->{$_});
}

my $connection_string = sprintf(
  "DBI:mysql:database=%s:host=%s;port=%s",
  $config->{db},
  $config->{host},
  $config->{port}
);

# connect to DB
my $dbc = DBI->connect(
  $connection_string,
  $config->{user},
  $config->{password},
  {'RaiseError' => 1}
);

print "Getting phenotype descriptions\n";

# get phenotype descriptions
my $sth = $dbc->prepare("SELECT phenotype_id, description FROM phenotype");
my ($id, $desc);
$sth->execute();
$sth->bind_columns(\$id, \$desc);

my %descs;
$descs{$id} = $desc while($sth->fetch());
$sth->finish();

# map a reverse of the hash
my %id_by_desc = map {$descs{$_} => $_} keys %descs;

# store a "normalised version"
my %normalised;
my %denormalised;

foreach my $string(keys %id_by_desc) {
  my $normal = lc($string);
  $normal =~ s/\W//g;
  
  $normalised{$string} = $normal;
  push @{$denormalised{$normal}}, $string;
}

# create a table to store mappings
$dbc->do(qq{
  CREATE TABLE IF NOT EXISTS `phenotype_map` (
    `old_phenotype_id` int(11) unsigned NOT NULL,
    `new_phenotype_id` int(11) unsigned NOT NULL,
    KEY `old_phenotype_id` (`old_phenotype_id`),
    KEY `new_phenotype_id` (`new_phenotype_id`)
  )
});
# $dbc->do(qq{TRUNCATE TABLE phenotype_map});

$sth = $dbc->prepare("INSERT INTO phenotype_map(old_phenotype_id, new_phenotype_id) VALUES (?, ?)");

my $count = 0;

print "Finding semantic matches\n";

foreach my $normal(keys %denormalised) {
  if(scalar @{$denormalised{$normal}} > 1) {
    
    # create sort structure
    my @sort;
    
    foreach my $string(@{$denormalised{$normal}}) {
      my $c = $string =~ tr/[A-Z]/[A-Z]/;
      
      push @sort, {
        uc => $c,
        length => length($string),
        string => $string
      };
    }
    
    # sort
    @sort = map {$_->{string}} sort {$a->{uc} <=> $b->{uc} || $a->{length} <=> $b->{length}} @sort;
    my $main = shift @sort;
    
    my $new_id = $id_by_desc{$main};
    
    for(@sort) {
      $count++;
      $sth->execute($new_id, $id_by_desc{$_});
    }
  }
}

$sth->finish();

if($count) {
  my $a;
  
  print "\nFound $count semantic duplicate entries in phenotype\n";
  
  print "Backing up phenotype_feature entries to phenotype_feature_bak_map\n";

  # create backups of the entries we're going to change
  $dbc->do(qq{
    CREATE TABLE IF NOT EXISTS `phenotype_feature_bak_map`
    LIKE phenotype_feature
  });
  $a = $dbc->do(qq{
    INSERT IGNORE INTO phenotype_feature_bak_map
    SELECT pf.*
    FROM phenotype_feature pf, phenotype_map pm
    WHERE pm.old_phenotype_id = pf.phenotype_id
  });
  print "Backed up $a entries\n";
  
  print "Backing up phenotype entries to phenotype_bak_map\n";

  $dbc->do(qq{
    CREATE TABLE IF NOT EXISTS `phenotype_bak_map`
    LIKE phenotype
  });
  $a = $dbc->do(qq{
    INSERT IGNORE INTO phenotype_bak_map
    SELECT p.*
    FROM phenotype p, phenotype_map pm
    WHERE pm.old_phenotype_id = p.phenotype_id
  });
  print "Backed up $a entries\n";
  
  print "Updating entries in phenotype_feature\n";

  # alter the phenotype_feature entries
  $a = $dbc->do(qq{
    UPDATE phenotype_feature pf, phenotype_map pm
    SET pf.phenotype_id = pm.new_phenotype_id
    WHERE pf.phenotype_id = pm.old_phenotype_id
  });
  print "Updated $a entries\n";
  
  print "Removing entries from phenotype\n";

  # delete from phenotype
  $a = $dbc->do(qq{
    DELETE FROM p
    USING phenotype p, phenotype_map pm
    WHERE p.phenotype_id = pm.old_phenotype_id
  });

  print "Removed $a entries\n";
}

else {
  print "\n Found no duplicates, removing phenotype_map table\n";
  $dbc->do(qq{DROP TABLE phenotype_map});
}

print "\nAll done!\n";



sub usage {
  print qq{#---------------------------#
# rationalise_phenotypes.pl #
#---------------------------#

This script looks for semantic duplicates in the phenotype table using the
description column. Any phenotypes that match in all alphanumeric characters
(spaces and punctuation not considered) will be merged into one entry.

Removed phenotypes are backed up to phenotype_bak_map
Changed phenotype_feature entries are backed up to phenotype_feature_bak_map

Usage:
perl filter_vep.pl -h [host] -u [user] -p [pass] -P [port] -d [database]

};
}
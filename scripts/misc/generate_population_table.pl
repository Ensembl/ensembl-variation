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

my $evidence_icon_prefix = '/i/val/evidence_';
my $evidence_icon_suffix = '.png';
my $evidence_doc_url  = '#evidence_status';

my $table_header = qq{
  <tr>
    <th style="width:200px">Name</th>
    <th>Size</th>
    <th>Description</th>
  </tr>
};

my %pops = ('1000 Genomes'                   => { 'order'      => 1,
                                                  'species'    => 'Human',
                                                  'term'       => '1000GENOMES:phase_1_%',
                                                  'constraint' => 'size is not null',
                                                  'url'        => 'http://www.1000genomes.org',
                                                  'evidence'   => '1000Genomes'
                                                },
            'HapMap'                         => { 'order'    => 2,
                                                  'species'  => 'Human',
                                                  'term'     => 'CSHL-HAPMAP:%',
                                                  'url'      => 'http://hapmap.ncbi.nlm.nih.gov/index.html.en',
                                                  'evidence' => 'HapMap'
                                                },
            'Exome Sequencing Project (ESP)' => { 'order'    => 3,
                                                  'species'  => 'Human',
                                                  'term'     => 'ESP6500:%',
                                                  'url'      => 'http://evs.gs.washington.edu/EVS/',
                                                  'evidence' => 'ESP'
                                                }
           );


foreach my $project (sort{ $pops{$a}{'order'} <=> $pops{$b}{'order'} } keys(%pops)) {
  
  my $bg = '';

  my $term = $pops{$project}{'term'};
  my $spe  = $pops{$project}{'species'};
  my $constraint  = ($pops{$project}{'constraint'}) ? $pops{$project}{'constraint'}.' AND ' : '';
  my $project_word = ($project =~ /project/i) ? '' : ' Project';
  my $url = $pops{$project}{'url'};
  my $evidence = $pops{$project}{'evidence'};
  
  print qq{<h4>Populations from the <a href="$url" target="_blank" style="text-decoration:none">$project$project_word</a> ($spe)</h4>\n};

  my $project_id = $project;
  $project_id =~ s/ /_/g;
  $project_id = lc $project_id;
  print qq{<table id="$project_id" class="ss" style="margin-bottom:4px">\n  $table_header\n};

  my $stmt = qq{ SELECT population_id, name, size, description FROM population WHERE $constraint name like ? ORDER BY name};
  my $sth = $dbVar->prepare($stmt);
  $sth->execute($term);

  while(my @data = $sth->fetchrow_array) {
    $data[2] = '-' if (!$data[2]);
    my ($pop_name,$pop_suffix) = split(":", $data[1]);
    $pop_name .=  ":<b>$pop_suffix</b>" if ($pop_suffix);

    $data[3] = '-' if (!$data[3]);
    my @desc = split(/\.,/, $data[3]);
    $data[3] = "$desc[0]. $desc[1]." if scalar(@desc > 1);

    my $size = ($data[2] && $data[2] ne '-' ) ? $data[2] : get_size($data[0]);

    print qq{  <tr$bg>\n    <td>$pop_name</td>\n    <td style="text-align:right">$size</td>\n    <td>$data[3]</td>\n  </tr>\n};

    if ($bg eq '') { $bg = ' class="bg2"'; }
    else { $bg = ''; }
  }

  print "</table>\n";
  $sth->finish;

  # Evidence status
  if ($pops{$project}{'evidence'}) {
    my $evidence = $pops{$project}{'evidence'};
    my $evidence_img = "$evidence_icon_prefix$evidence$evidence_icon_suffix";
    my $margin = ($pops{$project}{'order'} eq scalar(keys(%pops))) ? '20px' : '25px';
    print qq{
<p style="margin-bottom:$margin">
Variants which have been discovered in this project have the "evidence status" <a href="$evidence_doc_url"><b>$evidence</b></a>.
On the website this corresponds to the icon <a href="$evidence_doc_url"><img class="_ht" src="$evidence_img" title="$evidence" style="vertical-align:bottom"/></a>.
</p>
    };
  }
}



sub get_size {
  my $pop_id = shift;

  my $stmt = qq{ SELECT count(*) FROM individual_population WHERE population_id=?};
  my $sth = $dbVar->prepare($stmt);
  $sth->execute($pop_id);
  my $size = ($sth->fetchrow_array)[0];
  $sth->finish;

  return ($size == 0) ? '-' : $size;
}

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


# Script to generate tables to display the clinical significance tables

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

# Load the registry from db
$registry->load_registry_from_db(
    -host => $host,
    -user => 'ensro',
    -db_version => $db_version,
);

my $vdb = $registry->get_DBAdaptor($species,'variation');
my $dbVar = $vdb->dbc->db_handle;


my %info = (
             'label'     => ['ClinVar','DGVa'],
             'clin_sign' => [
                             qq{ SELECT DISTINCT clinical_significance FROM variation WHERE variation_id NOT IN (SELECT variation_id FROM failed_variation) AND clinical_significance is not NULL},
                             qq{ SELECT DISTINCT clinical_significance FROM structural_variation 
                                 WHERE structural_variation_id NOT IN (SELECT structural_variation_id FROM failed_structural_variation)
                                 AND clinical_significance is not NULL
                               }
                            ],
             'query'     => [
                             qq{ SELECT name FROM variation
                                 WHERE FIND_IN_SET(?,clinical_significance)
                                   AND variation_id NOT IN (SELECT variation_id FROM failed_variation) 
                                 LIMIT 1
                               },
                             qq{ SELECT v1.variation_name FROM structural_variation v1, structural_variation v2, structural_variation_association vas
                                 WHERE v1.is_evidence=0
                                   AND FIND_IN_SET(?,v2.clinical_significance)
                                   AND v2.structural_variation_id=vas.supporting_structural_variation_id
                                   AND v2.is_evidence=1
                                   AND v1.structural_variation_id=vas.structural_variation_id
                                   AND v1.structural_variation_id NOT IN 
                                   (SELECT structural_variation_id FROM failed_structural_variation)
                                 LIMIT 1
                                }
                             ],
             'link'      => [ qq{/Homo_sapiens/Variation/Explore?v=},qq{/Homo_sapiens/StructuralVariation/Evidence?sv=}],
           );

my $html;
my $bg = '';
my $icon_path = '/i/val/clinsig_';
my $border_left = qq{ style="border-left:1px solid #BBB"};

# Clinical significance terms
my %clin_sign;
foreach my $type_stmt (@{$info{'clin_sign'}}) {
  my $sth = $dbVar->prepare($type_stmt);
  $sth->execute();
  while (my ($vals) = $sth->fetchrow_array()){
    foreach my $val (split(',',$vals)) {
      $clin_sign{$val} = 1;
    }
  }  
  $sth->finish;
}

# Clinical significance examples
my $html_content = add_table_header($info{'label'});

my $count = 0;
my $cs_term_count = scalar (keys %clin_sign);
foreach my $cs_term (sort(keys %clin_sign)) {
  $count ++;
  print STDERR qq{Term "$cs_term" done ($count/$cs_term_count)\n};
  my $icon_label = $cs_term;
     $icon_label =~ s/ /-/g;
  my $icon_col = qq{<td style="text-align:center"><img src="$icon_path$icon_label.png" title="$cs_term"/></td>};
  my $examples;
  for (my $i=0; $i < scalar(@{$info{'query'}});$i++) {
    $examples .= get_variant_example($i,$cs_term);
  }
  $html_content .= qq{  <tr$bg>$icon_col<td>$cs_term</td>$examples</tr>\n};
  $bg = set_bg();
}


## CONTENT ##
$html = qq{
<table class="ss" style="width:auto">
  $html_content
</table>
};
print $html;



#############
## METHODS ##
#############

sub set_bg {
  return ($bg eq '') ? ' class="bg2"' : '';
}


sub execute_stmt_one_result {
  my $stmt  = shift;
  my $value = shift;

  my $sth = $dbVar->prepare($stmt);
  $sth->execute($value);

  return $sth->fetchrow_array;
}


sub get_variant_example {
  my $order = shift;
  my $value = shift;
  
  my $var = (execute_stmt_one_result($info{'query'}->[$order],$value))[0];
  my $example = (defined($var)) ? sprintf (qq{<a href="%s%s">%s</a>},$info{'link'}->[$order],$var,$var) : '-';

  return qq{<td$border_left>$example</td>};
}

sub add_table_header {
  my $labels = shift;
  my $icon_column = qq{<th><span class="_ht ht" title="Icons designed by Ensembl">Icon</span></th>};

  my $eg_columns;
  foreach my $label (@$labels) {
    $eg_columns .= qq{<th$border_left>$label example</th>};
  }
  return qq{  <tr>$icon_column<th>Value</th>$eg_columns</tr>\n};
}


sub usage {
    
    print STDOUT qq{
Usage:

  $0 SPECIES DB_HOST DB_VERSION
  
Description:

  Prints html code for a table containing the different clinical significances for a species. The species 
  has to be specified on the command line as the first argument and the host database has to be 
  specified as the second argument
         
};
    
    exit(0);
}

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


my $list_stmt = qq{ SELECT value FROM attrib WHERE attrib_type_id IN 
                    (select attrib_type_id from attrib_type where code=?)
                  };
my $desc_stmt = qq{ SELECT name,description FROM attrib_type WHERE code=? };

my %types = (
  'dbsnp_clin_sig' => { 'query' => qq{ SELECT name FROM variation
                                     WHERE FIND_IN_SET(?,clinical_significance)
                                     AND variation_id NOT IN (SELECT variation_id FROM failed_variation) 
                                     LIMIT 1},
                        'link' => qq{/Homo_sapiens/Variation/Explore?v=},
                        'icon' => 1
                      },
  'dgva_clin_sig' => { 'query' => qq{ SELECT v1.variation_name FROM structural_variation v1, structural_variation v2, structural_variation_association vas, attrib a
                                    WHERE v2.structural_variation_id=vas.supporting_structural_variation_id
                                    AND v2.clinical_significance_attrib_id=a.attrib_id
																		AND a.value=?
                                    AND v1.structural_variation_id=vas.structural_variation_id
                                    AND v1.structural_variation_id NOT IN 
                                    (SELECT structural_variation_id FROM failed_structural_variation)
                                    LIMIT 1},
                       'link' => qq{/Homo_sapiens/StructuralVariation/Evidence?sv=}
                     },
);

my $html;
my $bg = '';
my $icon_path = '/i/val/clinsig_';


# Types
foreach my $type (keys %types) {
  my @list_val;
  
  my ($name,$desc) = execute_stmt_one_result($desc_stmt,$type);
  
  my $sth = $dbVar->prepare($list_stmt);
  $sth->execute($type);
  while (my $val = ($sth->fetchrow_array)[0]){
    push (@list_val,$val);
  }  
  $sth->finish;

  my $content = add_table_header($type);
  foreach my $value (sort(@list_val)) {
    my $example = get_variant_example($type,$value);
    my $label = ucfirst($value);
    my $icon_col = ($types{$type}{'icon'}) ? qq{<td style="text-align:center"><img src="$icon_path$value.png" title="$label"/></td>} : '';
    $content .= qq{    <tr$bg>$icon_col<td>$value$example</td></tr>\n};
    $bg = set_bg();
  }
  
  $html .= set_div($content, $type, $name, $desc);  
}
$html .= qq{\n<div style="clear:both" />};


## CONTENT ##
print $html;



sub set_div {
  my $content = shift;
  my $id   = shift;
  my $name = shift;
  my $desc = shift;
  
  $desc = (defined($desc)) ? qq{<p>$desc.</p>\n} : ''; 
  $name =~ /^(.+)(\s+clinical\s+significance)/;
  my $label = ($1 =~ /dbsnp/i) ? 'ClinVar' : $1;
  my $second_label = ($2) ? $2 : '';
  $name = qq{<span style="color:#333">$label</span>$second_label};
  $desc =~ s/dbSNP/ClinVar and dbSNP/;
  my $div = qq{
<div id="$id" style="float:left;margin-right:100px;">
  <span style="font-weight:bold">$name</span><br />
  $desc
  <table class="ss" style="width:auto">
   $content
  </table>
</div>};

  return $div;
}


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
  my $type  = shift;
  my $value = shift;

  return '' if (!defined($types{$type}));
  
  my $example = qq{</td><td>};
  
  my $var = (execute_stmt_one_result($types{$type}{query},$value))[0];
  $example .= (defined($var)) ? sprintf (qq{<a href="%s%s">%s</a>},$types{$type}{link},$var,$var) : '';
  
  return $example;
}

sub add_table_header {
  my $type = shift;
  my $icon_column = ($types{$type}{'icon'}) ? qq{<th><span class="_ht ht" title="Icons designed by Ensembl">Icon</span></th>} : '';
  my $ex_column = ($types{$type}) ? qq{<th>Example</th>} : '';
  return qq{ <tr>$icon_column<th>Value</th>$ex_column</tr>\n};
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

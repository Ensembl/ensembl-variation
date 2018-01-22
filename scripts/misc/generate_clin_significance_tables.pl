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


# Script to generate tables to display the clinical significance tables

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

my $registry = 'Bio::EnsEMBL::Registry';

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my ($species, $host, $port, $db_version, $output_file, $help);

GetOptions(
  'v=i'         => \$db_version,
  'o=s'         => \$output_file,
  'host=s'      => \$host,
  'port=i'      => \$port,
  'species|s=s' => \$species,
  'help!'       => \$help
);

usage ("Species, host, port, version and output_file must be specified") unless ($species && $host && $port && $db_version && $output_file);

# Load the registry from db
$registry->load_registry_from_db(
    -host => $host,
    -port => $port,
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


my %star_ranking = ( 'status' => { 'not classified by submitter'       => { 'term' => 'no assertion',             'stars' => 0 },
                                   'classified by single submitter'    => { 'term' => 'single submitter',         'stars' => 1 },
                                   'classified by multiple submitters' => { 'term' => 'multiple submitters',      'stars' => 2 },
                                   'reviewed by expert panel'          => { 'term' => 'reviewed by expert panel', 'stars' => 3 },
                                   'practice guideline'                => { 'term' => 'practice guideline',       'stars' => 4 }
                                 },
                     'query'  => [qq{ SELECT pf.object_id FROM phenotype_feature pf, phenotype_feature_attrib pfa, attrib_type a 
                                     WHERE pf.phenotype_feature_id=pfa.phenotype_feature_id 
                                       AND pfa.attrib_type_id=a.attrib_type_id
                                       AND a.code='review_status'
                                       AND pf.type='Variation'
                                       AND pfa.value LIKE ?
                                   }],
                     'link'   => [qq{/Homo_sapiens/Variation/Phenotype?v=}]
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
  my $icon_label = $cs_term;
     $icon_label =~ s/ /-/g;
  my $icon_col = qq{<td style="text-align:center"><img src="$icon_path$icon_label.png" title="$cs_term"/></td>};
  my $examples;
  for (my $i=0; $i < scalar(@{$info{'query'}});$i++) {
    $examples .= get_variant_example($i,$cs_term,\%info);
  }
  $html_content .= qq{  <tr$bg>$icon_col<td>$cs_term</td>$examples</tr>\n};
  $bg = set_bg();
  print STDERR qq{Term "$cs_term" done ($count/$cs_term_count)\n};
}

# Four-star rating
my $html_star_content;
foreach my $review_status (sort {$star_ranking{'status'}{$a}{'stars'} <=> $star_ranking{'status'}{$b}{'stars'}} keys(%{$star_ranking{'status'}})) {
  my $count_stars = $star_ranking{'status'}{$review_status}{'stars'};
  my $search_term = $star_ranking{'status'}{$review_status}{'term'};
  my $stars = qq{<span class="_ht" title="$review_status">};
  for (my $i=1; $i<5; $i++) {
    my $star_color = ($i <= $count_stars) ? 'gold' : 'grey';
    $stars .= qq{<img style="vertical-align:top" src="/i/val/$star_color\_star.png" alt="$star_color"/>};
  }
  $stars .= qq{</span>};
  my $star_example = get_variant_example(0,'%'.$search_term.'%',\%star_ranking);
  $html_star_content .= qq{  <tr$bg>\n    <td>$stars</td>\n    <td>$review_status</td>\n    $star_example\n  </tr>};
  $bg = set_bg();
}


## CONTENT ##
$html = qq{
<table class="ss" style="width:auto">
  $html_content
</table>
<p>Further explanations about the clinical significance terms are available on the <a href="http://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/">ClinVar website</a>.</p>
<h3>ClinVar rating</h3>
<p>We use the <a href="http://www.ncbi.nlm.nih.gov/clinvar/docs/details/#interpretation">ClinVar "four-star" rating</a> system to indicate the quality of classification/validation of the variant:</p>
<table class="ss" style="width:auto">
  <tr><th>Rating</th><th>Description</th><th$border_left>Example</th></tr>
$html_star_content
</table>
};


open  OUT, "> $output_file" or die $!;
print OUT $html;
close(OUT);




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
  my $data  = shift;
  
  my $var = (execute_stmt_one_result($data->{'query'}->[$order],$value))[0];
  my $example = (defined($var)) ? sprintf (qq{<a href="%s%s">%s</a>},$data->{'link'}->[$order],$var,$var) : '-';

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
  my $msg = shift;
  print qq{
  $msg
  Usage: perl generate_clin_significance_tables.pl [OPTION]
  
  Update the clinical significance tables in "data_description.html" (under public-plugins/ensembl/htdocs/info/genome/variation/).
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -o              An HTML output file name (Required)
    -host           Host of the human database (Required)
    -port           Human database port (Required)
    -species        Species name (Required)
  } . "\n";
  exit(0);
}

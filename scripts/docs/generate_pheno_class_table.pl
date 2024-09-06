#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
  <http://www.ensembl.org/Help/Contact>.

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

my %phenotype_class = ( 'class' => {'trait'       => {
                                      source =>'ClinVar', order=> 1},
                                   'tumour'       => {
                                      source => 'COSMIC', order=> 2},
                                   'non_specified' =>  {
                                      source => 'ClinVar', order=> 3} },
                     'query'  => qq{ SELECT pf.object_id, p.description FROM phenotype_feature pf, phenotype p, attrib a, attrib_type att, source s
                                     WHERE pf.phenotype_id=p.phenotype_id
                                       AND pf.source_id = s.source_id
                                       AND p.class_attrib_id=a.attrib_id
                                       AND a.attrib_type_id = att.attrib_type_id
                                       AND att.code='phenotype_type'
                                       AND pf.type='Variation'
                                       AND a.value LIKE ?
                                       AND s.name LIKE ?
                                       LIMIT 1
                                   },
                     'link'   => qq{/Homo_sapiens/Variation/Phenotype?v=},
                   );

my %source_url = (
  'ClinVar'    => 'https://www.ncbi.nlm.nih.gov/clinvar/',
  'COSMIC'     => 'https://cancer.sanger.ac.uk/cosmic',
);


my $html;
my $bg = '';
my $border_left = qq{ style="border-left:1px solid #BBB"};

# Phenotype classes variant examples
my $html_content = qq{<tr> <th>Classes</th> <th $border_left> Example phenotypes</th> <th $border_left>Example variants</th></tr>};

foreach my $pheno_class (sort {$phenotype_class{'class'}{$a}{order} <=> $phenotype_class{'class'}{$b}{'order'}} keys(%{$phenotype_class{'class'}})) {
  my $source = $phenotype_class{'class'}{$pheno_class}{source};

  my $class = qq{<span class="_ht" title="$pheno_class">};
  $class .= qq{</span>};

  my $class_example_data = get_variant_example('%'.$pheno_class.'%',$source,\%phenotype_class);
  my $class_example = $class_example_data->{html};
  my $pheno_desc = get_pheno_desc($pheno_class, $class_example_data, \%phenotype_class);

  $html_content .= qq{  <tr$bg>\n    <td><b>$pheno_class</b></td>\n    $pheno_desc\n    $class_example\n  </tr>};
  $bg = set_bg();
}


## CONTENT ##
$html = qq{
<table class="ss" style="width:auto">
  $html_content
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
  my $source = shift;

  my $sth = $dbVar->prepare($stmt);
  $sth->bind_param(1, $value);
  $sth->bind_param(2, $source);
  $sth->execute();

  return $sth->fetchrow_array;
}


sub get_variant_example {
  my $value = shift;
  my $source = shift;
  my $data  = shift;

  my @var = (execute_stmt_one_result($data->{'query'},$value, $source));

  my $example = (defined($var[0])) ? sprintf (qq{<a href="%s%s">%s</a>},$data->{'link'},$var[0],$var[0]) : '-';

  my %result = ('var' => \@var, 'html' => qq{<td$border_left>$example</td>} );
  return \%result;
}


sub get_pheno_desc {
  my $class = shift;
  my $class_example = shift;
  my $data = shift;

  my $pheno = $class_example->{var}[1];
  $pheno = 'not specified' if ($pheno eq "ClinVar: phenotype not specified");
  my $source = $data->{class}{$class}{source};
  my $desc_ex = (defined($class)) ? sprintf (qq{'$pheno' from <a href="%s">%s</a>},$source_url{$source}, $source) : '-';

  return qq{<td$border_left>$desc_ex</td>};
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

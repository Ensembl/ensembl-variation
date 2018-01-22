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


# Script to dump out a table of variation sets that can be used in the documentation

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::OntologyTermAdaptor;
use Bio::EnsEMBL::Variation::Utils::Constants qw(%VARIATION_CLASSES);
use Getopt::Long;

my $registry = my $onto_registy = 'Bio::EnsEMBL::Registry';

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my ($species, $host, $port, $ohost, $db_version, $output_file, $help);

GetOptions(
  'v=i'         => \$db_version,
  'o=s'         => \$output_file,
  'host=s'      => \$host,
  'port=i'      => \$port,
  'ohost=s'     => \$ohost,
  'species|s=s' => \$species,
  'help!'       => \$help
);

usage ("Species, host, port, version and output_file must be specified") unless ($species && $host && $port && $db_version && $output_file);
usage("Host providing the ontology database must be specified") unless ($ohost);

# Load the registry from db
$registry->load_registry_from_db(
    -host => $host,
    -port => $port,
    -user => 'ensro',
    -db_version => $db_version
);

my $internal_link = '/i/16/internal_link.png';

my %colour = (
    'copy_number_variation'          => '#8601AF',
    'insertion'                      => '#000099',
    'copy_number_gain'               => '#0000CC',
    'copy_number_loss'               => '#CC0000',
    'loss_of_heterozygosity'         => '#8601AF',
    'inversion'                      => '#C530FF',
    'complex_structural_alteration'  => '#007FFF',
    'complex_substitution'           => '#007FFF',
    'tandem_duplication'             => '#732E00',
    'short_tandem_repeat_variation'  => '#732E00',
    'mobile_element_insertion'       => '#000080',
    'novel_sequence_insertion'       => '#000080',
    'Alu_insertion'                  => '#000080',
    'translocation'                  => '#C3A4FF',
    'interchromosomal_translocation' => '#C3A4FF',
    'intrachromosomal_translocation' => '#C3A4FF',
    'deletion'                       => '#CC0000',
    'mobile_element_deletion'        => '#CC0000',
    'duplication'                    => '#0000CC',
    'indel'                          => 'gold',
    'probe'                          => '#4682B4',
);
my $default_colour = '#000000';


my %type = (
  '1' => { 'label' => 'Variant', 'type' => 'var'},
  '2' => { 'label' => '<span class="_ht ht" title="Structural variant">SV</span>', 'type' => 'sv'},
  '3' => { 'label' => 'Variant</li><li style="margin:4px 0px 0px"><span class="_ht ht" title="Structural variant">SV</span>',
           'type'  => 'both'},
  '4' => { 'label' => 'CNV probe', 'type' => 'sv'},
);


my %example_urls = (
  'var' => qq{<a href="/Homo_sapiens/Variation/Explore?v=####NAME####" target="_blank" title="See a variant example"><img src="$internal_link" alt="Link"/></a>},
  'sv'  => qq{<a href="/Homo_sapiens/StructuralVariation/Explore?sv=####NAME####" target="_blank" title="See a structural variant example"><img src="$internal_link" alt="Link"/></a>}
);

my $vdb = $registry->get_DBAdaptor($species,'variation');
my $dbVar = $vdb->dbc->db_handle;

my ($onto_host,$onto_port) = split(':',$ohost);
$onto_registy->load_registry_from_db(
    -host => $onto_host,
    -port => $onto_port,
    -user => 'ensro',
    -db_version => $db_version
);

my $odb = $onto_registy->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );

my %var_class;
my %sv_class;
my %both_class;
my %cnv_probe_class;

# Variation classes
my $stmt1  = qq{ SELECT distinct a.value FROM variation v, attrib a WHERE a.attrib_id=v.class_attrib_id};

my $stmt1a = qq{ SELECT v.name FROM variation v, attrib a WHERE v.variation_id NOT IN (SELECT variation_id FROM failed_variation) AND a.attrib_id=v.class_attrib_id AND a.value=? LIMIT 1 };

my $stmt1b = qq{ SELECT v.name FROM variation v, attrib a WHERE a.attrib_id=v.class_attrib_id AND a.value=? LIMIT 1 };

my $sth1  = $dbVar->prepare($stmt1);
my $sth1a = $dbVar->prepare($stmt1a);
my $sth1b = $dbVar->prepare($stmt1b);

$sth1->execute;
while(my $v_class = ($sth1->fetchrow_array)[0]) {
  $var_class{$v_class}{'class'} = $VARIATION_CLASSES{$v_class};
  
  $sth1a->execute($v_class);
  my $v_name = ($sth1a->fetchrow_array)[0];
  
  if (!$v_name) {
    $sth1b->execute($v_class);
    $v_name = ($sth1b->fetchrow_array)[0];
  }
  
  $var_class{$v_class}{'example'} = $v_name if ($v_name);
}
$sth1->finish; 
$sth1a->finish;


# Structural variation classes + CNV probes
my $stmt2  = qq{ SELECT distinct a.value FROM structural_variation v, attrib a WHERE a.attrib_id=v.class_attrib_id};

my $stmt2a = qq{ SELECT v.variation_name FROM structural_variation v, attrib a WHERE v.structural_variation_id NOT IN (SELECT structural_variation_id FROM failed_structural_variation) AND v.is_evidence=0 AND a.attrib_id=v.class_attrib_id AND a.value=? LIMIT 1 };

my $stmt2b = qq{ 
SELECT sv1.variation_name 
FROM 
  structural_variation sv1,
  structural_variation sv2,
  structural_variation_association sva,
  attrib a
WHERE 
  sv1.is_evidence=0 AND sv2.is_evidence=1 AND
  sv1.structural_variation_id = sva.structural_variation_id AND
  sv2.structural_variation_id = sva.supporting_structural_variation_id AND
  a.attrib_id=sv2.class_attrib_id AND a.value=? LIMIT 1
};

my $sth2  = $dbVar->prepare($stmt2);
my $sth2a = $dbVar->prepare($stmt2a);
my $sth2b = $dbVar->prepare($stmt2b);

$sth2->execute;
while(my $sv_class = ($sth2->fetchrow_array)[0]) {

  $sth2a->execute($sv_class);
  my $sv_name = ($sth2a->fetchrow_array)[0];

  if (!$sv_name) {
    $sth2b->execute($sv_class);
    $sv_name = ($sth2b->fetchrow_array)[0];
  }

  if ($sv_class =~ /probe/) {
    $cnv_probe_class{$sv_class}{'class'}   = $VARIATION_CLASSES{$sv_class};
    $cnv_probe_class{$sv_class}{'example'} = $sv_name if ($sv_name);
    next;
  }
  if ($var_class{$sv_class}) {
    $both_class{$sv_class}{'class'} = $VARIATION_CLASSES{$sv_class};
    $both_class{$sv_class}{'sv_example'} = $sv_name if ($sv_name);
    $both_class{$sv_class}{'v_example'}  = $var_class{$sv_class}{'example'} if ($var_class{$sv_class}{'example'});
    delete($var_class{$sv_class});
    next;
  }
  $sv_class{$sv_class}{'class'}   = $VARIATION_CLASSES{$sv_class};
  $sv_class{$sv_class}{'example'} = $sv_name if ($sv_name);
} 
$sth2->finish;
$sth2a->finish;
$sth2b->finish;



my $html = qq{
<table id="variation_classes" class="ss">
  <tr>
    <th style="width:8px;padding-left:0px;padding-right:0px;text-align:center">*</th>
    <th>SO term</th>
    <th>SO description</th>
    <th>SO accession</th>
    <th colspan="2">Called for (e.g.)</th>
  </tr>
};

my $bg = '';

# Variation
foreach my $var (sort(keys(%var_class))) {
  print_line($var,$var_class{$var},1);
}

# Structural variation
foreach my $sv (sort(keys(%sv_class))) {
  print_line($sv,$sv_class{$sv},2);
}

# Both Variation and Structural variation
foreach my $b (sort(keys(%both_class))) {
  print_line($b,$both_class{$b},3);
}

# CNV probe
foreach my $cnv (sort(keys(%cnv_probe_class))) {
  print_line($cnv,$cnv_probe_class{$cnv},4);
}

$html .= qq{</table>\n};

$html .= qq{
<p>
<b>*</b> Corresponding colours for the Ensembl web displays (only for Structural variants). 
The colours were originally based on the <a rel="external" href="http://www.ncbi.nlm.nih.gov/dbvar/content/overview/">dbVar</a> displays.
<p>
};

$html .= get_var_class_piechart();

open  OUT, "> $output_file" or die $!;
print OUT $html;
close(OUT);


sub print_line {
  my $so_term = shift;
  my $data    = shift;
  my $type_id = shift;
  
  my $so_acc   = $data->{'class'}{SO_accession};
  my $som_term = $data->{'class'}{somatic_display_term};
  my $t_name   = $type{$type_id}{'label'};
  my $examples = get_examples($type_id,$data);

  my $so_desc;
  print STDERR "Class $so_term ... ";
  my $oterm = $odb->fetch_by_accession($so_acc);
  if ($oterm) {
    my $desc = $oterm->definition;
    $desc =~ /"(.+)"/;
    $so_desc = $1;
  }
  print STDERR "done\n";
  
  my $class_col = '';
  if ($colour{$so_term}) {
    $class_col = $colour{$so_term};
    $class_col = qq{;background-color:$class_col};
  } elsif ($type_id == 3 || $type_id == 2) {
    $class_col = qq{;background-color:$default_colour};
  }
  my $border = ($class_col eq '') ? '' : ';border-top:1px solid #FFF';
  
  $html .= qq{
  <tr$bg>
    <td style="padding:0px;margin:0px;$class_col$border"></td>
    <td style="font-weight:bold">$so_term</td>
    <td>$so_desc</td>
    <td><a rel="external" href="http://www.sequenceontology.org/miso/current_release/term/$so_acc">$so_acc</a></td>
    <td>
      <ul style="margin:0px;padding-left:1em">
        <li style="margin:0px">$t_name</li>
      </ul>
    </td>
    <td style="padding-left:0px;width:16px">$examples</td>
  </tr>};
  
  
  if ($bg eq '') { $bg = ' class="bg2"'; }  
  else { $bg = ''; }
}


sub get_examples {
  my $vtype = shift;
  my $data  = shift;
  
  my $example_type = $type{$vtype}{'type'};
  my $html = "";
  if ($example_type eq 'both') {
  
    my $var_url = '';
    if ($data->{'v_example'}) {
      my $v_name = $data->{'v_example'};
      
      $var_url = $example_urls{'var'};
      $var_url =~ s/####NAME####/$v_name/;
    }
    $html .= qq{<div>$var_url</div>};
    
    my $sv_url = '';
    if ($data->{'sv_example'}) {
      my $sv_name = $data->{'sv_example'};
      
      $sv_url = $example_urls{'sv'};
      $sv_url =~ s/####NAME####/$sv_name/;
    }
    $html .= qq{<div style="padding-top:1px">$sv_url</div>};
  }
  elsif ($example_type eq 'var') {
    my $var_url = '';
    if ($data->{'example'}) {
      my $v_name = $data->{'example'};
      
      $var_url = $example_urls{'var'};
      $var_url =~ s/####NAME####/$v_name/;
    }
    $html .= qq{<div>$var_url</div>};
  }
  elsif ($example_type eq 'sv') {
    my $sv_url = '';
    if ($data->{'example'}) {
      my $sv_name = $data->{'example'};
      
      $sv_url = $example_urls{'sv'};
      $sv_url =~ s/####NAME####/$sv_name/;
    }
    $html .= qq{<div>$sv_url</div>};
  }
  
  return $html;
}

sub get_var_class_piechart {

  print STDERR "Generate the pie chart of the variation class distribution ...\n";
  
  my $species_name = ($species =~ /homo|human/i) ? 'Human' : $species;
  my $stmt = qq{ SELECT distinct a.value, count(v.variation_id) av FROM variation v,attrib a 
                 WHERE a.attrib_id=v.class_attrib_id GROUP by a.attrib_id ORDER BY av DESC
               };
  my $sth = $dbVar->prepare($stmt);
  $sth->execute;
  my @data;
  while(my ($class,$count) = $sth->fetchrow_array) {
    push @data, "[$count,'$class']";
  }
  $sth->finish;
  
  my $html .= sprintf( 
  qq{
<div class="js_panel" id="pie_chart_panel" style="margin-top:15px;margin-bottom:20px">
  <input class="panel_type" type="hidden" value="Piechart">
  <input class="graph_config" type="hidden" name="legendpos" value="'east'" />
  <input class="graph_dimensions" type="hidden" value="[90,80,80]" />
  
  <input type="hidden" class="graph_data" value="[%s]" />
  
  <div class="pie_chart_classes" title="classes" style="width:400px;height:220px;border:1px solid #CCC;border-radius:8px;box-shadow:0 1px 3px #666">
     <h3 style="padding:10px 10px 5px">%s variant class distribution - Ensembl %s</h4>
     <div id="graphHolder0"></div>
  </div>
</div>
  },
  join(',',@data),
  $species_name,
  $db_version
  );
  
  return $html;
} 

sub usage {
  my $msg = shift;
  print qq{
  $msg
  Usage: perl generate_classes_table.pl [OPTION]
  
  Update the classes table in "data_description.html" (under public-plugins/ensembl/htdocs/info/genome/variation/).
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -o              An HTML output file name (Required)
    -host           Host of the human database (Required)
    -port           Human database port (Required)
    -ohost          Host name where the ontology database is stored, with the port, e.g. ensembldb.ensembl.org:1234 (Required)
    -species        Species name (Required)
  } . "\n";
  exit(0);
}

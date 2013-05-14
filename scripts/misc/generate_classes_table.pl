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


# Script to dump out a table of variation sets that can be used in the documentation

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Constants qw(%VARIATION_CLASSES);

my $registry = 'Bio::EnsEMBL::Registry';

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

my %colour = (
    'copy_number_variation'         => '#000000',
    'insertion'                     => '#FFCC00',
    'copy_number_gain'              => '#0000CC', 
    'copy_number_loss'              => '#CC0000',
    'inversion'                     => '#9933FF', 
    'complex_structural_alteration' => '#99CCFF',
    'tandem_duplication'            => '#732E00',
    'mobile_element_insertion'      => '#FFCC00',
    'translocation'                 => '#C3A4FF',
    'deletion'                      => '#CC0000',
    'duplication'                   => '#000000',
    'probe'                         => '#4682B4',
);
my $default_colour = '#B2B2B2';


my %type = (
	'1' => 'Variation',
	'2' => 'Structural variation', 
	'3' => 'Variation<br />Structural variation',
	'4' => 'CNV probe',
);

my $vdb = $registry->get_DBAdaptor($species,'variation');
my $dbVar = $vdb->dbc->db_handle;

my %var_class;
my %sv_class;
my %both_class;
my %cnv_probe_class;

# Variation classes
my $stmt1 = qq{ SELECT distinct a.value FROM variation v, attrib a WHERE a.attrib_id=v.class_attrib_id};
my $sth1  = $dbVar->prepare($stmt1);
$sth1->execute;
while(my $v_class = ($sth1->fetchrow_array)[0]) {
	$var_class{$v_class} = $VARIATION_CLASSES{$v_class};
} 
$sth1->finish;


# Structural variation classes + CNV probes
my $stmt2 = qq{ SELECT distinct a.value FROM structural_variation v, attrib a WHERE a.attrib_id=v.class_attrib_id};
my $sth2  = $dbVar->prepare($stmt2);
$sth2->execute;
while(my $sv_class = ($sth2->fetchrow_array)[0]) {
	if ($sv_class =~ /probe/) {
		$cnv_probe_class{$sv_class} = $VARIATION_CLASSES{$sv_class};
		next;
	}
	if ($var_class{$sv_class}) {
		$both_class{$sv_class} = $VARIATION_CLASSES{$sv_class};
		delete($var_class{$sv_class});
		next;
	}
	$sv_class{$sv_class} = $VARIATION_CLASSES{$sv_class};
} 
$sth2->finish;




my $html = qq{
<table id="variation_classes" class="ss">
  <tr>
    <th style="width:8px;padding-left:0px;padding-right:0px;text-align:center">*</th>
    <th>SO term</th>
    <th>SO description</th>
    <th>SO accession</th>
    <th>Ensembl term</th>
    <th>Called for</th>
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
<b>*</b> Corresponding colours for the Ensembl web displays (only for Structural variations). 
The colours are based on the <a rel="external" href="http://www.ncbi.nlm.nih.gov/dbvar/content/overview/">dbVar</a> displays.
<p>
};

$html.= get_var_class_piechart();

print $html;


sub print_line {
	my $so_term = shift;
	my $data    = shift;
	my $type_id = shift;
	
	my $e_class  = $data->{display_term};
	my $so_acc   = $data->{SO_accession};
	my $som_term = $data->{somatic_display_term};
	my $t_name   = $type{$type_id};
	
	my $so_desc;
	`wget http://www.sequenceontology.org/browser/current_release/export/term_only/csv_text/$so_acc`;
	
	
	if (-e $so_acc) {
		my $content = `grep -w $so_acc $so_acc`;
		$so_desc = (split("\t",$content))[2];
		`rm -f ./$so_acc`;
	}
	
	my $class_col = '';
	if ($colour{$so_term}) {
		$class_col = $colour{$so_term};
		$class_col = qq{;background-color:$class_col};
	} elsif ($type_id == 3 || $type_id == 2) {
	  $class_col = qq{;background-color:$default_colour};
	}
	my $border = ($class_col eq '') ? '' : ';border-top:1px solid #FFF';
	
	my $rowspan = ($so_term eq 'probe') ? '' : ' rowspan="2"';
	
	$html .= qq{
	<tr$bg>
		<td$rowspan style="padding:0px;margin:0px$class_col$border"></td>
		<td$rowspan>$so_term</td>
		<td$rowspan>$so_desc</td>
		<td$rowspan><a rel="external" href="http://www.sequenceontology.org/miso/current_release/term/$so_acc">$so_acc</a></td>
		<td>$e_class</td>
		<td$rowspan>$t_name</td>
	</tr>};
	
	if ($so_term ne 'probe') {
	  $html .= qq{
	<tr$bg>
		<td>$som_term</td>
	</tr>\n};
  }
	
	if ($bg eq '') { $bg = ' class="bg2"'; }	
	else { $bg = ''; }
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
  
	<div class="pie_chart_classes" title="classes" style="width:400px;height:220px;border:1px solid #AAA;border-radius:10px;-moz-border-radius:10px;-webkit-border-radius:10px;-webkit-box-shadow:0 1px 3px #666;">
   	<h3 style="padding:10px 10px 5px">%s Variation class distribution - Ensembl %s</h4>
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


#!/usr/bin/env perl

# Script to dump out a table of variation sets that can be used in the documentation

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Constants qw(%VARIATION_CLASSES);

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


my %colour = (
	'copy_number_variation'         => '#000000',
	'insertion'                     => '#FFCC00',
	'copy_number_gain'              => '#0000FF', 
	'copy_number_loss'              => '#FF0000',
	'inversion'                     => '#9933FF', 
	'complex_structural_alteration' => '#99CCFF',
	'tandem_duplication'            => '#732E00',
	'mobile_element_insertion'      => '#FFCC00',
	'translocation'                 => '#C3A4FF',
);

my %type = (
	'1' => 'Variation',
	'2' => 'Structural variation', 
	'3' => 'Variation<br />Structural variation'
);

my $vdb = $registry->get_DBAdaptor($species,'variation');
my $dbVar = $vdb->dbc->db_handle;

my %var_class;
my %sv_class;
my %both_class;

# Variation classes
my $stmt1 = qq{ SELECT distinct a.value FROM variation v, attrib a WHERE a.attrib_id=v.class_attrib_id};
my $sth1  = $dbVar->prepare($stmt1);
$sth1->execute;
while(my $v_class = ($sth1->fetchrow_array)[0]) {
	$var_class{$v_class} = $VARIATION_CLASSES{$v_class};
} 
$sth1->finish;


# Structural variation classes
my $stmt2 = qq{ SELECT distinct a.value FROM structural_variation v, attrib a WHERE a.attrib_id=v.class_attrib_id};
my $sth2  = $dbVar->prepare($stmt2);
$sth2->execute;
while(my $sv_class = ($sth2->fetchrow_array)[0]) {
	next if ($sv_class =~ /probe/);
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

$html .= qq{ </table>\n};

$html .= qq{
<p>
* Corresponding colours for the Ensembl web displays (only for Structural variations). 
The colours are based on the <a rel="external" href="http://www.ncbi.nlm.nih.gov/dbvar/content/overview/">dbVar</a> displays.
<p>
};

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
	`wget http://www.sequenceontology.org/browser/current_cvs/export/term_only/csv_text/$so_acc`;
	
	
	if (-e $so_acc) {
		my $content = `grep -w $so_acc $so_acc`;
		$so_desc = (split("\t",$content))[2];
		`rm -f ./$so_acc`;
	}
	
	my $class_col = '';
	if ($colour{$so_term}) {
		$class_col = $colour{$so_term};
		$class_col = qq {;background-color:$class_col};
	}
	
	$html .= qq{
		<tr$bg>
			<td rowspan="2" style="padding:0px;margin:0px$class_col"></td>
			<td rowspan="2">$so_term</td>
			<td rowspan="2">$so_desc</td>
			<td rowspan="2"><a rel="external" href="http://www.sequenceontology.org/miso/current_release/term/$so_acc">$so_acc</a></td>
			<td>$e_class</td>
			<td rowspan="2">$t_name</td>
		</tr>
		<tr$bg>
			<td>$som_term</td>
		</tr>
	};

	if ($bg eq '') { $bg = ' class="bg2"'; }	
	else { $bg = ''; }
}

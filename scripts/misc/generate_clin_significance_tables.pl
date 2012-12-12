#! perl -w

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
  'dbsnp_clin_sig' => { query => qq{ SELECT v.name FROM variation v, attrib a 
                                      WHERE v.clinical_significance_attrib_id=a.attrib_id
                                      AND a.value=? 
																			AND v.variation_id NOT IN (SELECT variation_id FROM failed_variation) 
																			LIMIT 1},
											  link => qq{/Homo_sapiens/Variation/Explore?v=}
											},
  'dgva_clin_sig' => { query => qq{ SELECT v.variation_name FROM structural_variation v, structural_variation_association vas, structural_variation_annotation va, attrib a 
                                     WHERE va.structural_variation_id=vas.supporting_structural_variation_id
												             AND va.clinical_attrib_id=a.attrib_id
																		 AND v.structural_variation_id=vas.structural_variation_id
                                     AND a.value=? 
																		 AND v.structural_variation_id NOT IN 
																		 (SELECT structural_variation_id FROM failed_structural_variation)
																		 LIMIT 1},
											 link => qq{/Homo_sapiens/StructuralVariation/Evidence?sv=}
										 },
);

my $html;
my $bg = '';



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
		$content .= qq{    <tr$bg><td>$value$example</td></tr>\n};
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
	
	$desc = (defined($desc)) ? qq{  <p>$desc.</p>\n} : ''; 
	$name =~ /^(.+)(\s+clinical\s+significance)/;
	
	$name = qq{<span style="color:#333">$1</span>$2};
	
	my $div = qq{
<div id="$id" style="float:left;margin-right:100px;">
  <b>$name</b><br />
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
	my $ex_column = ($types{$type}) ? qq{<th>Example</th>} : '';
	return qq{    <tr><th>Value</th>$ex_column</tr>\n};
}

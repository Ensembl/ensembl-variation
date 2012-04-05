# Script to generate an HTML page containing the variation sources of each species

use Bio::EnsEMBL::Registry;
use DBI;
use strict;
use POSIX;
use Getopt::Long;

###############
### Options ###
###############
my ($e_version,$html_file,$source,$s_version,$s_description,$s_url,$s_type,$hlist,$phost,$help);
## EG options
my ($site, $etype);

usage() if (!scalar(@ARGV));
 
GetOptions(
	   'v=s'   => \$e_version,
	   'o=s'   => \$html_file,
	   'help!' => \$help,
	   'hlist=s' => \$hlist,
	   'phost=s' => \$phost,
	   'site=s' => \$site,
	   'etype=s' =>  \$etype
);

if (!$e_version) {
	print "> Error! Please give an Ensembl version, using the option '-v' \n";
	usage();
}
if (!$html_file) {
	print "> Error! Please give an output file using the option '-o'\n";
	usage();
}
if (!$phost) {
	print "> Error! Please give host name where the previous databases are stored using the option '-phost'\n";
	usage();
}
if (!$hlist) {
	print "> Error! Please give the list of host names where the new databases are stored using the option '-hlist'\n";
	usage();
}

usage() if ($help);

my $server_name = 'http://static.ensembl.org';
my $ecaption = 'Ensembl';
my $previous_host = $phost;
my @hostnames = split /,/, $hlist;

if ($site) {
  $server_name = $site;
}
if ($etype) {
  $ecaption .= ' '.ucfirst($etype);
}
# Settings
my $database = "";
my $login = "ensro";
my $pswd = "";
my $sep = "\t";
my $start = 0;
my %colors = ( 'version' => '#090', 'source'  => '#00F' );

##############
### Header ###
##############
my $html_header = qq{
<html>
<head>
	<meta http-equiv="CONTENT-TYPE" content="text/html; charset=utf-8" />
	<title>Variation Sources</title>
</head>

<body>
<b>
<a href="index.html">About Ensembl Variation</a> |
<a href="data_description.html">Data Description</a> | 
<a href="predicted_data.html">Predicted Data</a> |
<a href="database.html">Database Description</a> | 
<a href="/info/docs/api/variation/index.html">Perl API</a> |
<a href="vep/index.html">Variant Effect Predictor</a>
</b>
<br /><br />

<h1>Ensembl Variation - Sources Documentation</h1>

<h2>List of Variation sources for each species - $ecaption $e_versio</h2>

};


##############
### Footer  ##
##############
my $html_footer = qq{
<div style="border:1px solid #000;padding:5px;width:400px;margin-top:40px">
	<b>Legend</b><br />
	<table>
		<tr>
			<td style="width:5px;background-color:#090"></td>
			<td style="color:#090">New version of the data source in this release for the species</td>
		</tr>
		<tr>
			<td style="width:5px;background-color:#00F"></td>
			<td style="color:#00F">New data source in this release for the species</td>
		</tr>
	</table>
</div>
</body>
</html>};


my $html_content = '';

foreach my $hostname (@hostnames) {

	my $sql = qq{SHOW DATABASES LIKE '%variation_$e_version%'};
	my $sth = get_connection_and_query($database, $hostname, $sql);

	# loop over databases
	while (my ($dbname) = $sth->fetchrow_array) {
		next if ($dbname =~ /^master_schema/);
		print $dbname;
		$dbname =~ /^(.+)_variation/;
		my $s_name = $1;

		if ($etype) { # EG site - need to filter out species
			my $img_thumb = sprintf qq{eg-plugins/%s/htdocs/img/species/thumb_%s.png}, $etype, ucfirst($s_name);
			#	print "- checking for $img_thumb ... ";
		  if (! -e $img_thumb) {
				print "\t... skipping \n";
				next;
		  } 
		}
		print "\n";
		# Get list of sources from the new databases
		my $sql2 = qq{SELECT name, version, description, url, type FROM source};
		my $sth2 = get_connection_and_query($dbname, $hostname, $sql2);
		$sth2->bind_columns(\$source,\$s_version,\$s_description,\$s_url,\$s_type);
		
		
		# Previous database (and sources)
		my $p_version = $e_version-1;
		my $sql3 = qq{SHOW DATABASES LIKE '%$s_name\_variation_$p_version%'};
		my $sth3 = get_connection_and_query($database, $previous_host, $sql3);
		my $p_dbname = $sth3->fetchrow_array;
		
		my %p_list;
		my $is_new_species = 0;
		if ($p_dbname) {
			my $sql4 = qq{SELECT name, version FROM source};
			my $sth4 = get_connection_and_query($p_dbname, $previous_host, $sql4);
			while (my @p = $sth4->fetchrow_array) {
				$p_list{$p[0]} = $p[1];
			}
		}
		else {
			$is_new_species = 1;
		}
		
		$html_content .= '<br />' if ($start == 1);
		$html_content .= source_table($s_name,$sth2,$is_new_species,\%p_list);
		
		$start = 1 if ($start == 0);
  }
}

## HTML/output file ##
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML $html_header."\n";
print HTML $html_content."\n";
print HTML $html_footer."\n";
close(HTML);


sub source_table {
	my $name   = shift;
	my $sth    = shift;
	my $is_new = shift;
	my $p_list = shift;
	my $species = $name;
	$species =~ s/_/ /;
 	$species =~ /^(\w)(.+)$/;
	$species = uc($1).$2;
	my $s_name = $species;
	$s_name =~ s/\s/_/g;
	my $new_species = ($is_new == 1) ? qq{<td style="padding-left:20px;color:#00F;font-weight:bold">New species!</td>} : '';
	my $html = qq{
	<table id="$name"><tr style="vertical-align:middle">
		<td style="padding-left:0px"><img src="${server_name}/img/species/thumb_$s_name.png" alt="$species" /></td>
		<td style="padding-left:10px"><h3>$species</h3></td>$new_species
	</tr></table>
	};
	
	$html .= qq{
			<table class="ss" style="width:60%">
				<tr><th colspan="2">Source</th><th>Version</th><th>Description</th><th></th></tr>
		 };
	
	my $bg = 1;
	my @p_sources = keys(%{$p_list});
	
	
	# Chip headers
	my $cbg = 1;
	my $chip_table;
	my $chip_header .= qq{<br />
			<table class="ss" style="width:60%">
				<tr><th colspan="2">Chip Source</th><th>Version</th><th>Description</th><th></th></tr>
		 };
		 
							 
	while (my @a = $sth->fetchrow_array) {
	
		# Check if the source or its version is new
		my $s_new      = '';
		my $s_new_type = '';
		my $s_header   = '<td style="width:4px;padding:0px;margin:0px';
		if (!grep {$_ eq $source} @p_sources) {
			$s_new_type = 'source';
		}
		elsif ($p_list->{$source} ne $s_version){
			$s_new_type = 'version';
		}
		
		if ($s_new_type) {
	  	$s_new = '<span style="color:'.$colors{$s_new_type}.'">New '.$s_new_type.'</span>' if ($s_new_type);
			$s_header .= ';background-color:'.$colors{$s_new_type};
		}
			
		$s_header .= '"></td>';

		
		# Display
		if ($s_url) {
			$source = qq{<a href="$s_url">$source</a>};
		}
		
		$s_version = format_version($s_version);
		
		
		# Is chip ?
		if ($s_type eq 'chip') {
			$chip_table .= $chip_header if (!$chip_table);
			$chip_table .= qq{
				<tr class="bg$cbg">
					$s_header
					<td>$source</td>
					<td style="text-align:center">$s_version</td>
					<td>$s_description</td>
					<td style="text-align:center;width:80px">$s_new</td>
				</tr>
			};
			if ($cbg == 1) { $cbg = 2; }
			else { $cbg = 1; }
		}
		else {
			$html .= qq{
				<tr class="bg$bg">
					$s_header
					<td>$source</td>
					<td style="text-align:center">$s_version</td>
					<td>$s_description</td>
					<td style="text-align:center;width:80px">$s_new</td>
				</tr>
			};
			if ($bg == 1) { $bg = 2; }
			else { $bg = 1; }
		}
	}
	
	$html .= qq{</table>};
	$html .= qq{$chip_table</table>} if ($chip_table);
	
	return $html;
}


sub format_version {
	my $version = shift;
	
	# e.g. 20110513
	if ($version =~ /(20\d{2})(\d{2})(\d{2})/) {
		$version = "$3/$2/$1";
	}
	# e.g. 201105
	elsif ($version =~ /(20\d{2})(\d{2})/) {
		$version = "$2/$1";
	}
	# e.g. 110408
	elsif ($version =~ /(\d{2})(\d{2})(\d{2})/) {
		$version = "$3/$2/20$1";
	}
	# e.g. 20121 (HGMD data version)
	elsif ($version =~ /^(20\d{2})(\d)$/) {
		$version = "$1.$2";
	}
	elsif ($version eq '') {
		$version = '-';
	}
	
	return $version;
}

# Connects and execute a query
sub get_connection_and_query {
	my $dbname = shift;
	my $hname   = shift;
	my $sql    = shift;
	
	my ($host, $port) = split /\:/, $hname;

	# DBI connection 
	my $dsn = "DBI:mysql:$dbname:$host:$port";
	my $dbh = DBI->connect($dsn, $login, $pswd) or die "Connection failed";

	my $sth = $dbh->prepare($sql);
	$sth->execute;
	return $sth;
}


sub usage {
	
  print qq{
  Usage: perl sources2html.pl [OPTION]
  
  Put all variation sources, for each species, into an HTML document.
	
  Options:

    -help           Print this message
      
    -v							Ensembl version, e.g. 65 (Required)
    -o              An HTML output file name (Required)			
    -phost          Host name where the previous databases are stored, e.g. ensembldb.ensembl.org	(Required)
    -hlist          The list of host names where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1, ensembldb.ensembl.org2 (Required)
    -site           The URL of the website (optional)
    -etype          The type of Ensembl, e.g. Plant (optional)
  } . "\n";
  exit(0);
}

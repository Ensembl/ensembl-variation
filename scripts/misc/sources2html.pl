# Script to generate an HTML page containing the variation sources of each species

use Bio::EnsEMBL::Registry;
use DBI;
use strict;
use POSIX;
use Getopt::Long;

###############
### Options ###
###############
my ($e_version,$html_file,$source,$s_version,$s_description,$s_url,$help);

usage() if (!scalar(@ARGV));
 
GetOptions(
    'v=i' => \$e_version,
    'o=s' => \$html_file,
		'help!' => \$help
);

if (!$e_version) {
	print "> Error! Please give an Ensembl version, using the option '-v' \n";
	usage();
}
if (!$html_file) {
	print "> Error! Please give an output file using the option '-o'\n";
	usage();
}

usage() if ($help);


my @hostnames = ('ens-staging1','ens-staging2');

# Settings
my $database = "";
my $port = "3306";
my $login = "ensro";
my $pswd = "";
my $sep = "\t";
my $start = 0;


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
<h1>Ensembl Variation Sources Documentation</h1>

<a href="/info/docs/variation/index.html">About Variation Data</a> | 
<a href="/info/docs/variation/database.html">Database Description</a> |
<a href="/info/docs/variation/sources_documentation.html">Variation Sources</a> |
<a href="/info/docs/variation/variation_schema.html">Variation Tables Description</a> |
<a href="/info/docs/api/variation/index.html">Perl API</a>

<h2>List of Variation sources for each species - Ensembl $e_version</h2>

};


##############
### Footer  ##
##############
my $html_footer = qq{
</body>
</html>};


my $html_content = '';

foreach my $hostname (@hostnames) {
	# DBI connection to get all phenotype annotation in th variation database
	my $dsn = "DBI:mysql:$database:$hostname:$port";
	my $dbh = DBI->connect($dsn, $login, $pswd) or die "Connection failed";

	my $sql = qq{SHOW DATABASES LIKE '%variation_$e_version%'};
	my $sth = $dbh->prepare($sql);
	$sth->execute;

	# loop over databases
	while (my ($dbname) = $sth->fetchrow_array) {
		next if ($dbname =~ /^master_schema/);
		print $dbname."\n";
		
		# SQL connection and query
		my $dsn2 = "DBI:mysql:$dbname:$hostname:$port";
		my $dbh2 = DBI->connect($dsn2, $login, $pswd) or die "Connection failed";
		
		my $sql2 = qq{SELECT name, version, description, url FROM source};
		my $sth2 = $dbh2->prepare($sql2);
		$sth2->execute;
		$sth2->bind_columns(\$source,\$s_version,\$s_description,\$s_url);
		
		$dbname =~ /^(.+)_variation/;
		my $s_name = $1;
		
		$html_content .= '<br />' if ($start == 1);
		$html_content .= source_table($s_name,$sth2);
		
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
	my $name = shift;
	my $sth = shift;
	
	my $species = $name;
  $species =~ s/_/ /;
 	$species =~ /^(\w)(.+)$/;
	$species = uc($1).$2;
	
	my $html  = qq{\n<h3 id="$name"># $species</h3>};
	   $html .= qq{
			<table class="ss" style="width:60%">
				<tr><th>Source</th><th>Version</th><th>Description</th></tr>
		 };
	
	my $bg = 1;
	
	while (my @a = $sth->fetchrow_array) {
		if ($s_url) {
			$source = qq{<a href="$s_url">$source</a>};
		}	
		$s_version = format_version($s_version);
		
		$html .= qq{
			<tr class="bg$bg">
				<td>$source</td>
				<td>$s_version</td>
				<td>$s_description</td>
			</tr>
		};
		if ($bg == 1) { $bg = 2; }
		else { $bg = 1; }
	}
	
	$html .= qq{</table>};
	
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
	elsif ($version eq '') {
		$version = '-';
	}
	
	return $version;
}


sub usage {
	
  print qq{
  Usage: perl source2html.pl [OPTION]
  
  Put all variation sources, for each species, into an HTML document.
	
  Options:

    -help           Print this message
      
		-v							Ensembl version, e.g. 63 (Required)
    -o              An HTML output file name (Required)							 
  } . "\n";
  exit(0);
}

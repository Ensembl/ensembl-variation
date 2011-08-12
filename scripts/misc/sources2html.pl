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
my $previous_host = 'ens-livemirror';

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

	my $sql = qq{SHOW DATABASES LIKE '%variation_$e_version%'};
  my $sth = get_connection_and_query($database, $hostname, $sql);

	# loop over databases
	while (my ($dbname) = $sth->fetchrow_array) {
		next if ($dbname =~ /^master_schema/);
		print $dbname."\n";
		
		# Get list of sources from the new databases
		my $sql2 = qq{SELECT name, version, description, url FROM source};
    my $sth2 = get_connection_and_query($dbname, $hostname, $sql2);
		$sth2->bind_columns(\$source,\$s_version,\$s_description,\$s_url);
		
		$dbname =~ /^(.+)_variation/;
		my $s_name = $1;
		
		# Previous database (and sources)
		my $p_version = $e_version-1;
		my $sql3 = qq{SHOW DATABASES LIKE '%$s_name\_variation_$p_version%'};
		my $sth3 = get_connection_and_query($database, $previous_host, $sql3);
		my $p_dbname = $sth3->fetchrow_array;
		
		#print "\tPREVIOUS: $p_dbname\n";
		my %p_list;
		if ($p_dbname) {
			my $sql4 = qq{SELECT name, version FROM source};
			my $sth4 = get_connection_and_query($p_dbname, $previous_host, $sql4);
			while (my @p = $sth4->fetchrow_array) {
				$p_list{$p[0]} = $p[1];
			}
		}
		
		$html_content .= '<br />' if ($start == 1);
		$html_content .= source_table($s_name,$sth2,\%p_list);
		
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
	my $p_list = shift;
	
	my $species = $name;
  $species =~ s/_/ /;
 	$species =~ /^(\w)(.+)$/;
	$species = uc($1).$2;
	
	my $s_name = $species;
	$s_name =~ s/\s/_/g;
	my $html = qq{
	<table id="$name"><tr style="vertical-align:middle">
		<td style="padding-left:0px"><img src="http://static.ensembl.org/img/species/thumb_$s_name.png" alt="$species" /></td>
		<td style="padding-left:10px"><h3>$species</h3></td>
	</tr></table>
	};
	
	$html .= qq{
			<table class="ss" style="width:60%">
				<tr><th>Source</th><th>Version</th><th>Description</th><th></th></tr>
		 };
	
	my $bg = 1;
	my @p_sources = keys(%{$p_list});
	
	my %colors = ( 'version' => '#009900',
								 'source'  => '#000066',
	             );
	while (my @a = $sth->fetchrow_array) {
	
		# Check if the source or the its version are new
		my $s_new      = '';
		my $s_new_type = '';
		
		if ($p_list->{$source} ne $s_version){
			$s_new_type = 'version';
		}
		elsif (!grep {$_ eq $source} @p_sources) {
			$s_new_type = 'source';
		}
	  $s_new = '<span style="color:'.$colors{$s_new_type}.'">New '.$s_new_type.'</span>' if ($s_new_type);
	
		# Display
		if ($s_url) {
			$source = qq{<a href="$s_url">$source</a>};
		}
		
		$s_version = format_version($s_version);
		
		
		
		$html .= qq{
			<tr class="bg$bg">
				<td>$source</td>
				<td>$s_version</td>
				<td>$s_description</td>
				<td style="text-align:center;width:80px">$s_new</td>
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

# Connects and execute a query
sub get_connection_and_query {
	my $dbname = shift;
	my $host   = shift;
	my $sql    = shift;
	
	# DBI connection 
	my $dsn = "DBI:mysql:$dbname:$host:$port";
	my $dbh = DBI->connect($dsn, $login, $pswd) or die "Connection failed";

	my $sth = $dbh->prepare($sql);
	$sth->execute;
	return $sth;
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

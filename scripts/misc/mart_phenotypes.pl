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


=begin
The script generates biomart tmp tables for variation phenotype data.
It stores gene phenotype associations in MTMP_phenotypes.
=end
=cut


use DBI;
use Getopt::Long;
use ImportUtils qw(load);
my $config = {};

GetOptions(
    $config,
	'tmpdir=s',
	'tmpfile=s',
	'host|h=s',
	'user|u=s',
	'port|P=i',
	'password|p=s',
	'db|d=s',
	'version|v=i',
	'pattern=s',
) or die "ERROR: Could not parse command line options\n";

# check options
for(qw(host user port tmpdir tmpfile)) {
	die("ERROR: $_ not defined, use --$_\n") unless defined $config->{$_};
}

my @db_list;

if(!defined($config->{db})) {
  @db_list = @{get_species_list($config, $config->{host})};
}
else {
  push @db_list, $config->{db};
}

die "ERROR: no suitable databases found on host ".$config->{host}."\n" unless scalar @db_list;

my $TMP_DIR = $config->{tmpdir};
my $TMP_FILE = $config->{tmpfile};

$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;

foreach my $db (@db_list) {
	my $dbc = DBI->connect(
		sprintf(
			"DBI:mysql(RaiseError=>1):host=%s;port=%s;db=%s",
			$config->{host}, $config->{port}, $db
		), $config->{user}, $config->{password},
	);
	$dbc->do(qq{DROP TABLE IF EXISTS MTMP_phenotype;});

	$dbc->do(qq{
		CREATE TABLE `MTMP_phenotype` (
		`mtmp_phenotype_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
		`gene_stable_id` varchar(255) NOT NULL,
		`description` varchar(255) NOT NULL,
		`source` varchar(24) NOT NULL,
		`p_value` double default NULL, 
		`strain_name`  varchar(255) default NULL,
		`strain_gender` varchar(255) default NULL,
		`external_id` varchar(255) default NULL,
		PRIMARY KEY (`mtmp_phenotype_id`),
		KEY `gene_stable_idx` (`gene_stable_id`)
		) ENGINE=MyISAM DEFAULT CHARSET=latin1;
	});
	my $strains;

	my $sth = $dbc->prepare(qq{
		SELECT individual_id, name, gender FROM individual;
	});
	$sth->execute();
	while (my $row = $sth->fetchrow_arrayref) {
		my ($individual_id, $name, $gender) = @$row;
		$strains->{$individual_id}->{name} = $name;
		$strains->{$individual_id}->{gender} = $gender;
	}
	$sth->finish();	

	my $sth = $dbc->prepare(qq{
		SELECT pf.phenotype_feature_id, pf.object_id, p.description, at.code, pfa.value, s.name
		FROM attrib_type at
		LEFT JOIN phenotype_feature_attrib as pfa ON (at.attrib_type_id = pfa.attrib_type_id)
		LEFT JOIN phenotype_feature as pf ON (pf.phenotype_feature_id = pfa.phenotype_feature_id)
		LEFT JOIN phenotype as p ON (pf.phenotype_id = p.phenotype_id)
		LEFT JOIN source as s on (s.source_id = pf.source_id)	
		WHERE pf.type = 'Gene'
		ORDER BY pf.phenotype_feature_id;
	});

	$sth->execute();
   	my ($pf_id, $gene_id, $desc, $code, $value, $source);
	$sth->bind_columns(\($pf_id, $gene_id, $desc, $code, $value, $source));	    
	my $current_pf_id;
	my $row;
    open OUT, ">$TMP_DIR/$TMP_FILE";

	while ($sth->fetch) {
		if ($pf_id ne $current_pf_id) {
			if ($current_pf_id) {
				print OUT join("\t", @{get_values($row)}), "\n";
				$row = {};
			}
		}
		$row->{pf_id}          = $pf_id;
		$row->{gene_stable_id} = $gene_id;
        $desc =~ s/\&|,//g;
		$row->{description}    = $desc;
		$row->{source}         = $source;
		if ($code eq 'strain_id') {
			$row->{strain_name}   = $strains->{$value}->{name};
			$row->{strain_gender} = $strains->{$value}->{gender}; 
		}
		if ($code eq 'p_value')	 {
			$row->{p_value} = $value;
		}
		if ($code eq 'external_id') {
			$row->{external_id} = $value;
		}
		$current_pf_id = $pf_id;
    }
        
    $sth->finish;
	# print last phenotype feature:
	print OUT join("\t", @{get_values($row)}), "\n";

	close OUT;
	load($dbc, 'MTMP_phenotype');
}

sub get_values {
	my $row = shift;
	$row->{strain_name} ||= '\N';		
	$row->{strain_gender} ||= '\N';		
	$row->{p_value} ||= '\N';
	$row->{external_id} ||= '\N';
	my @values = (
		$row->{pf_id}, 
		$row->{gene_stable_id}, 
		$row->{description}, 
		$row->{source}, 
		$row->{p_value}, 
		$row->{strain_name}, 
		$row->{strain_gender},
		$row->{external_id},);
	return \@values; 	
}

sub get_species_list {
	my $config = shift;
	my $host   = shift;
	
	# connect to DB
	my $dbc = DBI->connect(
		sprintf(
			"DBI:mysql(RaiseError=>1):host=%s;port=%s;db=mysql",
			$host,
			$config->{port}
		), $config->{user}, $config->{password}
	);
	
	my $version = $config->{version};
	
	my $sth = $dbc->prepare(qq{
		SHOW DATABASES LIKE '%\_variation\_$version%'
	});
	$sth->execute();
	
	my $db;
	$sth->bind_columns(\$db);
	
	my @dbs;
	push @dbs, $db while $sth->fetch;
	$sth->finish;
	
	# remove master and coreexpression
	@dbs = grep {$_ !~ /master|express/} @dbs;
	
	# filter on pattern if given
	my $pattern = $config->{pattern};
	@dbs = grep {$_ =~ /$pattern/} @dbs if defined($pattern);
	
	# remove version, build
	#$_ =~ s/^([a-z]+\_[a-z]+)(.+)/$1/ for @dbs;
	
	return \@dbs;
}

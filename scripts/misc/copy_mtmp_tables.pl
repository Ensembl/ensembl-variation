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


use DBI;
use Getopt::Long;

my $config = {};

GetOptions(
    $config,
	'dir|d=s',
	'tables|t=s',
	'pattern|p=s'
) or die "ERROR: Could not parse command line options\n";

# check host and user
die "ERROR: You are not logged in as user mysqlens\n" unless $ENV{USER} eq 'mysqlens';
#die "ERROR: You must run this script on the database server machine on which the MTMP table database resides\n" unless $ENV{HOST} eq 'ens-variation2';

# parse tables
$config->{tables} ||= 'MTMP_population_genotype,subsnp_map,tmp_sample_genotype_single_bp';
my %t = map {$_ => 1} split /\,/, $config->{tables};
$config->{tables} = \%t;

# set up to connect
$config->{user} = 'ensro';
$config->{password} = '';
$config->{port} = 3306;

# set up dir
$config->{dir} ||= '/mysql/data_3306/databases/production_archive/';

# check dir exists
die "ERROR: Directory ".$config->{dir}." not found, perhaps you need to specify with --dir?\n" unless -e $config->{dir};

# get file list
opendir DIR, $config->{dir};
my @file_list = readdir DIR;

# get species lists
my $hash = {};

foreach my $host(qw(ens-staging ens-staging2)) {

	my @list = @{get_species_list($config, $host)};
	
	foreach my $db(@list) {
		my $species = $db;
		$species =~ s/^([a-z]+\_[a-z]+)(.+)/$1/;
        next if ($species =~ /homo_sapiens/);	
		# connect to DB
		my $dbc = DBI->connect(
			sprintf(
				"DBI:mysql(RaiseError=>1):host=%s;port=%s;db=%s",
				$host,
				$config->{port},
				$db
			), $config->{user}, $config->{password}
		);
		
		# find existing tables
		my $sth = $dbc->prepare(qq{
			SHOW TABLES LIKE 'MTMP_%'
		});
		$sth->execute();
		
		my $table;
		$sth->bind_columns(\$table);
		
		my @tables;
		push @tables, $table while $sth->fetch;
		$sth->finish;
		
		$sth = $dbc->prepare(qq{
			SHOW TABLES LIKE '%subsnp_map%'
		});
		$sth->execute();
		
		$sth->bind_columns(\$table);
		
		push @tables, $table while $sth->fetch;
		$sth->finish;
		
		$sth = $dbc->prepare(qq{
			SHOW TABLES LIKE '%tmp_sample_genotype_single_bp%'
		});
		$sth->execute();
		
		$sth->bind_columns(\$table);
		
		push @tables, $table while $sth->fetch;
		$sth->finish;
		
		my %table_hash = map {$_ => 1} grep {$config->{tables}->{$_}} @tables;
		
		# find tables to copy
		my @species_files = grep {$_ =~ /^$species/} @file_list;
		my (@filtered, %map);
		
		foreach my $f(@species_files) {
			# strip off filetype
			$f =~ s/\..+$//g;
			
			# special case for tmp_gt, filename is too long for MySQL 64-char table name limit
			my $match = $f;
			$match =~ s/tmp_gt/tmp_sample_genotype_single_bp/;
			$match =~ s/pop_gt/population_genotype/;
			
			my $ok = 0;
			foreach my $t(keys %{$config->{tables}}) {
				if($match =~ /$t/) {
					$ok = 1;
          $f =~ /(.+\_variation_)(\d+)(.+)/;
          my $version = $2;
					$map{$t}{$version} = $f;
				}
			}
			
			next unless $ok;
			
			push @filtered, $f;
		}
		
		if(!scalar keys %map) {
			print "No table files found for $species, skipping\n";
			next;
		}
		
		foreach my $t(keys %map) {
			if($table_hash{$t}) {
				print "Table $t already exists on $host for $species, skipping\n";
				next;
			}
		  my $latest_version = (sort {$b <=> $a} keys %{$map{$t}})[0];	
			print "$species: Copying $t\n";
			
			for my $ext(qw(frm MYD MYI)) {
				system "scp ".$config->{dir}."/$map{$t}{$latest_version}\.$ext $host:/mysql/data_3306/databases/$db/$t\.$ext";
#				print "scp ".$config->{dir}."/$map{$t}{$latest_version}\.$ext $host:/mysql/data_3306/databases/$db/$t\.$ext\n";
			}
		}
	}
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

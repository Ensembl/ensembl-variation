#!/usr/bin/perl

use DBI;
use Getopt::Long;

my $config = {};

GetOptions(
    $config,
	'dir|d=s',
	'tables|t=s',
) or die "ERROR: Could not parse command line options\n";

# check host and user
die "ERROR: You are not logged in as user mysqlens\n" unless $ENV{USER} eq 'mysqlens';
die "ERROR: You must run this script on the database server machine on which the MTMP table database resides\n" unless $ENV{HOST} eq 'ens-variation2';

# parse tables
$config->{tables} ||= 'MTMP_allele,MTMP_population_genotype';
my %t = map {$_ => 1} split /\,/, $config->{tables};
$config->{tables} = \%t;

# set up to connect
$config->{user} = 'ensro';
$config->{password} = '';
$config->{port} = 3306;

# set up dir
$config->{dir} ||= '/mysql/data_3306/databases/mart_mtmp_tables/';

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
		
		my %table_hash = map {$_ => 1} grep {$config->{tables}->{$_}} @tables;
		
		# find tables to copy
		my @species_files = grep {$_ =~ /^$species/} @file_list;
		my (@filtered, %map);
		
		foreach my $f(@species_files) {
			
			# strip off filetype
			$f =~ s/\..+$//g;
			
			my $ok = 0;
			foreach my $t(keys %{$config->{tables}}) {
				if($f =~ /$t/) {
					$ok = 1;
					$map{$t} = $f;
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
			
			print "$species: Copying $t\n";
			
			for my $ext(qw(frm MYD MYI)) {
				system "scp ".$config->{dir}."/$map{$t}\.$ext $host:/mysql/data_3306/databases/$db/$t\.$ext\n";
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

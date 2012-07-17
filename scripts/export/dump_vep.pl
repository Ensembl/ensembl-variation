#!/usr/bin/perl

use Getopt::Long;
use DBI;

my $config = {};

my %special_options = (
	'homo_sapiens' => ' --sift b --polyphen b --regulatory',
	'mus_musculus' => ' --regulatory',
);

GetOptions(
    $config,
	'hosts|h=s',
	'user|u=s',
	'port|P=i',
	'password|p=s',
	'dir=s',
	'command=s',
	'mem=i',
	'pattern=s',
	'queue=s',
	'version=i',
	'overwrite',
) or die "ERROR: Could not parse command line options\n";

# set defaults
$config->{hosts}    ||= 'ens-staging,ens-staging2';
$config->{user}     ||= 'ensro';
$config->{port}     ||= 3306;
$config->{dir}      ||= $ENV{'HOME'}.'/.vep/';
$config->{command}  ||= 'perl variant_effect_predictor.pl --build all';
$config->{mem}      ||= 12000000;
$config->{queue}    ||= 'normal';

# check dir exists
die "ERROR: Directory ".$config->{dir}." does not exist\n" unless -e $config->{dir};

# check version defined
die "ERROR: No Ensembl DB version defined - use --version [version]\n" unless defined($config->{version});

# check command supplied looks sensible
die "ERROR: Supplied command doesn't look right, it should look something like:\nperl -I /include/perl/libs/ /path/to/variant_effect_predictor.pl --build all" unless
	$config->{command} =~ /\-build (all|\w+)/ &&
	$config->{command} =~ /variant_effect_predictor/;

# parse hosts
$config->{hosts} = [split /\,/, $config->{hosts}];

foreach my $host(@{$config->{hosts}}) {
	debug("Getting species list for host $host");
	
	my $species_list = get_species_list($config, $host);
	
	debug("Found ".(scalar @$species_list)." valid databases\n");
	
	foreach my $species(@$species_list) {
		debug("Running VEP dump for $species");
		dump_vep($config, $host, $species);
		
		debug("Compressing dump file");
		tar($config, $species);
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
		SHOW DATABASES LIKE '%\_core\_$version%'
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
	$_ =~ s/^([a-z]+\_[a-z]+)(.+)/$1/ for @dbs;
	
	return \@dbs;
}

sub dump_vep {
	my $config  = shift;
	my $host    = shift;
	my $species = shift;
	
	# check if dir exists
	if(!defined($config->{overwrite}) && -e $config->{dir}.'/'.$species.'/'.$config->{version}) {
		debug("Dump directory found for $species, skipping (use --overwrite to overwrite)\n");
		return;
	}
	
	my $command = join " ", (
		sprintf(
			'bsub -K -J %s -M %i -R"select[mem>%i] rusage[mem=%i]" -q %s -o %s -e %s',
			$species,
			$config->{mem},
			$config->{mem} / 1000,
			$config->{mem} / 1000,
			$config->{queue},
			$config->{dir}.'/'.$species.'_vep_dump.farmout',
			$config->{dir}.'/'.$species.'_vep_dump.farmerr',
		),
		$config->{command},
		defined $special_options{$species} ? $special_options{$species} : "",
		"--no_adaptor",
		"--species ".$species,
		"--host ".$host,
		"--user ".$config->{user},
		"--port ".$config->{port},
		"--dir ".$config->{dir},
		defined $config->{password} ? "--password ".$config->{pass} : "",
	);
	
	debug("Use \"tail -f ".$config->{dir}.'/'.$species.'_vep_dump.farmout'."\" to check progress");
	system($command);
}

sub tar {
	my $config  = shift;
	my $species = shift;
	
	my $tar_file = $config->{dir}."/".$species."_vep_".$config->{version}.".tar.gz";
	
	# check if tar exists
	if(!defined($config->{overwrite}) && -e $tar_file) {
		debug("Dump file found for $species, skipping (use --overwrite to overwrite)\n");
		return;
	}
	
	# check dir exists
	my $root_dir = $config->{dir};
	my $sub_dir  = $species."/".$config->{version};
	
	die("ERROR: VEP dump directory $root_dir/$sub_dir not found") unless -e $root_dir.'/'.$sub_dir;
	
	my $command = "tar -cz -C $root_dir -f $tar_file $sub_dir";
	system($command);
}

# gets time
sub get_time() {
    my @time = localtime(time());

    # increment the month (Jan = 0)
    $time[4]++;

    # add leading zeroes as required
    for my $i(0..4) {
        $time[$i] = "0".$time[$i] if $time[$i] < 10;
    }

    # put the components together in a string
    my $time =
        ($time[5] + 1900)."-".
        $time[4]."-".
        $time[3]." ".
        $time[2].":".
        $time[1].":".
        $time[0];

    return $time;
}

# prints debug output with time
sub debug {
    my $text = (@_ ? (join "", @_) : "No message");
    my $time = get_time;
    
    print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}
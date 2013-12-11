use strict;
use warnings;

use DBI;
use FileHandle;
use Getopt::Long;
use JSON;

my $config = {};
GetOptions(
	$config,
	'host=s',
	'user=s',
	'port=s',
    'version=i',
	'save_config_file=s',
	'db_name=s',	
	'debug',
) or die 'Failed to parse command line arguments';

$config->{user} ||= 'ensro';
$config->{port} ||= 3306;

die 'Argument \'host\' is missing' unless ($config->{host});
die 'Argument \'save_config_file\' is missing' unless ($config->{save_config_file});

my $types = {
	GENERIC         => 'generic',
	SVS             => 'structural_variations',
	FAILED          => 'failed',
	INC_CONS        => 'incl_consequences',
	SOMATIC         => 'somatic',
	SOMATIC_INC_CONS => 'somatic_incl_consequences',	
};

my $populations = {
	'homo_sapiens' => {
		'1000GENOMES:phase_1_AFR'   => '1000GENOMES-phase_1_AFR',
		'1000GENOMES:phase_1_AMR'   => '1000GENOMES-phase_1_AMR',
		'1000GENOMES:phase_1_ASN'   => '1000GENOMES-phase_1_ASN',
		'1000GENOMES:phase_1_EUR'   => '1000GENOMES-phase_1_EUR',
		'CSHL-HAPMAP:HAPMAP-ASW'    => 'CSHL-HAPMAP-HAPMAP-ASW',
		'CSHL-HAPMAP:HAPMAP-CHB'    => 'CSHL-HAPMAP-HAPMAP-CHB',
		'CSHL-HAPMAP:HAPMAP-CHD'    => 'CSHL-HAPMAP-HAPMAP-CHD',
		'CSHL-HAPMAP:HAPMAP-GIH'    => 'CSHL-HAPMAP-HAPMAP-GIH',
		'CSHL-HAPMAP:HAPMAP-LWK'    => 'CSHL-HAPMAP-HAPMAP-LWK',
		'CSHL-HAPMAP:HAPMAP-MEX'    => 'CSHL-HAPMAP-HAPMAP-MEX',
		'CSHL-HAPMAP:HAPMAP-MKK'    => 'CSHL-HAPMAP-HAPMAP-MKK',
		'CSHL-HAPMAP:HAPMAP-TSI'    => 'CSHL-HAPMAP-HAPMAP-TSI',
		'CSHL-HAPMAP:HapMap-CEU'    => 'CSHL-HAPMAP-HapMap-CEU',
		'CSHL-HAPMAP:HapMap-HCB'    => 'CSHL-HAPMAP-HapMap-HCB',
		'CSHL-HAPMAP:HapMap-JPT'    => 'CSHL-HAPMAP-HapMap-JPT',
		'CSHL-HAPMAP:HapMap-YRI'    => 'CSHL-HAPMAP-HapMap-YRI',
		'ESP6500:African_American'  => 'ESP6500-African_American',
		'ESP6500:European_American' => 'ESP6500-European_American',
	},
};

my $individuals = {
	'homo_sapiens' => {
		'Venter' => 'Venter',
		'Watson' => 'Watson',
	},
};

my $sets = {
	'homo_sapiens' => {
		'clinically associated'             => 'clinically_associated',
		'all phenotype-associated variants' => 'phenotype_associated',
	},	
};

get_species_with_variation_db($config);

if ($config->{debug}) {
	my $all_species = $config->{species};
	foreach my $species (keys %$all_species) {
		print STDERR $species, "\n";	
	}
} else {
	my $species_with_sift_predictions = get_species_with_sift_predictions($config);	
	my $species_with_aa               = get_species_with_aa($config);
	my $species_with_svs_data         = get_species_with_svs_data($config);
#my $species_with_phenotype_data = get_species_with_phenotype_data($config);


	my $dump_configuration = {};

	my $species = $config->{species};
	foreach my $species (keys %$species) {
		my $parameter = 'evidence_values';
		my $protein_function_prediction = 'incl_consequences';
		if ($species_with_aa->{$species}) {
			$parameter .= ',ancestral_allele';	
		}
		if ($species_with_sift_predictions->{$species}) {
			$protein_function_prediction .= ',sift';	
		}
		if ($species eq 'homo_sapiens') {
			$parameter .= ',global_maf,clinical_significance';
			$protein_function_prediction .= ',polyphen';

			$dump_configuration->{$species}->{individuals} = $individuals->{$species};
			$dump_configuration->{$species}->{populations} = $populations->{$species};
			$dump_configuration->{$species}->{sets} = $sets->{$species};
			$dump_configuration->{$species}->{$types->{SOMATIC}} = $parameter . ',somatic';
			$dump_configuration->{$species}->{$types->{SOMATIC_INC_CONS}} = $protein_function_prediction . ',somatic';
		}
		$dump_configuration->{$species}->{$types->{GENERIC}} = $parameter; 
		$dump_configuration->{$species}->{$types->{INC_CONS}} = $protein_function_prediction;
		$dump_configuration->{$species}->{$types->{FAILED}} = 'failed'; 

		if ($species_with_svs_data->{$species}) {
			$dump_configuration->{$species}->{$types->{SVS}} = 'structural_variations';
		}
	}


	my $fh_config = FileHandle->new($config->{save_config_file}, 'w');
	my $json = JSON->new->allow_nonref;
	print $fh_config $json->encode($dump_configuration);
	$fh_config->close();
}
sub get_species_with_variation_db {
	my $cofig = shift;
	my $db_name_2_species_name = {};
	if ($config->{db_name}) {
		foreach my $db_name (split(',', $config->{db_name})) {
			die "$db_name is not a valid name for a variation database" unless(is_variation_database($db_name));
			if ($config->{version}) {
				my $version = $config->{version};
				if ($db_name =~ m/$version/) {
					$db_name_2_species_name->{$db_name} = get_species_name($db_name);
				} else {
					warn "Version $version differs for $db_name";
				}
			} else {
				$db_name_2_species_name->{$db_name} = get_species_name($db_name);
			}
		}
	} else {
		my $user = $config->{user};
		my $port = $config->{port};
		foreach my $host (split(',', $config->{host})) {
			my $dbc = DBI->connect("DBI:mysql:host=$host;port=$port;user=$user", {RaiseError => 1});
			my $sth = $dbc->prepare(qq{show databases like '%variation%'});
			$sth->execute();
			my $db;
			$sth->bind_columns(\$db);
			while ($sth->fetch) {
				if (is_variation_database($db)) {
					if ($config->{version}) {
						my $version = $config->{version};
						if ($db =~ m/$version/) {
							$db_name_2_species_name->{$db} = get_species_name($db);
						}
					} else {
						$db_name_2_species_name->{$db} = get_species_name($db);
					}
				}
			}	
			$sth->finish();
		}
	}
	$config->{db2species} = $db_name_2_species_name;
	my %species = map {$_ => 1 } values %$db_name_2_species_name;
	$config->{species} = \%species;
}

sub is_variation_database {
	my $db_name = shift;
	return ($db_name =~ m/(\D+\_\D+)(\_variation\_)(\d+\_\d+)/); # match e.g. homo_sapiens_variation_74_37
}

sub get_species_name {
	my $db_name = shift;
	if ($db_name =~ m/(\D+\_\D+)(\_variation\_)(\d+\_\d+)/) { # match e.g. homo_sapiens_variation_74_37
		return $1;
	} else {
		die "$db_name not a valid name for a variation database";
	}
}

sub get_species_with_sift_predictions {
	my $config = shift;
	my $query = qq{select count(*) from protein_function_predictions;};
	return query_database($config, $query);
}

sub get_species_with_svs_data {
	my $config = shift;
	my $query = qq{select count(*) from structural_variation;};
	return query_database($config, $query);
}

sub get_species_with_phenotype_data {
	my $config = shift;
	my $query = qq{select count(*) from phenotype_feature;};
	return query_database($config, $query);
}

sub get_species_with_aa {
	my $config = shift;
	my $query = qq{select variation_id from variation where ancestral_allele is not null limit 1;};
	return query_database($config, $query);
}

sub query_database {
	my $config = shift;
	my $query = shift;
	my $hosts = $config->{host};
	my $user = $config->{user};
	my $port = $config->{port};
	my $species_names = {};	
	foreach my $host (split(',', $hosts)) {
		foreach my $db_name (keys %{$config->{db2species}})	 {
			my $dbc = DBI->connect("DBI:mysql:database=$db_name;host=$host;port=$port;user=$user", {RaiseError => 1});
			my $sth = $dbc->prepare($query);
			$sth->execute();
			while (my @row = $sth->fetchrow_array) {
				my $count = $row[0];
				if ($count > 0) {
					my $species_name = $config->{db2species}->{$db_name};
					$species_names->{$species_name} = 1;				
				}
			}
			$sth->finish();
		}
	}
	return $species_names;
}

# structural_variation
# somatic
# incl consequences: protein info: sift, polyphen
# evidence, clinical_significance, ancestral_allele, minor_allele_freq, validation_status
# populations:
# individuals:
# sets: phenotypes, clinically_associated

#!/usr/bin/perl

=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME
import_vcf.pl - imports variations from a VCF file into an Ensembl variation DB

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Getopt::Long;
use FileHandle;
use Data::Dumper;
use Time::HiRes qw(gettimeofday tv_interval);
use ImportUtils qw(debug load);
use FindBin qw( $Bin );

use constant DISTANCE => 100_000;
use constant MAX_SHORT => 2**16 -1;

my %Printable = ( "\\"=>'\\', "\r"=>'r', "\n"=>'n', "\t"=>'t', "\""=>'"' );



# COMMAND LINE OPTIONS
######################

# get command-line options
my ($in_file, $species, $registry_file, $help, $host, $database, $port, $user, $password, $source, $source_desc, $population, $flank_size, $TMP_DIR, $TMP_FILE, $skip_multi, $use_gp, $sample_prefix, $variation_prefix, $disable_keys, $include_tables, $merge_vfs, $skip_tables, $compressed_only, $only_existing, $merge_alleles, $new_var_name, $chrom_regexp, $check_synonyms, $force_multi, $force_no_var);

my $args = scalar @ARGV;

GetOptions(
	'input_file=s'   => \$in_file,
	'species=s'      => \$species,
	'registry=s'     => \$registry_file,
	'host=s'         => \$host,
	'database|db=s'  => \$database,
	'user=s'         => \$user,
	'password=s'     => \$password,
	'port=i'         => \$port,
	'help'           => \$help,
	'source=s'       => \$source,
	'source_desc=s'  => \$source_desc,
	'population=s'   => \$population,
	'flank=s'        => \$flank_size,
	'tmpdir=s'       => \$TMP_DIR,
	'tmpfile=s'      => \$TMP_FILE,
	'skip_multi'     => \$skip_multi,
	'gp'             => \$use_gp,
	'ind_prefix=s'   => \$sample_prefix,
	'var_prefix=s'   => \$variation_prefix,
	'disable_keys'   => \$disable_keys,
	'tables=s'       => \$include_tables,
	'skip_tables=s'  => \$skip_tables,
	'merge_vfs'      => \$merge_vfs,
	'only_existing'  => \$only_existing,
	'merge_alleles'  => \$merge_alleles,
	'create_name'    => \$new_var_name,
	'chrom_regexp=s' => \$chrom_regexp,
	'check_synonyms' => \$check_synonyms,
	'force_multi'    => \$force_multi,
	'force_no_var'   => \$force_no_var,
);


# print usage message if requested or no args supplied
if(defined($help) || !$args) {
	&usage;
	exit(0);
}

# set defaults
$species ||= "human";
$flank_size ||= 200;
$port ||= 3306;

# set default list of tables to write to
my $tables = {
	'variation'                       => 1,
	'variation_feature'               => 1,
	'flanking_sequence'               => 1,
	'allele'                          => 1,
	'population_genotype'             => 1,
	'compressed_genotype_var'         => 1,
	'individual_genotype_multiple_bp' => 0,
	'compressed_genotype_region'      => 0,
	'sample'                          => 1,
	'population'                      => 1,
	'individual'                      => 1,
	'individual_population'           => 1,
};

# override this with options if provided
if(defined($include_tables)) {
	
	# reset
	$tables->{$_} = 0 foreach keys %$tables;
	
	# set include tables
	foreach my $table(split /\,/, $include_tables) {
		$tables->{$table} = 1 if defined($tables->{$table});
	}
}

if(defined($skip_tables)) {
	
	# set skip tables
	foreach my $table(split /\,/, $skip_tables) {
		$tables->{$table} = 0 if defined($tables->{$table});
	}
}

$compressed_only = 1 if scalar (grep {$tables->{$_}} keys %$tables) == 1 && $tables->{compressed_genotype_region};
my $get_var_from_vf = 1;

# special case for sample, we also want to set the other sample tables to on
if($tables->{sample}) {
	$tables->{$_} = 1 for qw/population individual individual_population/;
}

# force population if user wants allele or population_genotype
if($tables->{allele} || $tables->{population_genotype} || $tables->{compressed_genotype_region} || $tables->{individual_genotype_multiple_bp} || $tables->{compressed_genotype_var}) {
	$tables->{population} = 1;
	$tables->{sample} = 1;
}

# check that at least one has been set
die "ERROR: no tables left included\n" unless grep {$tables->{$_}} keys %$tables;



die "ERROR: tmpdir not specified\n" if !defined $TMP_DIR && $tables->{compressed_genotype_region};
$TMP_FILE ||= 'compress.txt';
$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;



## DB CONNECTION
################

my $dbVar;

if(defined($database)) {
	$dbVar = DBI->connect( "DBI:mysql(RaiseError=>1):host=$host;port=$port;db=$database", $user, $password );
}
else {
	
	# get registry
	my $reg = 'Bio::EnsEMBL::Registry';
	
	if(defined($host) && defined($user)) {
		$reg->load_registry_from_db(-host => $host, -user => $user, -pass => $password);
	}
	
	else {
		if(-e $registry_file) {
			$reg->load_all($registry_file);
		}
		else {
			die "ERROR: could not read from registry file $registry_file\n";
		}
	}

	# connect to DB
	my $vdba = $reg->get_DBAdaptor($species,'variation')
		|| usage( "Cannot find variation db for $species in $registry_file" );
	$dbVar = $vdba->dbc->db_handle;

	debug("Connected to database ", $vdba->dbc->dbname, " on ", $vdba->dbc->host, " as user ", $vdba->dbc->username);
}



# INPUT FILE HANDLE
###################

# define the filehandle to read input from
my $in_file_handle = new FileHandle;

if(defined($in_file)) {

	# check defined input file exists
	die("ERROR: Could not find input file ", $in_file, "\n") unless -e $in_file;
	
	if ($in_file =~ /\.gz$/){
		$in_file_handle->open("zcat ". $in_file . " | " ) or die("ERROR: Could not read from input file ", $in_file, "\n");
	}
	elsif ($in_file =~ /\.vcf$/){
		$in_file_handle->open( $in_file ) or die("ERROR: Could not read from input file ", $in_file, "\n");
	}
	else{
		die "ERROR: Not sure how to handle file type of ", $in_file, "\n";
	}

	debug("Reading from file ", $in_file);
}

# no file specified - try to read data off command line
else {
	$in_file_handle = 'STDIN';
	debug("Attempting to read from STDIN");
}



# DB PREP
#########

# get seq_region_id hash
my $seq_region_ids = &get_seq_region_ids($dbVar);
die("ERROR: seq_region not populated\n") unless scalar keys %$seq_region_ids;

# get/set source_id
die("ERROR: no source specified\n") if !(defined $source) && !$only_existing;
my $source_id = &get_source_id($dbVar, $source, $source_desc);

# now do population
my $pop_id;
if($tables->{population}) {
	die("ERROR: no population specified\n") unless defined $population;
	$pop_id = &get_population_id($dbVar, $population);
}

# disable keys if requested
if($disable_keys) {
	debug("Disabling keys");
	
	foreach my $table(grep {$tables->{$_}} keys %$tables) {
		$dbVar->do(qq{ALTER TABLE $table DISABLE KEYS;})
	}
}

# get genotype codes
my $gca = Bio::EnsEMBL::Registry->get_adaptor($species, 'variation', 'genotypecode');
my %gt_codes = map {(join "|", @{$_->genotype}) => $_->dbID} @{$gca->fetch_all()};

# get allele codes
my $sth = $dbVar->prepare(qq{SELECT allele_code_id, allele FROM allele_code});
$sth->execute;
my %allele_codes = map {$_->[1] => $_->[0]} @{$sth->fetchall_arrayref};
$sth->finish();



# SET UP VARIABLES
##################

my (
	%headers,
	$first_sample_col,
	$sample_ids,
	$region_start,
	$region_end,
	$prev_region_end,
	$prev_seq_region,
	$genotypes,
	$var_counter,
	$start_time,
	$skipped,
);

$start_time = time();



# MAIN FILE LOOP
################

# read the file
while(<$in_file_handle>) {
	chomp;
	
	# header lines
	next if /^##/;
	
	my @split = split /\t/;
	my $data;
		
	# set some variables we'll need later
	$data->{source_id} = $source_id;
	$data->{pop_id} = $pop_id;
	$data->{first_sample_col} = $first_sample_col;
	$data->{sample_ids} = $sample_ids;
	
	# column definition line
	if(/^#/) {
		debug("Parsing header line");
		
		$headers{$split[$_]} = $_ for(0..$#split);
		
		# do sample stuff if required
		if($tables->{sample}) {
			
			# set location of first sample col
			if(defined($headers{FORMAT})) {
				$first_sample_col = $headers{FORMAT} + 1;
				$data->{first_sample_col} = $first_sample_col;
				
				# delete sample IDs
				undef $data->{sample_ids};
				
				# populate sample tables and get sample_ids
				&sample_tables($dbVar, $data, \@split);
				
				$sample_ids = $data->{sample_ids}
			}
			
			# if no sample data
			else {
				delete $tables->{$_} foreach qw(compressed_genotype_region compressed_genotype_var population_genotype);
			}
		}
	}
	
	# data
	else {
		
		$var_counter++;
		$skipped++;
		if($var_counter =~ /000$/) {
			debug "Processed $var_counter lines";
		}
		
		# parse into a hash
		$data->{$_} = $split[$headers{$_}] for keys %headers;
		
		# skip non-variant lines
		next if $data->{ALT} eq '.';
		
		# skip unwanted chromosomes
		next if defined($chrom_regexp) && $data->{'#CHROM'} !~ m/$chrom_regexp/;
		
		## VARIATION
		############
		
		# sometimes ID has many IDs separated by ";", just take the first
		$data->{ID} = (split /\;/, $data->{ID})[0];
		
		# make a var name if none exists
		if($data->{ID} eq '.' || $new_var_name) {
			$data->{ID} =
				($variation_prefix ? $variation_prefix : 'tmp').
				'_'.$data->{'#CHROM'}.'_'.$data->{POS};
		}
		
		# get/set the variation_id
		&variation($dbVar, $data, $tables->{variation}, $check_synonyms) unless $compressed_only;
		
		
		## COORDINATES
		##############
		
		&get_coordinates($data, $use_gp);
		
		
		
		## GET ALLELES
		##############
		
		if(
		   $tables->{variation_feature} ||
		   $tables->{allele} ||
		   $tables->{population_genotype} ||
		   $tables->{individual_genotype_multiple_bp} ||
		   $tables->{compressed_genotype_var} ||
		   $tables->{compressed_genotype_region}
		) {
			&get_alleles($data);
		}
		
		
		
		## VARIATION_FEATURE
		####################
		
		if($tables->{variation_feature} || $merge_vfs) {
			
			# merge variation_features if it's not an indel
			if($merge_vfs && !$data->{is_indel}) {
				&merge_variation_features($dbVar, $data);
			}
			
			# if var doesn't exist or merge failed
			if(
			   !$only_existing &&
			   defined($data->{seq_region}) &&
			   defined($data->{start}) &&
			   !$data->{variation_already_exists} &&
			   !$data->{merged}
			) {
				&variation_feature($dbVar, $data);
			}
		}
		
		# still no var ID
		&get_variation_id_from_variation_feature($dbVar, $data) if(!defined($data->{var_id}));
		
		# skip if no var ID found
		next if !defined($data->{var_id}) && !$force_no_var;
		
		
		## FLANKING_SEQUENCE
		####################
		
		if($tables->{flanking_sequence} && !defined($data->{variation_already_exists}) && !$data->{merged} && $data->{seq_region} && $data->{start}) {
			
			&flanking_sequence($dbVar, $data, $flank_size);
		}
		
		
		
		## PARSE INFO COLUMN
		####################
		
		if(
		   $tables->{allele} ||
		   $tables->{population_genotype} ||
		   $tables->{compressed_genotype_region} ||
		   $tables->{compressed_genotype_var} ||
		   $tables->{individual_genotype_multiple_bp}
		) {
			
			# gets frequencies and genotypes
			&parse_info($data, \@split);
		}
		
		
		
		## ALLELE TABLE
		###############
		
		if($tables->{allele}) {
			
			if($merge_alleles) {
				&merge_alleles($dbVar, $data);
			}
			
			else {
				&allele($dbVar, $data);
			}
		}
		
		
		
		## POPULATION_GENOTYPE
		######################
		
		if(scalar keys %{$data->{gt_freqs}} && $tables->{population_genotype}) {
			
			&population_genotype($dbVar, $data);
		}
		
		
		
		## GENOTYPES
		############
		
		# compressed by var
		if($tables->{compressed_genotype_var}) {
			&genotype_by_var($dbVar, $data);
		}
		
		# multi bp
		if($tables->{individual_genotype_multiple_bp} && $force_multi && @{$data->{genotypes}}) {
			&multi_bp_genotype($dbVar, $data);
		}
		
		# compressed by region
		if($tables->{compressed_genotype_region} && @{$data->{genotypes}}) {
			next unless defined($data->{seq_region}) && defined($data->{start});
			
			for my $i($data->{first_sample_col}..$#split) {
				my $sample_id = $data->{sample_ids}->[$i-$data->{first_sample_col}];
				
				my $gt = $data->{genotypes}->[$i-$data->{first_sample_col}];
				
				# otherwise add to compress hash for writing later
				if(defined($gt) && $gt !~ /\./) {
					
					if (!defined $genotypes->{$sample_id}->{region_start}){
						$genotypes->{$sample_id}->{region_start} = $data->{start};
						$genotypes->{$sample_id}->{region_end} = $data->{end};
					}
					
					# write previous data?
					#compare with the beginning of the region if it is within the DISTANCE of compression
					if (
						defined($genotypes->{$sample_id}->{genotypes}) &&
						(
							(abs($genotypes->{$sample_id}->{region_start} - $data->{start}) > DISTANCE()) ||
							(abs($data->{start} - $genotypes->{$sample_id}->{region_end}) > MAX_SHORT) ||
							(defined($prev_seq_region) && $data->{seq_region} != $prev_seq_region) ||
							($data->{start} - $genotypes->{$sample_id}->{region_end} - 1 < 0)
						)
					) {
						#snp outside the region, print the region for the sample we have already visited and start a new one
						print_file("$TMP_DIR/compressed_genotype_".$$.".txt",$genotypes, $prev_seq_region, $sample_id);
						delete $genotypes->{$sample_id}; #and remove the printed entry
						$genotypes->{$sample_id}->{region_start} = $data->{start};
					}
					
					if ($data->{start} != $genotypes->{$sample_id}->{region_start}){
						#compress information
						my $blob = pack ("w",$data->{start} - $genotypes->{$sample_id}->{region_end} - 1);
						$genotypes->{$sample_id}->{genotypes} .= escape($blob) .escape(pack("w", $data->{var_id} || 0)). escape(pack("w", gt_code($dbVar, $gt)));
					}
					else{
						#first genotype starts in the region_start, not necessary the number
						$genotypes->{$sample_id}->{genotypes} = escape(pack("w", $data->{var_id} || 0)).escape(pack("w", gt_code($dbVar, $gt)));
					}
					
					$genotypes->{$sample_id}->{region_end} = $data->{start};
				}
			}
		}
		
		$prev_seq_region = $data->{seq_region};
		$skipped--;
	}
}

# clean up remaining genotypes and import
if($tables->{compressed_genotype_region}) {

    debug("Importing compressed genotype data");
	
	print_file("$TMP_DIR/compressed_genotype_".$$.".txt",$genotypes, $prev_seq_region);
	&import_genotypes($dbVar);
}

# re-enable keys if requested
if($disable_keys) {
	debug("Re-enabling keys");
	
	foreach my $table(grep {$tables->{$_}} keys %$tables) {
		$dbVar->do(qq{ALTER TABLE $table ENABLE KEYS;})
	}
}

debug("Skipped $skipped variations in the file\n");
debug("Took ", time() - $start_time, "s to run\n");







## SUB-ROUTINES
###############

# prints usage message
sub usage {
	my $usage =<<END;
Usage:
perl import_vcf.pl [arguments]

Options
-h | --help           Display this message and quit


-i | --input_file     Input file - if not specified, attempts to read from STDIN
--tmpdir              Temporary directory to write genotype dump file. Required if
                      writing to compressed_genotype_region [default: no default]
--tmpfile             Name for temporary file [default: compress.txt]

--species             Species to use [default: "human"]
--source              Name of source [required]
--population          Name of population [required]
--ind_prefix          Prefix added to sample names [default: not used]
--var_prefix          Prefix added to constructed variation names [default: not used]
--create_name         Always create a new variation name i.e. don't use ID column
                      [default: not used]
--chrom_regexp        Limit processing to CHROM columns matching regexp
                      [default: not used]

-f | --flank          Size of flanking sequence [default: 200]
--gp                  Use GP tag from INFO column to get coords [default: not used]

--tables              Comma-separated list of tables to include when writing to DB
                      [default: all tables included]
--skip_tables         Comma-separated list of tables to exclude when writing to DB.
                      Takes precedence over --tables (i.e. any tables named in --tables
                      and --skip_tables will be skipped) [default: not used]

--check_synonyms      When looking up variants, also check in variation_synonym.name
                      [default: not used]
--merge_vfs           Attempt to merge VCF variants with existing variation features.
                      Default behaviour is to create a new variation feature entry.
                      [default: not used]
--only_existing       Only write to tables when an existing variant is found. Existing
                      can be a variation with the same name, or from a successful merge
                      using --merge_vfs [default: not used]

--disable_keys        Disable MySQL keys during inserts [default: not used]

-d | --db_host        Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user           Database username [default: "anonymous"]
--password            Database password [default: not used]
-r | --registry       Registry file to use defines DB connections [default: not used]
                      Defining a registry file overrides above connection settings.
END

	print $usage;
}



# gets time
sub getTime() {
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
	my $time = getTime;
	
	print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
} 



# dumps compressed data from hash to temporary file
sub print_file{
    my $file = shift;
    my $genotypes = shift;
    my $seq_region_id = shift;
    my $sample_id = shift;

    open( FH, ">>$file") or die "Could not add compressed information: $!\n";
    if (!defined $sample_id){
		#new chromosome, print all the genotypes and flush the hash
		foreach my $sample_id (keys %{$genotypes}){
			print FH join("\t",
				$sample_id,
				$seq_region_id,
				$genotypes->{$sample_id}->{region_start},
				$genotypes->{$sample_id}->{region_end},
				1,
				$genotypes->{$sample_id}->{genotypes}) . "\n";
		}
    }
    else{
		#only print the region corresponding to sample_id
		print FH join("\t",
			$sample_id,
			$seq_region_id,
			$genotypes->{$sample_id}->{region_start},
			$genotypes->{$sample_id}->{region_end},
			1,
			$genotypes->{$sample_id}->{genotypes}) . "\n";
    }
    close FH;
}



# $special_characters_escaped = printable( $source_string );
sub escape ($) {
	local $_ = ( defined $_[0] ? $_[0] : '' );
	s/([\r\n\t\\\"])/\\$Printable{$1}/sg;
	return $_;
}



# imports genotypes from tmp file to compressed_genotype_region
sub import_genotypes{
    my $dbVar = shift;
    my $call = "mv $TMP_DIR/compressed_genotype_".$$.".txt $TMP_DIR/$TMP_FILE";
    system($call);
    load($dbVar,qw(compressed_genotype_region sample_id seq_region_id seq_region_start seq_region_end seq_region_strand genotypes));
}



# gets seq_region_id to chromosome mapping from DB
sub get_seq_region_ids{
	my $dbVar = shift;
	
	my ($seq_region_id, $chr_name, %seq_region_ids);
	my $sth = $dbVar->prepare(qq{SELECT seq_region_id, name FROM seq_region});
	$sth->execute;
	$sth->bind_columns(\$seq_region_id, \$chr_name);
	$seq_region_ids{$chr_name} = $seq_region_id while $sth->fetch;
	$sth->finish;
	
	return \%seq_region_ids;
}



# gets source_id - retrieves if name already exists, otherwise inserts
sub get_source_id{
	my $dbVar  = shift;
	my $source = shift;
	my $desc   = shift;
	
	my $source_id;
	
	# check existing
	my $sth = $dbVar->prepare(qq{select source_id from source where name = ?});
	$sth->execute($source);
	$sth->bind_columns(\$source_id);
	$sth->fetch;
	$sth->finish;
	
	if(!defined($source_id)) {
		$sth = $dbVar->prepare(qq{insert into source(name, description) values(?,?)});
		$sth->execute($source, $desc);
		$sth->finish();
		$source_id = $dbVar->last_insert_id(undef, undef, qw(source source_id));
	}
	
	return $source_id;
}



# gets population_id - retrieves if already exists, otherwise inserts
sub get_population_id{
	my $dbVar = shift;
	my $population = shift;
	
	my $pop_id;
	
	# check existing
	my $sth = $dbVar->prepare(qq{select sample_id from sample where name = ?});
	$sth->execute($population);
	$sth->bind_columns(\$pop_id);
	$sth->fetch;
	$sth->finish;
	
	if(!defined($pop_id)) {
		# insert into sample
		$sth = $dbVar->prepare(qq{insert into sample(name) values(?)});
		$sth->execute($population);
		$sth->finish;
		$pop_id = $dbVar->last_insert_id(undef, undef, qw(sample sample_id));
		
		# insert into population
		$sth = $dbVar->prepare(qq{insert ignore into population(sample_id) values(?)});
		$sth->execute($pop_id);
		$sth->finish;
	}
	
	return $pop_id;
}



# inserts sample names if they don't exist and gets sample_ids
sub sample_tables {
	my $dbVar = shift;
	my $data = shift;
	my $split_ref = shift;
	
	my @split = @$split_ref;
	
	# populate sample-type tables
	my $sth = $dbVar->prepare(qq{INSERT INTO sample(name) VALUES(?)});
	my $sth2 = $dbVar->prepare(qq{INSERT INTO individual_population(population_sample_id, individual_sample_id) VALUES(?,?)});
	my $sth3 = $dbVar->prepare(qq{INSERT INTO individual(sample_id, individual_type_id) VALUES(?,?)});
	my $sth4 = $dbVar->prepare(qq{select sample_id from sample where name = ?});
	
	for my $i($data->{first_sample_col}..$#split) {
		
		my $sample_name = $sample_prefix.$split[$i];
		my $sample_id;
		$sth4->execute($sample_name);
		$sth4->bind_columns(\$sample_id);
		$sth4->fetch;
		
		if(!$sample_id) {
			$sth->execute($sample_name);
			$sample_id = $dbVar->last_insert_id(undef, undef, qw(sample sample_id));
			$sth2->execute($data->{pop_id}, $sample_id);
			$sth3->execute($sample_id, 3);
		}
		
		push @{$data->{sample_ids}}, $sample_id;
	}
}



# gets variation_id - retrieves if already exists, otherwise inserts
sub variation {
	my $dbVar = shift;
	my $data = shift;
	my $add_new_variations = shift;
	my $check_synonyms = shift;
	
	# check if variation exists
	my $sth = $dbVar->prepare(qq{SELECT variation_id FROM variation WHERE name = ?;});
	$sth->execute($data->{ID});
	$sth->bind_columns(\$data->{var_id});
	$sth->fetch;
	$sth->finish;
	
	if(!defined($data->{var_id}) && $check_synonyms) {
		$sth = $dbVar->prepare(qq{SELECT variation_id FROM variation_synonym WHERE name = ?;});
		$sth->execute($data->{ID});
		$sth->bind_columns(\$data->{var_id});
		$sth->fetch;
		$sth->finish;
	}
	
	# exists
	if(defined($data->{var_id})) {
		$data->{variation_already_exists} = 1;
	}
	
	# does not exist
	elsif($add_new_variations) {
		$sth = $dbVar->prepare(qq{INSERT INTO variation(source_id, name) VALUES (?,?);});
		$sth->execute($data->{source_id}, $data->{ID});
		$sth->finish;
		$data->{var_id} = $dbVar->last_insert_id(undef, undef, qw(variation variation_id));
	}
}



# gets coordinates from CHROM, POS or GP INFO fields
sub get_coordinates {
	my $data = shift;
	my $use_gp = shift;
	
	my $skip;
	
	if($use_gp) {
		foreach my $pair(split /\;/, $data->{INFO}) {
			my ($key, $value) = split /\=/, $pair;
			if($key eq 'GP') {
				($data->{chr}, $data->{start}) = split /\:/, $value;
				$data->{seq_region} = $seq_region_ids->{$data->{chr}};
				$data->{end} = $data->{start};
			}
		}
		
		unless(defined($data->{seq_region}) and defined($data->{start})) {
			warn "Could not determine coordinates from GP INFO field for ", $data->{ID};
			$skip = 1;
		}
	}
	
	else {
		($data->{start}, $data->{end}) = ($data->{POS}, $data->{POS});
		$data->{seq_region} = $seq_region_ids->{$data->{'#CHROM'}};
	}
	
	return $skip;
}



# gets alleles
sub get_alleles {
	my $data = shift;
	
	# work out if this is an indel or multi-bp		
	$data->{is_multi} = 1 if length($data->{REF}) > 1;
	foreach my $alt(split /\,/, $data->{ALT}) {
		$data->{is_multi} = 1 if length($alt) > 1;
	}
	
	# adjust end coord
	$data->{end} += (length($data->{REF}) - 1);
	
	# find out if any of the alt alleles make this an insertion or a deletion
	foreach my $alt_allele(split /\,/, $data->{ALT}) {
		$data->{is_indel} = 1 if $alt_allele =~ /D|I/;
		$data->{is_indel} = 1 if length($alt_allele) != length($data->{REF});
	}
	
	# multiple alt alleles?
	if($data->{ALT} =~ /\,/) {
		if($data->{is_indel}) {
			
			my @alts;
			
			if($data->{ALT} =~ /D|I/) {
				foreach my $alt_allele(split /\,/, $data->{ALT}) {
					# deletion (VCF <4)
					if($alt_allele =~ /D/) {
						push @alts, '-';
					}
					
					elsif($alt_allele =~ /I/) {
						$alt_allele =~ s/^I//g;
						push @alts, $alt_allele;
					}
				}
			}
			
			else {
				$data->{REF} = substr($data->{REF}, 1);
				$data->{REF} = '-' if $data->{REF} eq '';
				$data->{start}++;
				
				foreach my $alt_allele(split /\,/, $data->{ALT}) {
					$alt_allele = substr($alt_allele, 1);
					$alt_allele = '-' if $alt_allele eq '';
					push @alts, $alt_allele;
				}
			}
			
			$data->{ALT} = join "/", @alts;
		}
		
		else {
			# for substitutions we just need to replace ',' with '/' in $alt
			$data->{ALT} =~ s/\,/\//;
		}
	}
	
	else {
		if($data->{is_indel}) {
			# deletion (VCF <4)
			if($data->{ALT} =~ /D/) {
				my $num_deleted = $data->{ALT};
				$num_deleted =~ s/\D+//g;
				$data->{end} += $num_deleted - 1;
				$data->{ALT} = "-";
				$data->{REF} .= ("N" x ($num_deleted - 1)) unless length($data->{REF}) > 1;
			}
			
			# insertion (VCF <4)
			elsif($data->{ALT} =~ /I/) {
				$data->{REF} = '-';
				$data->{ALT} =~ s/^I//g;
				$data->{start}++;
			}
			
			# insertion or deletion (VCF 4+)
			else {
				# chop off first base if they match
				if(substr($data->{REF}, 0, 1) eq substr($data->{ALT}, 0, 1)) {
					$data->{REF} = substr($data->{REF}, 1);
					$data->{ALT} = substr($data->{ALT}, 1);
					
					$data->{start}++;
				}
				
				if($data->{REF} eq '') {
					# make ref '-' if no ref allele left
					$data->{REF} = '-';
				}
				
				# make alt '-' if no alt allele left
				$data->{ALT} = '-' if $data->{ALT} eq '';
			}
		}
	}
	
	@{$data->{alleles}} = ($data->{REF}, split /\,|\//, $data->{ALT});
}



# populates flanking_sequence
sub flanking_sequence {
	my $dbVar = shift;
	my $data = shift;
	my $flank_size = shift;
	
	die "ERROR: Attempting to insert into flanking_sequence with missing data" unless
		defined($data->{var_id}) &&
		defined($data->{start}) &&
		defined($data->{end}) &&
		defined($data->{seq_region});
	
	my $sth = $dbVar->prepare(qq{
		INSERT INTO flanking_sequence(
			variation_id,
			up_seq_region_start,
			up_seq_region_end,
			down_seq_region_start,
			down_seq_region_end,
			seq_region_id,
			seq_region_strand
		)
		VALUES (?,?,?,?,?,?,?);
	});
	
	$sth->execute(
		$data->{var_id},
		$data->{start} - $flank_size,
		$data->{start} - 1,
		$data->{end} + 1,
		$data->{end} + $flank_size,
		$data->{seq_region},
		1
	);
	$sth->finish();
}



# populates variation_feature
sub variation_feature {
	my $dbVar = shift;
	my $data = shift;
	
	die "ERROR: Attempting to insert into variation_feature with missing data" unless
		defined($data->{var_id}) &&
		defined($data->{start}) &&
		defined($data->{end}) &&
		defined($data->{alleles}) &&
		scalar @{$data->{alleles}} &&
		defined($data->{source_id}) &&
		defined($data->{seq_region});
	
	my $sth = $dbVar->prepare(qq{
		INSERT INTO variation_feature(
			variation_id,
			seq_region_id,
			seq_region_start,
			seq_region_end,
			seq_region_strand,
			variation_name,
			allele_string,
			map_weight,
			source_id
		)
		VALUES (?,?,?,?,?,?,?,?,?);
	});
	
	$sth->execute(
		$data->{var_id},
		$data->{seq_region},
		$data->{start},
		$data->{end},
		1,
		$data->{ID},
		(join "/", @{$data->{alleles}}),
		1,
		$data->{source_id} #needs changing
	);
	$sth->finish;
}



# merges variation feature with existing
sub merge_variation_features {
	my $dbVar = shift;
	my $data = shift;
	
	die "ERROR: Attempting to merge variation_features with missing data" unless
		defined($data->{start}) &&
		defined($data->{end}) &&
		defined($data->{alleles}) &&
		defined($data->{seq_region});
		
	my $sth = $dbVar->prepare(qq{
		SELECT variation_feature_id, variation_id, allele_string, variation_name
		FROM variation_feature
		WHERE seq_region_id = ?
		AND seq_region_start = ?
		AND seq_region_end = ?
		AND seq_region_strand = 1
		AND map_weight = 1
	});
	
	my ($vf_id, $variation_id, $allele_string, $variation_name, $new_allele_string);
	$sth->execute($data->{seq_region}, $data->{start}, $data->{end});
	$sth->bind_columns(\$vf_id, \$variation_id, \$allele_string, \$variation_name);
	
	my ($row_count);
	
	while($sth->fetch) {
		next if $allele_string =~ /\-/;
		
		$row_count++;
		
		my @existing_alleles = split /\//, $allele_string;
		
		# compare ref alleles - we don't want to merge if they differ
		return unless $existing_alleles[0] eq $data->{alleles}->[0];
		
		my %new_alleles = ();
		$new_alleles{$_}++ for @existing_alleles;
		$new_alleles{$_}++ for @{$data->{alleles}};
		
		if(scalar keys %new_alleles != scalar @existing_alleles) {
			$new_allele_string =
				$allele_string.
				'/'.
				(join /\//, grep {$new_alleles{$_} == 1} @{$data->{alleles}});
		}
	}
	$sth->finish;
	
	# merge OK
	if($row_count == 1) {
		
		$data->{merged} = 1;
		
		# delete original entry from variation
		unless($data->{variation_already_exists}) {
			$sth = $dbVar->prepare(qq{DELETE FROM variation WHERE variation_id = ?});
			$sth->execute($data->{var_id});
			$sth->finish;
		}
		
		# update the variation_id
		$data->{var_id} = $variation_id;
		
		# new allele string
		if(defined($new_allele_string)) {
			$sth = $dbVar->prepare(qq{
				UPDATE variation_feature
				SET allele_string = ?
				WHERE variation_feature_id = ?
			});
			$sth->execute($new_allele_string, $vf_id);
			$sth->finish;
		}
		
		# add entry to variation_synonym
		if ($variation_name ne $data->{ID}) {
			$sth = $dbVar->prepare(qq{
				INSERT IGNORE INTO variation_synonym(
					variation_id,
					source_id,
					name
				)
				VALUES(?, ?, ?)
			});
			
			$sth->execute($data->{var_id}, $data->{source_id}, $data->{ID});
			$sth->finish;
		}
	}
}


# gets variation ID from variation_feature overlap
sub get_variation_id_from_variation_feature {
	my $dbVar = shift;
	my $data = shift;
	
	unless(
		defined($data->{start}) &&
		defined($data->{end}) &&
		defined($data->{seq_region})
	) {
		warn "WARNING: Attempting to retrieve variation ID with missing data";
		return;
	}
		
	my $sth = $dbVar->prepare(qq{
		SELECT variation_id, allele_string, seq_region_strand
		FROM variation_feature
		WHERE seq_region_id = ?
		AND seq_region_start = ?
		AND seq_region_end = ?
	});
	
	$sth->execute($data->{seq_region}, $data->{start}, $data->{end});
	
	my ($var_id, $allele_string, $strand);
	$sth->bind_columns(\$var_id, \$allele_string, \$strand);
	
	while($sth->fetch) {
		my %vf_alleles = map {$_ => 1} split /\//, $allele_string;
		
		if($strand == -1) {
			my %new_alleles;
			
			foreach(keys %vf_alleles) {
				reverse_comp(\$_);
				$new_alleles{$_} = 1;
			}
			
			%vf_alleles = %new_alleles;
		}
		
		next unless grep {$vf_alleles{$_}} grep {$_ ne 'N'} @{$data->{alleles}};
		
		$data->{var_id} = $var_id;
		
		last;
	}
}

# parses info column and gets frequencies
sub parse_info {
	my $data = shift;
	my $split_ref = shift;
	
	my @split = @$split_ref;
	
	# parse info column
	my %info;	
	foreach my $chunk(split /\;/, $data->{INFO}) {
		my ($key, $val) = split /\=/, $chunk;
		$info{$key} = $val;
	}
	
	# get GMAF if available
	$data->{GMAF} = $info{GMAF} if defined $info{GMAF};
	
	# deal with frequencies
	my (@freqs, %gt_freqs, @genotypes);
	
	# may be defined in INFO column
	if(defined($info{AF})) {
		@freqs = split /\,/, $info{AF};
		my $total_alt_freq = 0;
		$total_alt_freq += $_ for @freqs;
		unshift @freqs, 1 - $total_alt_freq;
	}
	
	my $total_count;
	my %counts;
	my %allele_counts;
	
	# get counts
	if(defined($data->{first_sample_col})) {
		for my $i($data->{first_sample_col}..$#split) {
			my $gt;
			my @bits;
			my $gt = (split /\:/, $split[$i])[0];
			foreach my $bit(split /\||\/|\\/, $gt) {
				push @bits, ($bit eq '.' ? '.' : $data->{alleles}->[$bit]);
			}
			
			if(scalar @bits) {
				$allele_counts{$_}++ for @bits;
				
				$gt = join '|', sort @bits;
				
				# store genotypes for later
				push @genotypes, $gt;
				
				unless($gt =~ /\./) {
					$counts{$gt}++;
					$total_count++;
				}
			}
		}
		
		# now calculate frequencies
		for my $i(0..(scalar @{$data->{alleles}} - 1)) {
			my ($c_aa, $c_ab);
			my $a = $data->{alleles}->[$i];
			
			$c_aa = $counts{$a.'|'.$a};
			
			foreach my $gt(keys %counts) {
				$c_ab += $counts{$gt} if $gt =~ /$a/ and $gt ne $a.'|'.$a;
				$gt_freqs{$gt} = ($total_count ? $counts{$gt}/$total_count : 0);
			}
			
			push @freqs, ($total_count ? (((2*$c_aa) + $c_ab) / (2*$total_count)) : 0) unless defined($info{AF});
		}
	}
	
	$data->{freqs} = \@freqs;
	$data->{genotypes} = \@genotypes;
	$data->{gt_freqs} = \%gt_freqs;
	$data->{gt_counts} = \%counts;
	$data->{allele_counts} = \%allele_counts;
}



# populates allele table
sub allele {
	my $dbVar = shift;
	my $data = shift;
	
	die "ERROR: Attempting to populate allele table with missing data" unless
		defined($data->{var_id}) &&
		defined($data->{alleles}) &&
		defined($data->{freqs}) &&
		defined($data->{pop_id});
	
	my $sth = $dbVar->prepare(qq{
		INSERT IGNORE INTO allele(variation_id, subsnp_id, allele_code_id, frequency, sample_id, count)
		VALUES (?,NULL,?,?,?,?)
	});
	
	$sth->execute(
		$data->{var_id},
		allele_code($dbVar, $data->{alleles}->[$_]),
		$data->{freqs}->[$_],
		$data->{pop_id},
		$data->{allele_counts}->{$data->{alleles}->[$_]} || 0
	) for (0..(scalar @{$data->{alleles}} - 1));
	
	$sth->finish;
}



# attempts to merge in entries with existing in allele table
sub merge_alleles {
	my $dbVar = shift;
	my $data = shift;
	
	die "ERROR: Attempting to populate allele table with missing data" unless
		defined($data->{var_id}) &&
		defined($data->{alleles}) &&
		defined($data->{freqs}) &&
		defined($data->{pop_id});
		
	my $sth = $dbVar->prepare(qq{
		UPDATE allele
		SET frequency = ?, count = ?
		WHERE variation_id = ?
		AND sample_id = ?
		AND allele_code_id = ?
	});
	
	my $rows_affected;
	
	for my $i(0..(scalar @{$data->{alleles}} - 1)) {
		my $allele = $data->{alleles}->[$i];
		
		$rows_affected = $sth->execute(
			$data->{freqs}->[$i],
			$data->{allele_counts}->{$allele} || 0,
			$data->{var_id},
			$data->{pop_id},
			allele_code($allele),
		);
	}
	
	$sth->finish;
}



# populates individual_genotype_multiple_bp table
sub multi_bp_genotype {
	my $dbVar = shift;
	my $data = shift;
	
	die "ERROR: Attempting to populate individual_genotype_multiple_bp table with missing data" unless
		defined($data->{var_id}) &&
		defined($data->{genotypes}) &&
		defined($data->{sample_ids});
	
	if(scalar @{$data->{genotypes}}) {
		
		my @multi_rows;
		
		for my $i(0..(scalar @{$data->{sample_ids}}) - 1) {
			my $sample_id = $data->{sample_ids}->[$i];
			my ($a1, $a2) = split /\|/, $data->{genotypes}->[$i];
			next if $a1 eq '.' || $a2 eq '.';
			
			push @multi_rows,
				"(".
					(join ",",
						(
							$data->{var_id},
							"'$a1'",
							"'$a2'",
							$sample_id
						)
					).
				")";
		}
		
		my $vals = join ",", @multi_rows;
		my $sth = $dbVar->prepare(qq{
			INSERT IGNORE INTO individual_genotype_multiple_bp(
				variation_id,
				allele_1,
				allele_2,
				sample_id
			)
			values$vals
		});
		$sth->execute;
		$sth->finish;
	}
}


# genotype by variation
sub genotype_by_var {
	my $dbVar = shift;
	my $data = shift;
	
	die "ERROR: Attempting to populate compressed_genotype_var table with missing data" unless
		defined($data->{var_id}) &&
		defined($data->{genotypes}) &&
		defined($data->{sample_ids});

	if(scalar @{$data->{genotypes}}) {
		
		my $genotype_string = '';
		
		for my $i(0..(scalar @{$data->{sample_ids}}) - 1) {
			my $sample_id = $data->{sample_ids}->[$i];
			my $gt = $data->{genotypes}->[$i];
			next if $gt =~ /\./;
			
			$genotype_string .= pack("ww", $sample_id, &gt_code($dbVar, $gt));
		}
		if($genotype_string ne '') {
			my $sth = $dbVar->prepare(qq{
				INSERT INTO compressed_genotype_var(
					variation_id,
					subsnp_id,
					genotypes
				)
				VALUES(?, NULL, ?)
			});
			$sth->execute($data->{var_id}, $genotype_string);
			$sth->finish;
		}
	}
}

# gets/sets genotype codes
sub gt_code {
	my $dbVar = shift;
	my $gt = shift;
	
	if(!defined($gt_codes{$gt})) {
		my $new_code = (sort {$a <=> $b} values %gt_codes)[-1] + 1;
		create_genotype_code($dbVar, $gt, $new_code);
		$gt_codes{$gt} = $new_code;
	}
	
	return $gt_codes{$gt};
}

# creates a new genotype_code entry
sub create_genotype_code {
	my $dbVar = shift;
	my $genotype = shift;
	my $gt_code = shift;
	my $hap_id = 1;
	
	my $sth = $dbVar->prepare("INSERT INTO genotype_code VALUES(?, ?, ?)");
	
	foreach my $allele(split /\|/, $genotype) {
		my $allele_code = allele_code($dbVar, $allele);
		$sth->execute($gt_code, $allele_code, $hap_id++);
	}
	
	$sth->finish;
	return $gt_code;
}

# getter/setter for allele code
sub allele_code {
	my $dbVar = shift;
	my $allele = shift;
	
	# create if doesn't exist
	if(!defined($allele_codes{$allele})) {
		my $sth = $dbVar->prepare("INSERT INTO allele_code(allele) VALUES(?)");
		$sth->execute($allele);
		$sth->finish();
		
		$allele_codes{$allele} = $dbVar->last_insert_id(undef, undef, qw(allele_code allele_code_id));
	}
	
	return $allele_codes{$allele};
}


# populates population_genotype table
sub population_genotype {
	my $dbVar = shift;
	my $data = shift;
	
	die "ERROR: Attempting to populate population_genotype table with missing data" unless
		defined($data->{var_id}) &&
		defined($data->{gt_freqs}) &&
		defined($data->{pop_id});
	
	
	my $sth = $dbVar->prepare(qq{
		INSERT IGNORE INTO population_genotype(
			variation_id,
			subsnp_id,
			genotype_code_id,
			frequency,
			sample_id,
			count
		)
		VALUES (?,NULL,?,?,?,?)
	});
	
	foreach my $gt(keys %{$data->{gt_freqs}}) {
		$sth->execute(
			$data->{var_id},
			gt_code($dbVar, $gt),
			$data->{gt_freqs}->{$gt},
			$data->{pop_id},
			$data->{gt_counts}->{$gt}
		);
	}
	
	$sth->finish;
}